module CRThreeBodyProblem

using LinearAlgebra

"""
Initial value problem with a system of autonomous differential equations (no explicit time dependence).
"""
struct ADE_IVProblem
    time_int # time interval [a, b]
    y0 # initial state
    f # state derivative function y' = f(y)
end

"""
    sol = IVP_solution(t, y, problem)

Solution of an initial value problem. Contains the time steps `t`, the state values `y` and the original problem.
"""
struct IVP_solution
    t # times at which the solution was computed
    y # state approximations y[i] ≈ y(t[i])
    problem::ADE_IVProblem
end

"""
    params = ParamsRK4(n)

Parameters for explicit 4th order 4-stage Runge-Kutta integrator.
- `n`: number of steps
"""
struct ParamsRK4
    n::Integer # number of steps
end

"""
    sol = solve(prob::ADE_IVProblem, params::ParamsRK4)

Solve the initial value problem `prob` with the 4th order 4-stage Runge-Kutta method.
"""
function solve(prob::ADE_IVProblem, params::ParamsRK4)
    n = params.n
    y = zeros(n+1, length(prob.y0))
    y[1,:] = prob.y0
    x = LinRange(prob.time_int..., n+1)
    h = x[2] - x[1]
    for i=1:n
        k1 = h * prob.f(y[i, :])
        k2 = h * prob.f(y[i, :] + k1/2)
        k3 = h * prob.f(y[i, :] + k2/2)
        k4 = h * prob.f(y[i, :] + k3)
        y[i+1, :] = y[i, :] + (k1 + 2*k2 + 2*k3 + k4) / 6
    end
    return IVP_solution(x, y, prob)
end


"""
    params = ParamsDOPRI5(ε, h0, η)

Parameters for the DOPRI5 integrator.
- `ε`: maximum local error on unit interval (usually good approximate of max global error)
- `h0`: initial step size
- `η`: step damping factor (0 < η < 1), default: 0.9
"""
@kwdef struct ParamsDOPRI5
    ε # max local error on unit interval
    h0 # initial step size (NOTE: could try b - a?)
    η = 0.9 # step dampening factor
end

"""
    sol = solve(prob::ADE_IVProblem, params::ParamsDOPRI5)

Solve the initial value problem `prob` with the 5th order Dormand-Prince method.
"""
function solve(problem::ADE_IVProblem, par::ParamsDOPRI5)
    h = par.h0
    t, b = problem.time_int
    y = problem.y0
    times = [t]
    ys = [y]
    f = problem.f

    while t < b
        k1 = h*f(y)
        k2 = h*f(y + k1/5)
        k3 = h*f(y + k1*0.075 + k2*0.225)
        k4 = h*f(y+44*k1/45-56*k2/15+32*k3/9)
        k5 = h*f(y+19372*k1/6561-25360*k2/2187+64448*k3/6561-212*k4/729)
        k6 = h*f(y+9017*k1/3168-355*k2/33+46732*k3/5247+49*k4/176-5103*k5/18656)
        k7 = h*f(y+35*k1/384+500*k3/1113+125*k4/192-2187*k5/6784+11*k6/84)
        le = 71*k1/57600-71*k3/16695+71*k4/1920-17253*k5/339200+22*k6/525-k7/40;
        nle = norm(le); # NOTE: does this norm generalize well to higher dimensions?
        
        if nle < par.ε * h
            # accept step
            y = y + (5179*k1/57600+7571*k3/16695+393*k4/640-92097*k5/339200+187*k6/2100+k7/40);
            t += h
            push!(times, t)
            push!(ys, y)

            # compute new step size
            h = par.η * h*(par.ε*h/nle)^(1/5)
            if t + h > b
                h = b - t
            end
        else
            # reduce step
            h /= 2
        end
    end
    IVP_solution(times, hcat(ys...)', problem) 
end


"""
    problem = CR3BP(m1, m2)

Circular Restricted Three-Body Problem.
"""
struct CR3BP
    rel_mass::Float64 # relative mass of the lighter body
    pos1::Vector{Float64} # position of the primary (heavier) body
    pos2::Vector{Float64} # position of the secondary (lighter) body

    """
        problem = CR3BP(m1, m2)

    Create a new Circular Restricted Three-Body Problem with given masses of the primary and secondary bodies.
    - m1, m2: masses of the two bodies (default: Earth and Moon masses in kg)
    """
    function CR3BP(m1::Float64 = 5.972e24, m2::Float64 = 7.342e22)
        μ = min(m1, m2) / (m1 + m2)
        return new(
            μ,
            [-μ, 0.0, 0.0],
            [1-μ, 0.0, 0.0]
        )
    end
end


"""
    sol = solve_CR3BP(problem, init_pos, init_vel, duration, integrator_params)

Solve the Circular Restricted Three-Body Problem for t ∈ [0, `duration`] using non-dimensional equations of motion.
States `sol.y[i]` are 6D vectors [x, y, z, vx, vy, vz] with the position and velocity of the tertiary body.
Parameters:
- `problem`: CR3BP problem
- `init_pos`: initial position of the tertiary body.
              Use `problem.pos1` and `problem.pos2` to set initial position
              relative to the primary or secondary bodies.
              Note that one unit of distance is the distance between the two primary bodies.
- `init_vel`: initial velocity of the tertiary body.
- `duration`: duration of the simulation.
              Note that one unit of time is the secondary masses (e.g. moon's) orbital period divided by 2π.
- `integrator_params`: may be either of type `ParamsRK4` or `ParamsDOPRI5`
                       (this will determine the integrator used).
"""
function solve_CR3BP(problem::CR3BP; init_pos::Vector{Float64}, init_vel::Vector{Float64}, duration::Float64,
                     integrator_params)
    
    if length(init_pos) != 3 || length(init_vel) != 3
        throw(ArgumentError("Initial position and velocity must be 3D vectors"))
    end
    μ = problem.rel_mass
    μ1 = 1 - μ

    function state_deriv(state)
        x, y, z, vx, vy, vz = state
        y2 = y^2
        z2 = z^2
        rx = x + μ
        Rx = x - μ1
        
        μ1R3 = μ1 / (rx^2 + y2 + z2)^(3/2)
        μr3 = μ / (Rx^2 + y2 + z2)^(3/2)

        ax = x + 2vy - μ1R3*rx - μr3*Rx
        ay = y - 2vx - (μ1R3 + μr3)*y
        az = (-μ1R3 - μr3)*z
        return [vx, vy, vz, ax, ay, az]
    end

    initial_state = [init_pos; init_vel] # [x, y, z, vx, vy, vz] ∈ phase space
    ivp = ADE_IVProblem([0.0, duration], initial_state, state_deriv)
    return solve(ivp, integrator_params)
end

"""
    L4, L5 = stable_Lagrange_points(problem)

Return the (nominally stable) equilateral Lagrange points L4, L5 for the CR3B problem.
"""
function stable_Lagrange_points(problem::CR3BP)
    L4 = [0.5 - problem.rel_mass, sqrt(3)/2, 0.0]
    L5 = [0.5 - problem.rel_mass, -sqrt(3)/2, 0.0]
    return L4, L5
end

export CR3BP, solve_CR3BP, stable_Lagrange_points,
       ADE_IVProblem, IVP_solution, ParamsRK4, ParamsDOPRI5

end # module ThreeBodyProblem
