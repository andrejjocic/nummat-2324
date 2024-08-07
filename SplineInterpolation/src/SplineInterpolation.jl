module SplineInterpolation

using SparseArrays, Plots

export ncs_system, interpolate, evaluate, plot_spline

"""
    M, k = ncs_matrix(x, y)

Return sparse matrix M for solving natural cubic spline interpolation,
and the right-hand side vector k. The solution to M*u = k gives the coefficients
u = [b1,c1,d1,b2,c2,d2,...,bn-1,cn-1,dn-1] of the interpolating polynomials.
"""
function ncs_system(x, y)
    n = length(x)
    if n != length(y)
        error("x and y must have the same length")
    end
    dx = diff(x)
    dx2 = dx.*dx
    dx3 = dx2.*dx
    dy = diff(y)

    M = spzeros(3n - 3, 3n - 3)
    k = zeros(3n - 3)

    # continuity of the interpolating function
    for i = 1:n-1 
        base = 3(i-1) + 1
        # b_i dx_i + c_i dx_i^2 + d_i dx_i^3 = dy_i
        M[base, base:base+2] .= [dx[i], dx2[i], dx3[i]]
        k[base] = dy[i]
    end

    # continuity of the first derivative
    for i = 1:n-2
        base = 3(i-1) + 1
        # b_i + 2c_i dx_i + 3d_i dx_i^2 = b_{i+1}
        M[base+1, base:base+3] .= [1, 2dx[i], 3dx2[i], -1] 
    end 

    # inflection point at the start: c_1 = 0
    M[end - 1, 2] = 1 # end - 1 is the unused row (no condition on first deriv at x_n)

    # continuity of the second derivative + inflection point at the end
    for i = 1:n-1 
        base = 3(i-1) + 1
        # c_i + 3d_i dx_i = c_{i+1}
        M[base+2, base+1:base+2] = [1, 3dx[i]]

        if i == n-1 # inflection point condition
            continue # (c_n is not even defined)
        end
        M[base+2, base+4] = -1 
    end

    return M, k
end

"""
Read-only struct holding the coefficients of a cubic spline interpolant, as well as the x values.
The interpolant between x[i] and x[i+1] is given by a[i] + b[i]*dx + c[i]*dx^2 + d[i]*dx^3.
"""
struct CubicSpline
    x::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

"""
    S = interpolate(x, y)

Interpolate a set of points (x,y) with a cubic spline.
Assumes sorted and pairwise distinct x values (x_1 < x_2 < ... < x_n).
"""
function interpolate(x, y)::CubicSpline
    if length(x) != length(y)
        error("x and y must have the same length")
    elseif length(x) < 2
        error("need at least 2 points to interpolate")
    end
    M, k = ncs_system(x, y)
    coeffs = M \ k # NOTE: M could be singular if x_i are not unique
    
    return CubicSpline(
        # no need to copy x,y (read-only struct)
        x, # we need all the xs for range checking at evaluation
        y, # last y is redundant (can be interpolated exactly)
        coeffs[1:3:end],
        coeffs[2:3:end],
        coeffs[3:3:end]
    )
end


"""
    y = evaluate(cs, x)
    
Evaluate the interpolating cubic spline at x, which
must be within the domain of the spline.
"""
function evaluate(cs::CubicSpline, x)
    # find the index i such that x ∈ [cs.x[i], cs.x[i+1]]
    if x < cs.x[1] || x > cs.x[end]
        error("x out of range")
    elseif x == cs.x[1]
        return cs.a[1] # avoid the case i = 0 below
    elseif x == cs.x[end]
        return cs.a[end] # no point evaluating last polynomial
        # NOTE: would it be faster if we skip this check and just evaluate it anyway?
    end

    i = searchsortedfirst(cs.x, x) - 1
    dx = x - cs.x[i]
    return cs.a[i] + cs.b[i]*dx + cs.c[i]*dx^2 + cs.d[i]*dx^3
end

(s::CubicSpline)(x) = evaluate(s, x) # syntactic sugar


"""
Plot the interpolating cubic spline on its domain.

The resolution parameter controls the number of points used to plot each polynomial segment.
"""
function plot_spline(cs::CubicSpline, resolution = 50; kwargs...)
    # scatter(cs.x, cs.a, marker = :circle, color = :gray, markersize = 4, label = "Data Points")
    for i in 1:length(cs.x)-1
        x_range = range(cs.x[i], stop = cs.x[i+1], length = resolution)
        dx = x_range .- cs.x[i]
        dx2 = dx.*dx
        dx3 = dx2.*dx
        y_range = cs.a[i] .+ cs.b[i].*dx .+ cs.c[i].*dx2 .+ cs.d[i].*dx3
        plot!(x_range, y_range, color = (i % 2 == 1 ? :red : :blue), primary = false; kwargs...)
    end
    return plot!()
end

end # module SplineInterpolation