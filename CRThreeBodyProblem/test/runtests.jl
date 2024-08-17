using Test, CRThreeBodyProblem, LinearAlgebra

@testset "stable Lagrange points, DOPRI5" begin
    earth_moon = CR3BP();
    dopri_param = ParamsDOPRI5(ε = 1e-9, h0 = 0.01)

    for lp in stable_Lagrange_points(earth_moon)
        path = solve_CR3BP(
            earth_moon,
            init_pos = lp,
            init_vel = zeros(3),
            duration = 2π * 100,
            integrator_params = dopri_param
        )
        max_dist = maximum([norm(path.y[i, 1:3] - lp) for i in 1:size(path.y, 1)])
        @test max_dist < dopri_param.ε * 10e2
    end
end

@testset "stable Lagrange points, RK4" begin
    earth_moon = CR3BP();
    rk4_param = ParamsRK4(10_000)

    for lp in stable_Lagrange_points(earth_moon)
        path = solve_CR3BP(
            earth_moon,
            init_pos = lp,
            init_vel = zeros(3),
            duration = 2π * 3,
            integrator_params = rk4_param
        )
        max_dist = maximum([norm(path.y[i, 1:3] - lp) for i in 1:size(path.y, 1)])
        @test max_dist < 1e-3
    end
end

