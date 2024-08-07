using Test, SplineInterpolation


@testset "Two points (dx = 1)" begin
    x = [0, 1]
    y = [1, 3]
    S = interpolate(x, y)
    @test S.x == x
    @test S.a == y
    # should be a straight line with slope 2
    @test S.b ≈ [2]
    @test S.c ≈ [0]
    @test S.d ≈ [0]

    @test S(0) ≈ 1
    @test S(0.5) ≈ 2
    @test S(0.75) ≈ 2.5
    @test S(1) ≈ 3
end

@testset "Two points" begin
    x = [-1, 1.5]
    y = [1.1, 0]
    S = interpolate(x, y)
    @test S.x == x
    @test S.a == y
    @test S.b ≈ [(y[2] - y[1]) / (x[2] - x[1])]
    @test S.c ≈ [0]
    @test S.d ≈ [0]
end

@testset "Two points with equal y" begin
    x = [0, 1]
    y = [1, 1]
    S = interpolate(x, y)
    @test S.b ≈ [0]
    @test S.c ≈ [0]
    @test S.d ≈ [0]
end


@testset "Three points" begin
    x = [0, 1.5, 19]
    y = [3, -2, 0.1]
    S = interpolate(x, y)
    @test S.x == x
    @test S.a == y
    @test length(S.b) == 2
    @test length(S.c) == 2
    @test length(S.d) == 2
end


@testset "Sine function, many samples" begin
    x = range(0, 2π, length=50)
    y = sin.(x)
    S = interpolate(x, y)
    maxerr = maximum(abs.(S.(x) .- y))
    @test maxerr < 1e-5
end