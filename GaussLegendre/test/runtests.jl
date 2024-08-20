using Test, GaussLegendre

@testset "analytically integrable functions" begin
    @test integrate_GLQ1((x) -> x^2, 0, 1, 100) ≈ 1/3
    @test integrate_GLQ1(sin, 0, π, 100) ≈ 2
end

@testset "Si(x)" begin
    @test isapprox(integrate_GLQ1((x) -> sin(x) / x, 0, 5, 100), 1.54993, atol=1e-5)
end

@testset "polynomials" begin
    # ensure exactness for polynomials of degree <= 3 (any number of sub-intervals)
    p3 = (x) -> 3x^3 - 2x^2 + 5x - 7
    a, b = -1, 0.3
    @test integrate_GLQ1(p3, a, b, 1) ≈ integrate_GLQ1(p3, a, b, 10)
end