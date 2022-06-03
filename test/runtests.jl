using Quadrotor
using Test

@testset "Quadrotor.jl" begin
    # Test syntax for using the QuadrotorParameters struct
    params = Quadrotor.QuadrotorParameters(1.0, 1.0, 0.2, 0.15, 9.81)
    @test isapprox(params.M, 1.0)
    @test isapprox(params.m, 1.0)
    @test isapprox(params.L, 0.20)
    @test isapprox(params.l, 0.15)
    @test isapprox(params.g, 9.81)
end
