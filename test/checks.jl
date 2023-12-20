using Test
using DormandPrince.DP5Impl:
    Options,
    check_max_allowed_steps,
    check_uround,
    check_beta,
    check_safety_factor

@testset "Test Checks" begin

    @testset "Check Max Allowed Steps" begin
        options = Options{Float64}(maximum_allowed_steps=-1)
        @test check_max_allowed_steps(options) == false

        options = Options{Float64}(maximum_allowed_steps=1)
        @test check_max_allowed_steps(options) == true
    end

    @testset "Check uround" begin
        options = Options{Float64}(uround=1e-36)
        @test check_uround(options) == false

        options = Options{Float64}(uround=1.1)
        @test check_uround(options) == false

        options = Options{Float64}(uround=1e-10)
        @test check_uround(options) == true
    end

    @testset "Check beta" begin
        options = Options{Float64}(beta=0.3)
        @test check_beta(options) == false

        options = Options{Float64}(beta=0.01)
        @test check_beta(options) == true
    end

    @testset "Check safety factor" begin
        options = Options{Float64}(safety_factor=1.1)
        @test check_safety_factor(options) == false

        options = Options{Float64}(safety_factor=1e-5)
        @test check_safety_factor(options) == false

        options = Options{Float64}(safety_factor=0.5)
        @test check_safety_factor(options) == true
    end

end
