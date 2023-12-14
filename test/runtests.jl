using Test
using LinearAlgebra
using DormandPrince: DP5Solver, dopri5, integrate

include("exact_evol_helpers.jl")


function fcn(x, y, f)
    g(x) = 2.2*2π*sin(2π*x)

    f[1] = -1im * g(x)/2 * y[2]
    f[2] = -1im * g(x)/2 * y[1]
end

# standalone solver test
@testset "Integration Test" begin
    solver = DP5Solver(
        fcn,
        0.0,
        ComplexF64[1.0, 0.0]
    )

    dopri5(solver, 2π)

    @test solver.y ≈ solution(2π)
end

# Test integrate() 
@testset "Integration Interface Test" begin

    # test iterator generation
    @testset "Iterator Interface" begin
        times = [0.1, 0.5, 1.1]
        exact_values = []
        dp5_values = []
        
        # get exact values
        for t in times
            push!(exact_values, solution(t))
        end

        # use iterator to get exact values
        solver = DP5Solver(
            fcn,
            0.0,
            ComplexF64[1.0, 0.0]
        )

        iter = integrate(solver, times)

        for (t,y) in iter
            push!(dp5_values, copy(y))
        end

        @test dp5_values ≈ exact_values
    end

    @testset "Callback Interface" begin
        times = [0.1, 0.5, 1.1]
        callback_times = []
        exact_values = []
        dp5_values = []
        
        # get exact values
        for t in times
            push!(exact_values, solution(t))
        end

        # use iterator to get exact values
        solver = DP5Solver(
            fcn,
            0.0,
            ComplexF64[1.0, 0.0]
        )

        integrate(solver, times) do t, y
            push!(callback_times, t)
            push!(dp5_values, copy(y))
        end

        @test dp5_values ≈ exact_values
    end
end


