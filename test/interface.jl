using Test
using LinearAlgebra
using DormandPrince

function evolution_operator(t::Float64)
    ϕ = 2.2 * sin(π * t)^2
    U = zeros(ComplexF64, 2,2)
    U[1,1] =  1 / sqrt(2)
    U[2,1] =  1 / sqrt(2)
    U[2,2] =  1 / sqrt(2)
    U[1,2] = -1 / sqrt(2)

    U * diagm(exp.([-im*ϕ, im*ϕ])) * U'
end

function solution(t)
    U = evolution_operator(t)
    return U * [1.0, 0.0]

end

function fcn(x, y, f)
    g(x) = 2.2*2π*sin(2π*x)

    f[1] = -1im * g(x)/2 * y[2]
    f[2] = -1im * g(x)/2 * y[1]
end

# standalone solver test
@testset "Integration Test" begin

    for SolverType in [DP5Solver, DP8Solver]
        solver = SolverType(
            fcn,
            0.0,
            ComplexF64[1.0, 0.0]
        )

        integrate!(solver, 2π)

        @test get_current_state(solver) ≈ solution(2π)
    end
    
end

# Test integrate() 
@testset "Integration Interface Test" begin

    # test iterator generation
    @testset "Iterator Interface" begin
        times = [0.1, 0.5, 1.1]
        exact_values = []
        
        # get exact values
        for t in times
            push!(exact_values, solution(t))
        end

        # use iterator to get exact values
        
        for SolverType in [DP5Solver, DP8Solver]
            values = []
            solver = SolverType(
                fcn,
                0.0,
                ComplexF64[1.0, 0.0]
            )
            
            iter = SolverIterator(solver, times)

            for (t,y) in iter
                push!(values, copy(y))
            end

            @test values ≈ exact_values
        end


    end

    @testset "Callback Interface" begin
        times = [0.1, 0.5, 1.1]
        exact_values = []

        # get exact values
        for t in times
            push!(exact_values, solution(t))
        end


        for SolverType in [DP5Solver, DP8Solver]
            values = []
            callback_times = []

            # use iterator to get exact values
            solver = SolverType(
                fcn,
                0.0,
                ComplexF64[1.0, 0.0]
            )

            integrate!(solver, times) do t, y
                push!(callback_times, t)
                push!(values, copy(y))
            end

            @test values ≈ exact_values
            @test callback_times == times
        end
    end

    @testset "integrate! exception" begin 
        solver = DP5Solver(
            fcn,
            0.0,
            ComplexF64[1.0, 0.0]
            ;
            maximum_allowed_steps=10

        )

        @test_throws DormandPrince.DPException integrate!(solver, 2π)
    end
end