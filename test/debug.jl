using LinearAlgebra
using DormandPrince: DP8Solver, DP5Solver, integrate

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


function run()
    solver = DP8Solver(
        fcn,
        0.0,
        ComplexF64[1.0, 0.0]
    )
    
    report = integrate(solver, 2π)
    
    if !(solver.y ≈ solution(2π))
        println("failed, $(solver.y) != $(solution(2π)), $(report.idid)")
    end
    @assert solver.y ≈ solution(2π)
    # println(report.num_rejected_steps)
    # println(norm(solver.y))
    solver.vars.h
    solver.y

end

for i in 1:1000
    run()
    println("-------------------------------------------------------------------------------")
end