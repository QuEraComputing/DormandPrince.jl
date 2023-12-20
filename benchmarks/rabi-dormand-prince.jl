using BenchmarkTools
using DormandPrince:DP5Solver,  core_integrator

function fcn(x, y, f)
    g(x) = 2.2*2π*sin(2π*x)

    f[1] = -1im * g(x)/2 * y[2]
    f[2] = -1im * g(x)/2 * y[1]
end

solver = DP5Solver(
    fcn,
    0.0,
    ComplexF64[1.0, 0.0]
)

 core_integrator(solver, 2π)

@benchmark  core_integrator(clean_solver, 2π) setup=(clean_solver = DP5Solver(fcn, 0.0, ComplexF64[1.0, 0.0])) samples=10000 evals=5 seconds=500