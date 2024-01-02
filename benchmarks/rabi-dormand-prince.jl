using BenchmarkTools
using DormandPrince:DP5Solver, DP8Solver,  integrate

function fcn(x, y, f)
    g(x) = 2.2*2π*sin(2π*x)

    f[1] = -1im * g(x)/2 * y[2]
    f[2] = -1im * g(x)/2 * y[1]
end

solver = DP8Solver(
    fcn,
    0.0,
    ComplexF64[1.0, 0.0]
)

integrate(solver, 2π)

@benchmark  integrate(clean_solver, 2π) setup=(clean_solver = DP8Solver(fcn, 0.0, ComplexF64[1.0, 0.0])) samples=10000 evals=5 seconds=500