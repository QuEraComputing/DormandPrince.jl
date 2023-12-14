using DifferentialEquations
using BenchmarkTools


function f(du, u, p, t) 
    g(t) = 2.2*2π*sin(2π*t)

    du[1] = -1im * g(t)/2 * u[2]
    du[2] = -1im * g(t)/2 * u[1]
end

u0 = ComplexF64[1.0, 0.0]
tspan = (0.0, 2π)
prob = ODEProblem(f, u0, tspan)
# get precompilation out of the way
sol = solve(prob, DP5(), reltol=1e-10, abstol=1e-10)

# terminate benchmark after maximum of 5 minutes
@benchmark solve(prob, DP5(), reltol=1e-10, abstol=1e-10) samples=10000 evals=5 seconds=300