using DifferentialEquations
call_count = 0

function f(u, p, t)
    global call_count += 1
    return -15.0 * u
end

#=-
@kwdef mutable struct CounterFunc{F}
    fun::F
    counter::Int64 = 0
end

function (f::CounterFunc)(x...)
    f.counter += 1
    return f.fun(x...)
end

counted_f = CounterFunc(fun = f)
=#

u0 = 1
tspan = (0.0, 0.5)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, DP5(), reltol = 1e-10, abstol = 1e-10)