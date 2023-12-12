include("../solver.jl")
include("../types.jl")

# fcn(x,y,f)
## f is the end result of the derivative itself
## x is the dependent variable
## y is y(x)
function fcn(x,y,f)
    f[1] = -15.0*y[1]
end

@kwdef mutable struct CounterFunc{F}
    fun::F
    counter::Int64 = 0
end

function (f::CounterFunc)(x...)
    f.counter += 1
    return f.fun(x...)
end

f = CounterFunc(fun = fcn)

solver = DP5Solver(
    f,
    0.0, # initial x 
    [1.0] # initial y at x
)

dopri5(solver, 0.5)

