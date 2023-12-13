# recommended by Phillip Weinberg
# for profiling the number of calls to the ODE function

@kwdef mutable struct CounterFunc{F}
    fun::F
    counter::Int64 = 0
end

function (f::CounterFunc)(x...)
    f.counter += 1
    return f.fun(x...)
end