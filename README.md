# DormandPrince.jl

[![codecov](https://codecov.io/github/QuEraComputing/DormandPrince.jl/graph/badge.svg?token=qYZ4V7m0JY)](https://codecov.io/github/QuEraComputing/DormandPrince.jl)

Julia-native implementation of the Dormand-Prince 5th and 8th order solvers

## Installation

<p>
DormandPrince is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install DormandPrince,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd>
    key in the REPL to use the package mode, then type the following command
</p>

```julia
pkg> add DormandPrince
```

## Usage 

### Single End-Time

```julia
julia> using DormandPrince

# define your ODE, in this case, dy/dx = 0.85y
julia> function fcn(x, y, f)
            f[1] = 0.85y[1]
        end

# Create a solver object. We will use the 5th order solver and
# start integrating at t = 0.0 with initial value of 19.0
julia> solver = DP5Solver(fcn, 0.0, [19.0])

# Begin integration up to t = 2.1
julia> integrate!(solver, 2.1)

# get_current_state will return the answer
julia> get_current_state(solver)
```

Both the `DP5Solver` and `DP8Solver`'s are stateful allowing for memory-efficient integration to future end times from the last integrated end point (e.g. if you chose t = 1.0 as your endpoint you can call `integrate!` again with t=2.0 and it will "carry forward" the work starting from t = 1.0 instead of requiring you to set things up all over again).

### Multiple End-Times

```julia
julia> using DormandPrince

# Define ODE
julia> function fcn(x, y, f)
            f[1] = 0.85y[1]
        end

# Define times of interest to analyze/perform actions on the solution
julia> times = [1.0, 1.1, 1.9, 2.4]

# Create the solver object, with integration starting from t = 0.0 and initial value of 19.0
julia> solver = DP5Solver(fcn, 0.0, [19.0])

# Empty array to store intermediate values
julia> intermediate_values = []

# Use callback to store intermediate values. The solver object will also be mutated to store the solution 
# found at the last time point.
julia> integrate!(solver, times) do time, val
            push!(intermediate_values, val[])
        end

```

You may also create a `SolverIterator` that can then be iterated over. Note that this will require a fresh solver

```julia
for (time, value) in SolverIterator(solver, times)
    println("Time: ", time, " ", "Value: ", value[])
end
```

## License

Apache License 2.0
