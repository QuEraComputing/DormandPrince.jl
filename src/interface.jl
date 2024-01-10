"""
    struct SolverIterator{T <: Real}
    SolverIterator(solver, times)
    
Given a solver and a vector of times, 
this iterator will return the state of the solver at each time.

The solver will be mutated to hold the solution from the last time given.

# Examples

```julia
julia> solver = DP5Solver(fcn, 0.0, [0.0])

julia> times = [1.0, 2.0, 3.0]

julia> intermediate_vals = []

julia> for (time, value) in SolverIterator(solver, times)
            push!(intermediate_values, value[])
        end

```

"""
struct SolverIterator{T <: Real}
    solver::AbstractDPSolver{T}
    times::AbstractVector{T}
end

# gets the first (t,y), return index which is the state
# here we choose 2 because 1 is the initial state which 
# is has been returned by the iterator
function Base.iterate(solver_iter::SolverIterator)
    length(solver_iter.times) == 0 && return nothing # empty iterator
    # integrate to first time
    integrate!(solver_iter.solver, first(solver_iter.times))
    # return value and index which is the state
    return (solver_iter.times[1], get_current_state(solver_iter.solver)), 2    
end

# gets the next (t,y), return index+! which is the updated state
function Base.iterate(solver_iter::SolverIterator, index::Int) 
    index > length(solver_iter.times) && return nothing # end of iterator
    # integrate to next time
    integrate!(solver_iter.solver, solver_iter.times[index])
    # return time and state
    return (solver_iter.times[index], get_current_state(solver_iter.solver)), index+1
end

# 3 modes of operation for integrate
# 1. integrate(solver, time) -> state (modify solver object in place)
# 2. integrate(solver, times) -> iterator
# 3. integrate(callback, solver, times) -> vector of states with callback applied

get_current_state(::AbstractDPSolver) = error("not implemented")
integrate_core!(::AbstractDPSolver{T}, ::T) where T = error("not implemented")

"""
    integrate!(solver, time)

Integrate the ODE problem that is part of solver to the end time `time`.
`solver` can be a `DP5Solver` or a `DP8Solver` type.

The `solver` is mutated to hold the solution and `solver` state at the time `time`, allowing for 
efficient integration to future end times of interest through subsequent calls to `integrate!` with 
later times.

# Examples

```julia
julia> function fcn(x, y, f)
            f[1] = 0.85y[1]
        end

julia> solver = DP5Solver(fcn, 0.0, [19.0])

# Integrate from initial time of 0.0 to 1.0
julia> integrate!(solver, 1.0)

# integrate from the last time of 1.0 to 2.0
julia> integrate!(solver, 2.0)

# Solver object now holds solution and state at 2.0
julia> get_current_state(solver)
```
"""
function integrate!(solver::AbstractDPSolver{T}, time::T) where T <: Real
    report = integrate_core!(solver, time)
    if report.idid != COMPUTATION_SUCCESSFUL
        error("integration failed at time $time with report $report")
    end
end

"""
    integrate!(callback, solver, times)

Integrate the ODE problem that is part of `solver` to the end times `times` and apply the callback on 
the time and solution at each time. 
`solver` can be a `DP5Solver` or a `DP8Solver` type.

At the end of all the times the solver holds the solution and solver state at the last time in `times`.

# Examples

```julia
julia> function fcn(x, y, f)
            f[1] = 0.85y[1]
        end

julia> times = [1.0, 1.1, 1.9, 2.4]

julia> solver = DP5Solver(fcn, 0.0, [19.0])

julia> intermediate_values = []

julia> integrate!(solver, times) do time, val
            push!(intermediate_values, val[])
        end
```
"""
function integrate!(callback, solver::AbstractDPSolver{T}, times::AbstractVector{T}; sort_times::Bool = true) where {T <: Real}
    times = sort_times ? sort(collect(times)) : times

    result = []
    for time in times
        integrate!(solver, time)

        push!(result, callback(time, get_current_state(solver)))
    end

    return result
end


