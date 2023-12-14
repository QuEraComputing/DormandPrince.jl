mutable struct DP5Iterator{T <: Real}
    solver::DP5Solver
    times::AbstractVector{T}
    index::Int

    function DP5Iterator(solver::DP5Solver, times::AbstractVector{T}) where {T <: Real}
        new{T}(solver, times, 1)
    end
end

# gets the first (t,y) so we DO NOT increment the index
function Base.iterate(dp5_iterator::DP5Iterator)
    # integrate to next time
    if dp5_iterator.index <= length(dp5_iterator.times)
        # integrate to next time
        dopri5(dp5_iterator.solver, dp5_iterator.times[dp5_iterator.index])
        # return time and state
        return (dp5_iterator.times[dp5_iterator.index], dp5_iterator.solver.y), nothing
    else
        return nothing
    end
    
end

# don't really need the state here because we can just acccess it from the solver
# gets subsequent (t, y) so we SHOULD increment the index
function Base.iterate(dp5_iterator::DP5Iterator, state) 
    if dp5_iterator.index < length(dp5_iterator.times)
        dp5_iterator.index += 1
        # integrate to next time
        dopri5(dp5_iterator.solver, dp5_iterator.times[dp5_iterator.index])
        # return time and state
        return (dp5_iterator.times[dp5_iterator.index], dp5_iterator.solver.y), nothing
    else
        return nothing
    end
end

# 3 modes of operation for integrate
# 1. integrate(solver, time) -> state (modify solver object in place)
# 2. integrate(solver, times) -> iterator
# 3. integrate(callback, solver, times) -> vector of states with callback applied

integrate(solver::DP5Solver, time::Real) = dopri5(solver, time)
integrate(solver::DP5Solver, times::AbstractVector{T}) where {T <: Real} = DP5Iterator(solver, times)


function integrate(callback, solver::DP5Solver{StateVec, T}, times::AbstractVector{T}; sort_times::Bool = true) where {StateVec <: AbstractVector, T <: Real}
    times = sort_times ? sort(collect(times)) : times

    result = []
    for time in times
        integrate(solver, time)
        push!(result, callback(time, solver.y))
    end

    return result
end


