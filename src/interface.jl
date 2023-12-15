using DormandPrince.DP5Core: DP5Solver, dopri5

struct DP5Iterator{T <: Real}
    solver::DP5Solver
    times::AbstractVector{T}
end

# gets the first (t,y), return index which is the state
# here we choose 2 because 1 is the initial state which 
# is has been returned by the iterator
function Base.iterate(dp5_iterator::DP5Iterator)
    length(dp5_iterator.times) == 0 && return nothing # empty iterator
    # integrate to first time
    integrate(dp5_iterator.solver, first(dp5_iterator.times))
    # return value and index which is the state
    return (dp5_iterator.times[dp5_iterator.index], dp5_iterator.solver.y), 2    
end

# gets the next (t,y), return index+! which is the updated state
function Base.iterate(dp5_iterator::DP5Iterator, index::Int) 
    index > length(dp5_iterator.times) && return nothing # end of iterator
    # integrate to next time
    integrate(dp5_iterator.solver, dp5_iterator.times[index])
    # return time and state
    return (dp5_iterator.times[index], dp5_iterator.solver.y), index+1
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


