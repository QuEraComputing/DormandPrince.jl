
struct SolverIterator{T <: Real}
    solver::AbstractDPSolver
    times::AbstractVector{T}
end

# gets the first (t,y), return index which is the state
# here we choose 2 because 1 is the initial state which 
# is has been returned by the iterator
function Base.iterate(solver_iter::SolverIterator)
    length(solver_iter.times) == 0 && return nothing # empty iterator
    # integrate to first time
    integrate(solver_iter.solver, first(solver_iter.times))
    # return value and index which is the state
    return (solver_iter.times[1], get_current_state(solver_iter.solver)), 2    
end

# gets the next (t,y), return index+! which is the updated state
function Base.iterate(solver_iter::SolverIterator, index::Int) 
    index > length(solver_iter.times) && return nothing # end of iterator
    # integrate to next time
    integrate(solver_iter.solver, solver_iter.times[index])
    # return time and state
    return (solver_iter.times[index], get_current_state(solver_iter.solver)), index+1
end

# 3 modes of operation for integrate
# 1. integrate(solver, time) -> state (modify solver object in place)
# 2. integrate(solver, times) -> iterator
# 3. integrate(callback, solver, times) -> vector of states with callback applied

get_current_state(::AbstractDPSolver) = error("not implemented")
integrate(solver::AbstractDPSolver{T}, times::AbstractVector{T}) where {T <: Real} = SolverIterator(solver, times)

function integrate(callback, solver::AbstractDPSolver{T}, times::AbstractVector{T}; sort_times::Bool = true) where {T <: Real}
    times = sort_times ? sort(collect(times)) : times

    result = []
    for time in times
        integrate(solver, time)
        push!(result, callback(time, get_current_state(solver)))
    end

    return result
end


