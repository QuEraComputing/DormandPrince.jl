
struct DP5Report{T <: Real}
    x_final::T
    idid::Int64

    num_func_evals::Int64
    num_computed_steps::Int64
    num_accepted_steps::Int64
    num_rejected_steps::Int64
end

@kwdef struct DP5Options{T <: Real}
    # originally in work[1] - work[7]
    uround::T = eps(T)
    safety_factor:: T = 0.9
    step_size_selection_one::T = 0.2
    step_size_selection_two::T = 10.0
    beta::T = 0.04
    maximal_step_size::T = 0.0
    initial_step_size::T = 0.0

    # originally in iwork[1] - iwork[4]
    maximum_allowed_steps::Int64 = 100000
    print_error_messages::Bool = true
    stiffness_test_activation_step::Int64 = 1000

    atol::Union{T, Vector{T}} = 1e-10
    rtol::Union{T, Vector{T}} = 1e-10
end

# should "dopri5" take in DP5Solver or should DP5Solver have some associated method
# attached to it? 
mutable struct DP5Solver{StateType <: AbstractVector, T <: Real}
    f::Function
    x::T
    y::StateType
    k1::StateType
    k2::StateType
    k3::StateType
    k4::StateType
    k5::StateType
    k6::StateType
    y1::StateType
    ysti::StateType
    options::DP5Options

    function DP5Solver(f::Function, x::T, y::StateType; kw...) where {StateType <: AbstractVector, T<:Real}

        k1 = copy(y)
        k2 = copy(y)
        k3 = copy(y)
        k4 = copy(y)
        k5 = copy(y)
        k6 = copy(y)
        y1 = copy(y)
        ysti = copy(y)
        new{StateType, T}(f, x, y, k1, k2, k3, k4, k5, k6, y1, ysti, DP5Options(;kw...))
    end
end

#=
function integrate(solver::DP5Solver{StateVec, T}, x_end::T) where {StateVec, T}
    # integrate from solver.x to x_end
    ...
end

# focus on first two interfaces for integrator, 
# 

solver = DP5Solver(fcn,0,y0)
# low level interface
integrate(solver, 0.1)
# do stuff here
integrate(solver, 0.2)
# do stuff here
integrate(solver, 0.3)

# iterator interface
# TODO: Check how to dispatch to Iterator of a specific type objects
iterator = integrate(solver, times)

for t,y in iterator:
    #do stuff
end

# callback for each time in times
integrate(callback, solver, times) -> Vector

results = integrate(solver, times) do t, state
    # do stuff with t and state
end



# 99% of use cases:
function integrate(callback, solver::DP5Solver{StateVec, T}, times::AbstractVector{T}; sort_times::Bool = true) where {StateVec, T}
    times = sort_times ? sorted(collect(times)) : times

    result = []
    for time in times
        integrate(solver, time)
        push!(result, callback(time, solver.y))
    end
end
=#


# DP5Options() = DP5Options(
#     2.3e-16, #uround
#     0.9,     # safety_factor
#     0.2,     # step_size_selection_one
#     10.0,    # step_size_selection_two
#     0.04,    # beta
#     0.0,     # maximal step size - default to 0.0, later set to xend - x
#     0.0,     # initial step size - default to 0.0, trigger hinit later
#     100000,  # maximum number of allowed steps 
#     true,    # whether or not error messages should be printed
#     1000 # stiffness test activated after step J * this number
# )