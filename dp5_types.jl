struct DP5Report 
    x_final::Float64
    idid::Int64

    num_func_evals::Int64
    num_computed_steps::Int64
    num_accepted_steps::Int64
    num_rejected_steps::Int64
end

@kwdef struct DP5Options{T <: Real}
    # originally in work[1] - work[7]
    uround::T = 2.3e-16
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
end


mutable struct DP5Solver{StateType <: AbstractVector, T <: Real}
    f::Function
    x::T
    x_end::T
    y::StateType
    k1::StateType
    k2::StateType
    k3::StateType
    k4::StateType
    k5::StateType
    k6::StateType
    options::DP5Options

    function DP5Solver(f::Function, x::T, x_end::T, y; kw...) where {StateType <: AbstractVector, T<:Real}

        k1 = copy(y)
        k2 = copy(y)
        k3 = copy(y)
        k4 = copy(y)
        k5 = copy(y)
        k6 = copy(y)
        new{StateType, T}(f, x, x_end, y, k1, k2, k3, k4, k5, k6, DP5Options(;kw...))
    end
end








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