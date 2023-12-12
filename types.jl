using Base.Iterators:repeated, Repeated

@enum Idid begin
    COMPUTATION_SUCCESSFUL = 1
    INPUT_NOT_CONSISTENT = -1 # use for check failures in the beginning of dopri5 call
    LARGER_NMAX_NEEDED = -2
    STEP_SIZE_BECOMES_TOO_SMALL = -3
end

@enum Checks begin
    CHECKS_SUCCESSFUL
    MAX_ALLOWED_STEPS_NEGATIVE
    UNSUPPORTED_UROUND
    CURIOUS_BETA 
    CURIOUS_SAFETY_FACTOR
end

struct DP5Report{T <: Real}
    x_final::T
    checks::Checks
    idid::Idid

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

    # should be either vector or repeated for type
    atol::Union{T, Vector{T}} = 1e-10
    rtol::Union{T, Vector{T}} = 1e-10

end

struct DP5Consts{T <: Real}
    expo1::T
    facc1::T
    facc2::T
    atol_iter::Union{Repeated{T}, Vector{T}}
    rtol_iter::Union{Repeated{T}, Vector{T}}

    function DP5Consts(options::DP5Options{T}) where {T <: Real}
        expo1 = 0.20-options.beta*0.75
        facc1 = 1.0/options.step_size_selection_one
        facc2 = 1.0/options.step_size_selection_two
        atol_iter = options.atol isa Number ? repeated(options.atol) : options.rtol
        rtol_iter = options.rtol isa Number ? repeated(options.rtol) : options.rtol
        new{T}(expo1, facc1, facc2,atol_iter,rtol_iter)
    end
end


@kwdef mutable struct DP5Vars{T<:Real}
    facold::T = 1e-4
    iasti::Int64 = 0
    nonsti::Int64 = 0
    hlamb::T = 0.0
    last::Bool = false
end

# should "dopri5" take in DP5Solver or should DP5Solver have some associated method
# attached to it? 
struct DP5Solver{StateType <: AbstractVector, T <: Real, F}
    f::F
    x::T
    current_h::T
    y::StateType
    k1::StateType
    k2::StateType
    k3::StateType
    k4::StateType
    k5::StateType
    k6::StateType
    y1::StateType
    ysti::StateType
    options::DP5Options{T}
    consts::DP5Consts{T}
    vars::DP5Vars{T}

    function DP5Solver(f::F, x::T, y::StateType; kw...) where {StateType <: AbstractVector, T<:Real, F}

        k1 = empty(y)
        k2 = empty(y)
        k3 = empty(y)
        k4 = empty(y)
        k5 = empty(y)
        k6 = empty(y)
        y1 = empty(y)
        ysti = empty(y)
        options = DP5Options(;kw...)
        consts = DP5Consts(options)
        vars = DP5Vars{T}()

        new{StateType, T, F}(f, x, options.initial_step_size, y, k1, k2, k3, k4, k5, k6, y1, ysti, options, consts, vars)
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