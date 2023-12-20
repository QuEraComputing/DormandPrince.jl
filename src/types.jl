
abstract type AbstractDPSolver{T <: Real, StateType <: AbstractVector, F} end

@enum Idid begin
    COMPUTATION_SUCCESSFUL = 1
    INPUT_NOT_CONSISTENT = -1 # use for check failures in the beginning of  core_integrator call
    LARGER_NMAX_NEEDED = -2
    STEP_SIZE_BECOMES_TOO_SMALL = -3
end

@enum Checks begin
    INPUT_CHECKS_SUCCESSFUL
    MAX_ALLOWED_STEPS_NEGATIVE
    UNSUPPORTED_UROUND
    CURIOUS_BETA 
    CURIOUS_SAFETY_FACTOR
end

struct Report{T <: Real}
    x_final::T
    checks::Checks
    idid::Idid

    num_func_evals::Int
    num_computed_steps::Int
    num_accepted_steps::Int
    num_rejected_steps::Int
end

@kwdef struct Options{T <: Real}
    # originally in work[1] - work[7]
    uround::T = eps(T)
    safety_factor:: T = 0.9
    step_size_selection_one::T = 0.2
    step_size_selection_two::T = 10.0
    beta::T = 0.04
    maximal_step_size::T = 0.0
    initial_step_size::T = 0.0

    # originally in iwork[1] - iwork[4]
    maximum_allowed_steps::Int = 100000
    print_error_messages::Bool = true
    stiffness_test_activation_step::Int = 1000

    # should be either vector or repeated for type
    atol::Union{T, Vector{T}} = 1e-10
    rtol::Union{T, Vector{T}} = 1e-10

end

struct Consts{T <: Real}
    expo1::T
    facc1::T
    facc2::T
    atol_iter::Union{Repeated{T}, Vector{T}}
    rtol_iter::Union{Repeated{T}, Vector{T}}

    function Consts(options::Options{T}) where {T <: Real}
        expo1 = 0.20-options.beta*0.75
        facc1 = 1.0/options.step_size_selection_one
        facc2 = 1.0/options.step_size_selection_two
        atol_iter = options.atol isa Number ? repeated(options.atol) : options.rtol
        rtol_iter = options.rtol isa Number ? repeated(options.rtol) : options.rtol
        new{T}(expo1, facc1, facc2,atol_iter,rtol_iter)
    end
end


@kwdef mutable struct Vars{T <: Real}
    x::T = zero(T)
    h::T = zero(T)
    facold::T = 1e-4
    iasti::Int = 0
    nonsti::Int = 0
    hlamb::T = zero(T)
    last::Bool = false
end


# should " core_integrator" take in DP5Solver or should DP5Solver have some associated method
# attached to it? 
struct DP5Solver{T, StateType ,F} <: AbstractDPSolver{T, StateType, F}
    f::F
    y::StateType
    k1::StateType
    k2::StateType
    k3::StateType
    k4::StateType
    k5::StateType
    k6::StateType
    y1::StateType
    ysti::StateType
    options::Options{T}
    consts::Consts{T}
    vars::Vars{T}

    function DP5Solver(
        f::F, 
        x::T, 
        y::StateType,
        k1::StateType,
        k2::StateType,
        k3::StateType,
        k4::StateType,
        k5::StateType,
        k6::StateType,
        y1::StateType,
        ysti::StateType; kw...) where {T <: Real, StateType <: AbstractVector, F}

        #TODO: check if y, k1, k2, k3, k4, k5, k6, y1, ysti have the same length

        options = Options{T}(;kw...)
        consts = Consts(options)
        vars = Vars{T}(;x=x, h=options.initial_step_size)

        new{T, StateType, F}(f, y, k1, k2, k3, k4, k5, k6, y1, ysti, options, consts, vars)
    end
end

function DP5Solver(
    f, 
    x::Real, 
    y::AbstractVector; kw...)
    k1 = copy(y)
    k2 = copy(y)
    k3 = copy(y)
    k4 = copy(y)
    k5 = copy(y)
    k6 = copy(y)
    y1 = copy(y)
    ysti = copy(y)
    DP5Solver(f, x, y, k1, k2, k3, k4, k5, k6, y1, ysti;kw...)
end


# should " core_integrator" take in DP5Solver or should DP5Solver have some associated method
# attached to it? 
struct DP8Solver{T, StateType ,F} <: AbstractDPSolver{T, StateType, F}
    f::F
    y::StateType
    k1::StateType
    k2::StateType
    k3::StateType
    k4::StateType
    k5::StateType
    k6::StateType
    k7::StateType
    k8::StateType
    k9::StateType
    k10::StateType
    y1::StateType
    options::Options{T}
    consts::Consts{T}
    vars::Vars{T}

    function DP8Solver(
        f::F, 
        x::T, 
        y::StateType,
        k1::StateType,
        k2::StateType,
        k3::StateType,
        k4::StateType,
        k5::StateType,
        k6::StateType,
        k7::StateType,
        k8::StateType,
        k9::StateType,
        k10::StateType,
        y1::StateType; 
        # overwrite default options with explicit kw
        beta::T = 0.0,
        step_size_selection_one::T = 0.333,
        step_size_selection_two::T = 6.0,
        kw...) where {T <: Real, StateType <: AbstractVector, F}

        #TODO: check if y, k1, k2, k3, k4, k5, k6, y1, ysti have the same length

        options = Options{T}(;
            beta=beta, 
            step_size_selection_one=step_size_selection_one, 
            step_size_selection_two=step_size_selection_two, 
            kw...
        )
        consts = Consts(options)
        vars = Vars{T}(;x=x, h=options.initial_step_size)

        new{T, StateType, F}(f, y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, y1, options, consts, vars)
    end
end

function DP8Solver(
    f, 
    x::Real, 
    y::AbstractVector; kw...)
    k1 = copy(y)
    k2 = copy(y)
    k3 = copy(y)
    k4 = copy(y)
    k5 = copy(y)
    k6 = copy(y)
    k7 = copy(y)
    k8 = copy(y)
    k9 = copy(y)
    k10 = copy(y)
    y1 = copy(y)
    DP8Solver(f, x, y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, y1;kw...)
end

