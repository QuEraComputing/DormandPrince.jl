
abstract type AbstractDPSolver{T <: Real, StateType <: AbstractVector, F} end

"""
    @enum Idid

Enum used to represent status of the integration process in the `Report` type.

# Values
- `COMPUTATION_SUCCESSFUL`: Integration completed successfully.
- `INPUT_NOT_CONSISTENT`: The options given to the solver violate acceptable ranges (see the 
`checks` field of the `Report` type for more information).
- `LARGER_NMAX_NEEDED`: The maximum number of allowed steps is too small.
- `STEP_SIZE_BECOMES_TOO_SMALL`: The step size becomes too small.
- `STEP_SIZE_BECOMES_NAN`: The step size becomes NaN.
"""
@enum Idid begin
    COMPUTATION_SUCCESSFUL = 1
    INPUT_NOT_CONSISTENT = -1 # use for check failures in the beginning of  core_integrator call
    LARGER_NMAX_NEEDED = -2
    STEP_SIZE_BECOMES_TOO_SMALL = -3
    STEP_SIZE_BECOMES_NAN = -4
end

""" 
    @enum Checks

Enum used to represent status of checks on user options to the solver in the `Report` type.

# Values
- `INPUT_CHECKS_SUCCESSFUL`: All checks on user options passed.
- `MAX_ALLOWED_STEPS_NEGATIVE`: The maximum allowed steps value is negative.
- `UNSUPPORTED_UROUND`: The uround value is either too small (<= 1e-35) or too large (>= 1.0).
- `CURIOUS_BETA`: The Beta value for stabilized step control is greater than 0.2.
- `CURIOUS_SAFETY_FACTOR`: The safety factor is either too large (>= 1.0) or too small (<= 1e-4).
"""
@enum Checks begin
    INPUT_CHECKS_SUCCESSFUL
    MAX_ALLOWED_STEPS_NEGATIVE
    UNSUPPORTED_UROUND
    CURIOUS_BETA 
    CURIOUS_SAFETY_FACTOR
end

"""
    struct Report{T <: Real}

Contains data on the integration process including after integration has completed and if an error was detected within 
options provided for integration. 

# Fields
- `x_final`: The final time point integration reached. Should be equal to the time point provided by the user with successful usage.
- `checks`: Status of checks performed on the options provided by the user, represented by a `Checks` element. Should be `INPUT_CHECKS_SUCCESSFUL` with successful usage.
- `idid`: Status of the integration process, represented by an `Idid` element. Should be `COMPUTATION_SUCCESSFUL` with successful usage.
- `num_func_evals`: Number of function evaluations performed during integration.
- `num_computed_steps`: Number of computed steps during integration. 
- `num_accpeted_steps`: Number of accepted steps during integration.
- `num_rejected_steps`: Number of rejected steps during integration. 
- 
"""
struct Report{T <: Real}
    x_final::T
    checks::Checks
    idid::Idid

    num_func_evals::Int
    num_computed_steps::Int
    num_accepted_steps::Int
    num_rejected_steps::Int
end

"""
    struct Options{T <: Real}

Holds the options used in the integration process by a solver object.

# Fields
- `atol`: Absolute Tolerance. Default is 1e-10
- `rtol`: Relative Tolerance. Default is 1e-10
- `uround`: Rounding unit, default is eps(T) where T <: Real.
- `safety_factor`: Safety factor in step size prediction, default is 0.9.
- Step Size Selection parameters, with the new step size being subject to the constraint: `step_size_selection_one <= new_step/old_step <= step_size_selection_two`
  - `step_size_selection_one`: Defaut is 0.2
  - `step_size_selection_two`: Default is 10.0
- `beta`: Beta for stabilized step size control. Large values of beta (<= 0.1) make step size control more stable.
- `maximal_step_size`: The largest the step size can be. Default is 0.0, which internally translates to `xend - x`.
- `initial_step_size`: Initial step size, default is 0.0. An initial guess is computed internally. 
- `maximum_allowed_steps`: Maximum number of allowed steps, default is 100000.
- `print_error_messages`: Whether to print error messages, default is true.
- `stiffness_test_activation_step`: After integer multiples of this number of steps, perform stiffness detection. Default is 1000.

"""
@option struct Options{T <: Real}
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

struct Consts{T <: Real, Tol <: Union{Repeated{T}, Vector{T}}}
    expo1::T
    facc1::T
    facc2::T
    atol_iter::Tol
    rtol_iter::Tol
end

function Consts(expo1::T, options::Options{T}) where T
    facc1 = 1.0/options.step_size_selection_one
    facc2 = 1.0/options.step_size_selection_two
    atol_iter = options.atol isa T ? repeated(options.atol) : collect(T, options.atol)
    rtol_iter = options.rtol isa T ? repeated(options.rtol) : collect(T, options.rtol)
    Consts{T, typeof(atol_iter)}(expo1, facc1, facc2, atol_iter, rtol_iter)
end


@option mutable struct Vars{T <: Real}
    x::T = zero(T)
    xph::T = zero(T)
    h::T = zero(T)
    facold::T = 1e-4
    iasti::Int = 0
    nonsti::Int = 0
    hlamb::T = zero(T)
    last::Bool = false
end


"""
    struct DP5Solver

A 5th Order Dormand-Prince solver object that contains:
- the ODE problem of interest
- The solution to the ODE at the last integrated-to time point
- intermediate arrays used in the Runge-Kutta method
- constants and variables used by the solver
- user-provided options for the solver

The storage of such intermediate information allows for efficient integration from a previously integrated-to 
time point and a future time point.
"""
struct DP5Solver{T, StateType , F, OptionsType, ConstsType, VarsType} <: AbstractDPSolver{T, StateType, F}
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
    options::OptionsType
    consts::ConstsType
    vars::VarsType

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
        expo1 = 0.20 - options.beta * 0.75
        consts = Consts(expo1, options)
        vars = Vars{T}(;x=x, h=options.initial_step_size)

        new{T, StateType, F, typeof(options), typeof(consts), typeof(vars)}(
            f, y, k1, k2, k3, k4, k5, k6, y1, ysti, options, consts, vars
        )
    end
end

""" 
    DP5Solver(f, x, y; kw...)

Create a 5th Order Dormand-Prince solver object to solve the ODE problem `y' = f(x, y)`.

# Examples

```julia
julia> fcn(x, y, f)
            f[1] = 0.85y[1]
        end

julia> solver = DP5Solver(fcn, 0.0, [19.0]; atol = 1e-10, rtol = 1e-10)
```

# Arguments
- `f`: The function representing the ODE, should be in the form `f(x, y, f)`.
- `x`: The starting time point of the ODE problem.
- `y`: The initial value of the ODE in vector form

# Keyword Arguments
- `atol`: Absolute Tolerance. Default is 1e-10
- `rtol`: Relative Tolerance. Default is 1e-10
- `uround`: Rounding unit, default is eps(T) where T <: Real.
- `safety_factor`: Safety factor in step size prediction, default is 0.9.
- Step Size Selection parameters, with the new step size being subject to the constraint: `step_size_selection_one <= new_step/old_step <= step_size_selection_two`
  - `step_size_selection_one`: Defaut is 0.2
  - `step_size_selection_two`: Default is 10.0
- `beta`: Beta for stabilized step size control. Large values of beta (<= 0.1) make step size control more stable.
- `maximal_step_size`: The largest the step size can be. Default is 0.0, which internally translates to `xend - x`.
- `initial_step_size`: Initial step size, default is 0.0. An initial guess is computed internally. 
- `maximum_allowed_steps`: Maximum number of allowed steps, default is 100000.
- `print_error_messages`: Whether to print error messages, default is true.
- `stiffness_test_activation_step`: After integer multiples of this number of steps, perform stiffness detection. Default is 1000.
"""
function DP5Solver(
    f, 
    x::Real, 
    y::AbstractVector; kw...)
    k1 = similar(y)
    k2 = similar(y)
    k3 = similar(y)
    k4 = similar(y)
    k5 = similar(y)
    k6 = similar(y)
    y1 = similar(y)
    ysti = copy(y)
    DP5Solver(f, x, y, k1, k2, k3, k4, k5, k6, y1, ysti;kw...)
end


# should " core_integrator" take in DP5Solver or should DP5Solver have some associated method
# attached to it? 
"""
    struct DP8Solver

An 8th Order Dormand-Prince solver object that contains:
- the ODE problem of interest
- The solution to the ODE at the last integrated-to time point
- intermediate arrays used in the Runge-Kutta method
- constants and variables used by the solver
- user-provided options for the solver

The storage of such intermediate information allows for efficient integration from a previously integrated-to 
time point and a future time point.
"""
struct DP8Solver{T, StateType ,F, OptionsType, ConstsType, VarsType} <: AbstractDPSolver{T, StateType, F}
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
    options::OptionsType
    consts::ConstsType
    vars::VarsType

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
        expo1 = 0.125 - options.beta * 0.2
        consts = Consts(expo1, options)
        vars = Vars{T}(;x=x, h=options.initial_step_size)

        new{T, StateType, F, typeof(options), typeof(consts), typeof(vars)}(
            f, y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, y1, options, consts, vars
        )
    end
end

""" 
    DP8Solver(f, x, y; kw...)

Create an 8th Order Dormand-Prince solver object to solve the ODE problem `y' = f(x, y)`.

# Examples

```julia
julia> fcn(x, y, f)
            f[1] = 0.85y[1]
        end

julia> solver = DP8Solver(fcn, 0.0, [19.0]; atol = 1e-10, rtol = 1e-10)
```

# Arguments
- `f`: The function representing the ODE, should be in the form `f(x, y, f)`.
- `x`: The starting time point of the ODE problem.
- `y`: The initial value of the ODE in vector form

# Keyword Arguments
- `atol`: Absolute Tolerance. Default is 1e-10
- `rtol`: Relative Tolerance. Default is 1e-10
- `uround`: Rounding unit, default is eps(T) where T <: Real.
- `safety_factor`: Safety factor in step size prediction, default is 0.9.
- Step Size Selection parameters, with the new step size being subject to the constraint: `step_size_selection_one <= new_step/old_step <= step_size_selection_two`
  - `step_size_selection_one`: Defaut is 0.2
  - `step_size_selection_two`: Default is 10.0
- `beta`: Beta for stabilized step size control. Large values of beta (<= 0.1) make step size control more stable.
- `maximal_step_size`: The largest the step size can be. Default is 0.0, which internally translates to `xend - x`.
- `initial_step_size`: Initial step size, default is 0.0. An initial guess is computed internally. 
- `maximum_allowed_steps`: Maximum number of allowed steps, default is 100000.
- `print_error_messages`: Whether to print error messages, default is true.
- `stiffness_test_activation_step`: After integer multiples of this number of steps, perform stiffness detection. Default is 1000.
"""
function DP8Solver(
    f, 
    x::Real, 
    y::AbstractVector; kw...)
    k1 = similar(y)
    k2 = similar(y)
    k3 = similar(y)
    k4 = similar(y)
    k5 = similar(y)
    k6 = similar(y)
    k7 = similar(y)
    k8 = similar(y)
    k9 = similar(y)
    k10 = similar(y)
    y1 = similar(y)
    DP8Solver(f, x, y, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, y1;kw...)
end

"""
    get_current_state(solver::DP5Solver)

Gets the current state/solution held by the `DP5Solver` object. 
"""
get_current_state(solver::DP5Solver) = solver.y
"""
    get_current_state(solver::DP8Solver)

Gets the current state/solution held by the `DP8Solver` object. 
"""
get_current_state(solver::DP8Solver) = solver.y