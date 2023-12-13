using Base.Iterators:repeated
include("types.jl")
include("checks.jl")
include("helpers.jl")

function dopri5(
   solver,
   xend
)

    # check nmax, uround, safety factor, beta, safety_factor
    # just accept solver.options and handle accessing attributes internally
    #=
    arret = check_max_allowed_steps(solver.options) ||
            check_uround(solver.options) ||
            check_beta(solver.options) ||
            check_safety_factor(solver.options)
    =#
    check_max_allowed_steps(solver.options) || return DP5Report(solver.vars.x, MAX_ALLOWED_STEPS_NEGATIVE, INPUT_NOT_CONSISTENT, 0, 0, 0, 0)
    check_uround(solver.options) || return DP5Report(solver.vars.x, UNSUPPORTED_UROUND, INPUT_NOT_CONSISTENT, 0, 0, 0, 0)
    check_beta(solver.options) || return DP5Report(solver.vars.x, CURIOUS_BETA, INPUT_NOT_CONSISTENT, 0, 0, 0, 0)
    check_safety_factor(solver.options) || return DP5Report(solver.vars.x, CURIOUS_SAFETY_FACTOR, INPUT_NOT_CONSISTENT, 0, 0, 0, 0)

    ###### nstiff -  parameters for stiffness detection
    # nstiff = solver_options.stiffness_test_activation_step

    if solver.options.stiffness_test_activation_step < 0
        solver.options.stiffness_test_activation_step = solver.options.maximum_allowed_steps + 10
    end 

    ####### maximal step size
    hmax = if iszero(solver.options.maximal_step_size)
        xend-solver.vars.x
    else
        solver.options.maximal_step_size
    end

    ####### initial step size
    #h = solver_options.initial_step_size
    h = solver.vars.h

    #=
    Total Storage Requirement check used to live here but
    is no longer needed because the DP5Solver constructor
    ensures enough space is allocated based on the size of
    y :D 
    =#


    h, dp5_report = dopcor(
        solver, # contains x, y, k1, k2, k3, k4, k5, k6, y1, ysti, options
        xend, 
        hmax, 
        h,
    )

    # update with final h
    solver.vars.h = h

    return dp5_report

end

function dopcor(
    solver, # contains f, x, y, k1, k2, k3, k4, k5, k6, y1, ysti, options
    xend, 
    hmax,
    h, 
)
    ##### Initializations
    # replace sign with Julia-native Base.sign
    # posneg = sign(1.0, xend-solver.vars.x)
    posneg = 1.0 * Base.sign(xend-solver.vars.x)

    ###### Initial Preparations
    nfcn   = 0
    nstep  = 0
    naccpt = 0
    nrejct = 0

    solver.f(solver.vars.x, solver.y, solver.k1)
    hmax = abs(hmax)
    iord = 5

    # may be considered type unstable
    if iszero(h)
        # solver contains, fcn, x, y, k1, k2, k3, atol, rtol
        h = hinit(solver, posneg, iord, hmax)
    end

    nfcn += 2
    reject = false
    
    ###### Basic Integration Step
    while true
        if nstep > solver.options.maximum_allowed_steps
            # GOTO 78
            # println(" MORE THAN NMAX = ", solver.options.maximum_allowed_steps, " STEPS ARE NEEDED")
            return h, DP5Report(solver.vars.x, CHECKS_SUCCESSFUL, LARGER_NMAX_NEEDED , 0, 0, 0, 0)
        end
        
        if (0.10 * abs(h)) <= abs(solver.vars.x)*solver.options.uround 
            # GOTO 77
            # println("STEP SIZE TOO SMALL, H = ", h)
            return h, DP5Report(solver.vars.x, CHECKS_SUCCESSFUL, STEP_SIZE_BECOMES_TOO_SMALL, 0, 0, 0, 0)
        end

        if ((solver.vars.x+1.01*h-xend)*posneg) > 0.0
            h = xend-solver.vars.x
            solver.vars.last = true
        end
        
        nstep += 1

        do_step!(solver, h)

        nfcn += 6

        ###### Error Estimation
        err = error_estimation(solver)

        ###### Computation of hnew
        fac11 = err^solver.consts.expo1
        ###### Lund-Stabilization
        fac = fac11/(solver.vars.facold^solver.options.beta)
        ###### we require fac1 <= hnew/h <= fac2
        fac = max(solver.consts.facc2, min(solver.consts.facc1, fac/solver.options.safety_factor)) # facc1, facc2, fac must be Float64 
        hnew = h/fac
        if err <= 1.0 
            ###### Step is accepted
            solver.vars.facold = max(err, 1e-4)
            naccpt += 1
            ###### Stiffness Detection
            stiffness_detection!(solver, naccpt, h)

            solver.k1 .= solver.k2
            solver.y  .= solver.y1

            solver.vars.x += h

            ###### Normal Exit
            if solver.vars.last 
                h = hnew
                return h, DP5Report(
                    solver.vars.x, 
                    CHECKS_SUCCESSFUL,
                    COMPUTATION_SUCCESSFUL, 
                    nfcn, 
                    nstep, 
                    naccpt, 
                    nrejct
                )
            end

            if(abs(hnew) > hmax)
                hnew = posneg*hmax
            end

            if reject
                hnew = posneg*min(abs(hnew), abs(h))
            end
            reject = false
        else
            ###### Step is rejected
            hnew = h/min(solver.consts.facc1, fac11/solver.options.safety_factor)
            reject = true
            if naccpt > 1
                nrejct += 1
            end
            solver.vars.last = false

        end
        h = hnew
    end

end

function hinit(
    solver, 
    posneg, 
    iord,
    hmax
    # f0 arg is k1 from dopcor
    # f1 arg is k2 from dopcor
    # y1 arg is k3 from dopcor
)
    #=
    Compute a first guess for explicit euler as
        h = 0.01 * norm (y0) / norm (f0)
    the increment for explicit euler is small
    compared to the solution
    =#
    h, dnf = euler_first_guess(solver, hmax, posneg)

    ###### Perform an explicit step
    #y1 = y + h*f0
    #fcn(n, x+h, y1, f1)
    copyto!(solver.y1, solver.y + h*solver.k1)
    solver.f(solver.vars.x + h, solver.k3, solver.k2)

    ###### Estimate the second derivative of the solution
    der2 = estimate_second_derivative(solver, h)

    ##### Step size is computed such that
    ##### H**IORD * MAX ( NORM(F0), NORM(F1), DER2 ) = 0.01
    der12 = max(abs(der2), sqrt(dnf))
    if der12 <= 1e-15
        h1 = max(1.0e-6, abs(h)*1.0e-3)
    else
        h1 = (0.01/der12)^(1.0/iord)
    end
    h = min(100*abs(h), h1, hmax)
    return h * Base.sign(posneg)

end