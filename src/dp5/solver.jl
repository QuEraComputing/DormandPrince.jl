# using Base.Iterators:repeated
#include("types.jl")
#include("checks.jl")
#include("helpers.jl")

function  DormandPrince.integrate(
   solver::DP5Solver{T},
   xend::T
) where T

    # check nmax, uround, safety factor, beta, safety_factor
    # just accept solver.options and handle accessing attributes internally
    #=
    arret = check_max_allowed_steps(solver.options) ||
            check_uround(solver.options) ||
            check_beta(solver.options) ||
            check_safety_factor(solver.options)
    =#
    check_max_allowed_steps(solver.options) || return Report(solver.vars.x, DormandPrince.MAX_ALLOWED_STEPS_NEGATIVE, DormandPrince.INPUT_NOT_CONSISTENT, 0, 0, 0, 0)
    check_uround(solver.options) || return Report(solver.vars.x, DormandPrince.UNSUPPORTED_UROUND, DormandPrince.INPUT_NOT_CONSISTENT, 0, 0, 0, 0)
    check_beta(solver.options) || return Report(solver.vars.x, DormandPrince.CURIOUS_BETA, DormandPrince.INPUT_NOT_CONSISTENT, 0, 0, 0, 0)
    check_safety_factor(solver.options) || return Report(solver.vars.x, DormandPrince.CURIOUS_SAFETY_FACTOR, DormandPrince.INPUT_NOT_CONSISTENT, 0, 0, 0, 0)

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


    h, report = dopri5(
        solver, # contains x, y, k1, k2, k3, k4, k5, k6, y1, ysti, options
        xend, 
        hmax, 
        h,
    )

    # update with final h
    solver.vars.h = h

    # reset the necessary vars 
    solver.vars.last = false

    return report

end

function dopri5(
    solver::DP5Solver{T}, # contains f, x, y, k1, k2, k3, k4, k5, k6, y1, ysti, options
    xend::T, 
    hmax::T,
    h::T, 
) where T
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
    
    idid = DormandPrince.LARGER_NMAX_NEEDED

    ###### Basic Integration Step
    for _ in 1:solver.options.maximum_allowed_steps
        # if nstep > solver.options.maximum_allowed_steps
        #     # GOTO 78
        #     # println(" MORE THAN NMAX = ", solver.options.maximum_allowed_steps, " STEPS ARE NEEDED")
        #     return h, Report(solver.vars.x, DormandPrince.INPUT_CHECKS_SUCCESSFUL, DormandPrince.LARGER_NMAX_NEEDED , 0, 0, 0, 0)
        # end
        
        if (0.10 * abs(h)) <= abs(solver.vars.x)*solver.options.uround 
            # GOTO 77
            # println("STEP SIZE TOO SMALL, H = ", h)
            # return h, Report(solver.vars.x, DormandPrince.INPUT_CHECKS_SUCCESSFUL, DormandPrince.STEP_SIZE_BECOMES_TOO_SMALL, 0, 0, 0, 0)

            idid = DormandPrince.STEP_SIZE_BECOMES_TOO_SMALL
            break
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
                idid = DormandPrince.COMPUTATION_SUCCESSFUL
                break
                # return h, Report(
                #     solver.vars.x, 
                #     DormandPrince.INPUT_CHECKS_SUCCESSFUL,
                #     DormandPrince.COMPUTATION_SUCCESSFUL, 
                #     nfcn, 
                #     nstep, 
                #     naccpt, 
                #     nrejct
                # )
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

    return h, Report(
        solver.vars.x, 
        DormandPrince.INPUT_CHECKS_SUCCESSFUL,
        idid, 
        nfcn, 
        nstep, 
        naccpt, 
        nrejct
    )
    
end
