using Base.Iterators:repeated
include("dp5_types.jl")
include("dp5_checks.jl")
#=
interface will be integrate(solver, x_end)
=#
function dopri5(
   solver::DP5Solver,
   xend
)
    arret = false


    ###### nmax - maximal number of steps
    if solver.options.maximum_allowed_steps < 0
        if solver.options.print_error_messages
            println("Maximum Allowed steps cannot be negative")
        end
        arret = true
    end

    ###### nstiff -  parameters for stiffness detection
    # nstiff = solver_options.stiffness_test_activation_step

    if solver.options.stiffness_test_activation_step < 0
        solver.options.stiffness_test_activation_step = solver.options.maximum_allowed_steps + 10
    end 

    ###### uround - smallest number satisfying 1.0 + uround > 1.0

    if (solver.options.uround <= 1e-35) || (solver.options.uround >= 1.0)
        if solver.options.print_error_messages
            println("WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:", solver.options.uround)
        end
        arret = true
    end 

    ####### safety factor
    #safe = solver_options.safety_factor 
    if (solver.options.safety_factor >= 1.0) || (solver.options.safety_factor <= 1e-4)
        if solver.options.print_error_messages
            println("CURIOUS INPUT FOR SAFETY FACTOR WORK[2]=", solver.options.safety_factor)
        end
        arret = true
    end
   
    ###### beta - for step control stabilization

    if solver.options.beta > 0.2
        if solver.options.print_error_messages
            println("CURIOUS INPUT FOR BETA: ", solver.options.beta)
        end
        arret = true
    end

    ####### maximal step size

    if solver.options.maximal_step_size == 0.0
        hmax = xend-solver.x
    else
        hmax = solver.options.maximal_step_size
    end

    ####### initial step size
    #h = solver_options.initial_step_size
    h = solver.current_h

    #=
    Total Storage Requirement check used to live here but
    is no longer needed because the DP5Solver constructor
    ensures enough space is allocated based on the size of
    y :D 
    =#


    ####### When a fail has occured, return with idid = -1
    if arret
        return DP5REport(x, -1, 0, 0, 0, 0)
    end

    h, dp5_report = dopcor(
        solver, # contains x, y, k1, k2, k3, k4, k5, k6, y1, ysti, options
        xend, 
        hmax, 
        h,
    )

    # update with final h
    solver.current_h = h

    return dp5_report

end

function dopcor(
    solver, # contains f, x, y, k1, k2, k3, k4, k5, k6, y1, ysti, options
    xend, 
    hmax,
    h, 
)
    ##### Initializations

    coeffs = cdopri()
    c2 = coeffs[1]
    c3 = coeffs[2]
    c4 = coeffs[3]
    c5 = coeffs[4]
    a21 = coeffs[5]
    a31 = coeffs[6]
    a32 = coeffs[7]
    a41 = coeffs[8]
    a42 = coeffs[9]
    a43 = coeffs[10]
    a51 = coeffs[11]
    a52 = coeffs[12]
    a53 = coeffs[13]
    a54 = coeffs[14]
    a61 = coeffs[15]
    a62 = coeffs[16]
    a63 = coeffs[17]
    a64 = coeffs[18]
    a65 = coeffs[19]
    a71 = coeffs[20]
    a73 = coeffs[21]
    a74 = coeffs[22]
    a75 = coeffs[23]
    a76 = coeffs[24]
    e1  = coeffs[25]
    e3  = coeffs[26]
    e4  = coeffs[27]
    e5  = coeffs[28]
    e6  = coeffs[29]
    e7  = coeffs[30]

    facold = 1e-4
    expo1 = 0.20-solver.options.beta*0.75
    facc1 = 1.0/solver.options.step_size_selection_one
    facc2 = 1.0/solver.options.step_size_selection_two
    posneg = sign(1.0, xend-solver.x)

    ###### Initial Preparations
    nfcn   = 0
    nstep  = 0
    naccpt = 0
    nrejct = 0

    atol = solver.options.atol
    rtol = solver.options.rtol
    atol_iter = atol isa Number ? repeated(atol) : atol
    rtol_iter = rtol isa Number ? repeated(rtol) : rtol

    last = false
    hlamb = 0.0
    iasti = 0
    nonsti = 0
    # fcn(n, x, y, k1)
    solver.f(solver.x, solver.y, solver.k1)
    hmax = abs(hmax)
    iord = 5

    if h == 0.0
        # solver contains, fcn, x, y, k1, k2, k3, atol, rtol
        h = hinit(solver, posneg, iord, hmax)
    end

    nfcn += 2
    reject = false
    
    ###### Basic Integration Step
    xph = 0.0 # need to put this here because FORTRAN seems to be looser with scoping
    while true
        if nstep > solver.options.maximum_allowed_steps
            # GOTO 78
            println(" MORE THAN NMAX = ", solver.options.maximum_allowed_steps, " STEPS ARE NEEDED")
            return h, DP5Report(solver.x, -2, 0, 0, 0, 0)
        end
        
        if (0.10 * abs(h)) <= abs(solver.x)*solver.options.uround 
            # GOTO 77
            println("STEP SIZE TOO SMALL, H = ", h)
            return h, DP5Report(solver.x, -3, 0, 0, 0, 0)
        end

        if ((solver.x+1.01*h-xend)*posneg) > 0.0
            h = xend-solver.x
            last = true
        end
        
        nstep += 1

        ####### First 6 stages, just set to equality and should work bc everything is vector (no need for loops)
        # 22
        copyto!(solver.y1, solver.y + h * a21 * solver.k1)
        solver.f(solver.x + c2 * h, solver.y1, solver.k2)
        # 23
        copyto!(solver.y1, solver.y + h * (a31 * solver.k1 + a32 * solver.k2))
        solver.f(solver.x + c3 * h, solver.y1, solver.k3)
        # 24
        copyto!(solver.y1, solver.y + h * (a41 * solver.k1 + a42 * solver.k2 + a43 * solver.k3))
        solver.f(solver.x + c4 * h, solver.y1, solver.k4)
        # 25
        copyto!(solver.y1, solver.y + h * (a51 * solver.k1 + a52 * solver.k2 + a53 * solver.k3 + a54 * solver.k4))
        solver.f(solver.x + c5*h, solver.y1, solver.k5)
        # 26
        copyto!(solver.ysti, solver.y + h * (a61 * solver.k1 + a62 * solver.k2 + a63 * solver.k3 + a64 * solver.k4 + a65 * solver.k5))
        xph = solver.x + h
        solver.f(xph, solver.ysti, solver.k6)
        # 27
        copyto!(solver.y1, solver.y + h * (a71 * solver.k1 + a73 * solver.k3 + a74 * solver.k4 + a75 * solver.k5 + a76 * solver.k6))
        solver.f(xph, solver.y1, solver.k2)
        # 28
        copyto!(solver.k4, h * (e1 * solver.k1 + e3 * solver.k3 + e4 * solver.k4 + e5 * solver.k5 + e6 * solver.k6 + e7 * solver.k2))

        nfcn += 6

        ###### Error Estimation
        err = 0.0

        err = mapreduce(+, atol_iter, rtol_iter, solver.k4, solver.y, solver.ysti) do atoli, rtoli, k4i, yi, ystii
            sk = atoli + rtoli*max(abs(yi), abs(ystii))
            abs(k4i/sk)^2
        end

        err = sqrt(err/length(solver.y))

        ###### Computation of hnew
        fac11 = err^expo1
        ###### Lund-Stabilization
        fac = fac11/(facold^solver.options.beta)
        ###### we require fac1 <= hnew/h <= fac2
        fac = max(facc2, min(facc1, fac/solver.options.safety_factor)) # facc1, facc2, fac must be Float64 
        hnew = h/fac
        if err <= 1.0 
            ###### Step is accepted
            facold = max(err, 1e-4)
            naccpt += 1
            ###### Stiffness Detection
            if (mod(naccpt, solver.options.stiffness_test_activation_step) == 0) || (iasti > 0)
                stnum = 0.0
                stden = 0.0

                stnum, stden = mapreduce(.+, solver.k2, solver.k6, solver.y1, solver.ysti) do k2i, k6i, y1i, ystii
                    stnum = abs(k2i-k6i)^2
                    stden = abs(y1i-ystii)^2
                    stnum, stden
                end

                if stden > 0.0
                    hlamb = h*sqrt(stnum/stden)
                end

                
                if hlamb > 3.25
                    iasti += 1
                    if iasti == 15
                        println("THE PROBLEM SEEMS TO BECOME STIFF AT X = ", x)    
                    end
                else 
                    nonsti += 1
                    if nonsti == 6
                        iasti = 0
                    end
                end
            end

            copyto!(solver.k1, solver.k2)
            copyto!(solver.y, solver.y1)

            solver.x = xph

            ###### Normal Exit
            if last 
                h = hnew
                return h, DP5Report(
                    solver.x, 
                    1, 
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
            hnew = h/min(facc1, fac11/solver.options.safety_factor)
            reject = true
            if naccpt > 1
                nrejct += 1
            end
            last = false

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
    dnf = 0.0
    dny = 0.0

    atol = solver.options.rtol
    rtol = solver.options.atol

    atol_iter = atol isa Number ? repeated(atol) : atol
    rtol_iter = rtol isa Number ? repeated(rtol) : rtol

    dnf, dny = mapreduce(.+, atol_iter, rtol_iter, solver.k1, solver.y) do atoli, rtoli, f0i, yi
        sk = atoli + rtoli*abs(yi)
        dnf = abs(f0i/sk)^2
        dny = abs(yi/sk)^2
        dnf, dny
    end

   
    # problem with comparing ComplexF64 with Float64
    # take abs of dnf 
    if (dnf <= 1.0e-10) || (dny <= 1.0e-10)
        h = 1.0e-6
    else
        h = 0.01*sqrt(dny/dnf)
    end
    h = min(h, hmax)
    h = sign(h, posneg)

    ###### Perform an explicit step
    #y1 = y + h*f0
    #fcn(n, x+h, y1, f1)
    copyto!(solver.y1, solver.y + h*solver.k1)
    solver.f(solver.x + h, solver.k3, solver.k2)

    ###### Estimate the second derivative of the solution
    der2 = 0.0
        
    der2 = mapreduce(+, atol_iter, rtol_iter, solver.k2, solver.k1, solver.y) do atoli, rtoli, f1i, f0i, yi
        sk = atoli + rtoli*abs(yi)
        ((f1i-f0i)/sk)^2
    end

    der2 = sqrt(der2)/h

    ##### Step size is computed such that
    ##### H**IORD * MAX ( NORM(F0), NORM(F1), DER2 ) = 0.01
    der12 = max(abs(der2), sqrt(dnf))
    if der12 <= 1e-15
        h1 = max(1.0e-6, abs(h)*1.0e-3)
    else
        h1 = (0.01/der12)^(1.0/iord)
    end
    h = min(100*abs(h), h1, hmax)
    return sign(h, posneg)

end

function cdopri()

    return [
        0.2, #C2
        0.3, #C3
        0.8, #C4
        8.0/9.0, #C5
        0.2, #A21 
        3.0/40.0, #A31
        9.0/40.0, #A32
        44.0/45.0, #A4
        -56.0/15.0, #A42
        32.0/9.0, #A43 
        19372.0/6561.0, #A52
        -25360.0/2187.0, #
        64448.0/6561.0,
        -212.0/729.0,
        9017.0/3168.0,
        -355.0/33.0,
        46732.0/5247.0,
        49.0/176.0,
        -5103.0/18656.0,
        35.0/384.0,
        500.0/1113.0,
        125.0/192.0,
        -2187.0/6784.0,
        11.0/84.0,
        71.0/57600.0,
        -71.0/16695.0,
        71.0/1920.0,
        -17253.0/339200.0,
        22.0/525.0,
        -1.0/40.0
    ]
end

# FORTRAN sign - returns the value of A with the sign of B
function sign(a,b)
    if b >= 0
        sign = 1.0
    else
        sign = -1.0
    end
    
    return a*sign
end