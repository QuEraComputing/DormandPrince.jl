# returns (explicitly)
# IDID for error value, 
# computed x (should be the x_end value if successful)
# 

function dopri5(
    n,
    fcn, 
    x, # mutate
    y, # mutate
    xend, 
    rtol,
    atol, 
    itol,
    work,  
    iwork, # mutate
)

    arret = false

    ###### nmax - maximal number of steps
    if iwork[1] == 0 
        nmax = 100000
    else
        nmax = iwork[1]
        if nmax < 0
            println("WRONG INPUT IWORK[1]=", iwork[1])
            arret = true
        end
    end

    ###### nstiff -  parameters for stiffness detection
    nstiff = iwork[4]
    if nstiff == 0
        nstiff = 1000
    end
    if nstiff < 0
        nstiff = nmax + 10
    end

    ###### nrdens - number of dense output components
    nrdens = iwork[5]
    if (nrdens < 0) || (nrdens > n) 
        println("CURIOUS INPUT IWORK[5]=", iwork[5])
        arret = true
    else
        if nrdens == n
            for i in range(1, nrdens)
                iwork[20 + i] = i
            end
        end
    end

    ###### uround - smallest number satisfying 1.0 + uround > 1.0
    
    if work[1] == 0.0
        uround = 2.3e-16
    else
        uround = work[1]
        if (uround <= 1e-35) || (uround >= 1.0)
            println("WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:", work[1])
            arret = true
        end
    end

    ####### safety factor
    
    if work[2] == 0.0 
        safe = 0.9
    else
        safe = work[2]
        if (safe >= 1.0) || (safe <= 1e-4)
            println("CURIOUS INPUT FOR SAFETY FACTOR WORK[2]=", work[2])
            arret = true
        end
    end

    ###### fac1, fac2 - parameters for step size selection
    if work[3] == 0.0
        fac1 = 0.2
    else
        fac1 = work[3]
    end

    if work[4] == 0.0
        fac2 = 10.0
    else
        fac2 = work[4]
    end

    ###### beta - for step control stabilization
    if work[5] == 0.0
        beta = 0.04
    else
        if work[5] < 0.0
            beta = 0.0
        else 
            beta = work[5]
            if beta > 0.2
                println("CURIOUS INPUT FOR BETA: WORK[5]=", work[5])
                arret = true
            end
        end
    end

    ####### maximal step size
    if work[6] == 0.0
        hmax = xend-x
    else
        hmax = work[6]
    end

    ####### initial step size
    h = work[7]

    ####### prepare entry-points for arrays in work
    iey1 = 21
    iek1 = iey1 + n
    iek2 = iek1 + n
    iek3 = iek2 + n
    iek4 = iek3 + n
    iek5 = iek4 + n
    iek6 = iek5 + n
    ieys = iek6 + n
    # ieco = ieys + n

    ####### total storage requirement
    istore = ieys+5*nrdens-1
    if istore > length(work)
        println("INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", istore)
        arret = true
    end

    icomp = 21 
    istore = icomp + nrdens - 1
    if istore > length(iwork)
        println("INSUFFICIENT STORAGE FOR IWORK, MIN. LIWORK=", istore)
        arret = true
    end

    ###### When a fail has occurred, return with idid = -1
    if arret
        # (idid, x)
        return (-1, x)
    end

    # indices to work start at the starting locations but for a view we need 
    x, idid = dopcor(
        n, fcn, x, y, xend, hmax, h, rtol, atol, itol, nmax, uround, nstiff,
        safe, beta, fac1, fac2, 
        view(work, iey1:iey1+n-1), 
        view(work, iek1:iek1+n-1), 
        view(work, iek2:iek2+n-1), 
        view(work, iek3:iek3+n-1),
        view(work, iek4:iek4+n-1), 
        view(work, iek5:iek5+n-1), 
        view(work, iek6:iek6+n-1),
        view(work, ieys:ieys+n-1),
        nrdens
    )
    
    

end

# return idid, x
function dopcor(
    n,
    fcn, 
    x, 
    y, 
    xend,
    hmax,
    h, 
    rtol, 
    atol, 
    itol, 
    nmax,
    uround, 
    nstiff, 
    safe, 
    beta, 
    fac1, 
    fac2, 
    y1, 
    k1, 
    k2, 
    k3, 
    k4, 
    k5, 
    k6, 
    ysti,
    nrd
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
    expo1 = 0.20-beta*0.75
    facc1 = 1.0/fac1
    facc2 = 1.0/fac2
    posneg = sign(1.0, xend-x)

    ###### Initial Preparations
    nfcn   = 0
    nstep  = 0
    naccpt = 0
    nrejct = 0

    atoli = atol[1] # works with integers or arrays
    rtoli = rtol[1] # works with integers or arrays
    last = false
    hlamb = 0.0
    iasti = 0
    fcn(n, x, y, k1)
    hmax = abs(hmax)
    iord = 5
    
    if h == 0.0
        h = hinit(n, fcn, x, y, xend, posneg, k1, k2, k3, iord, hmax, atol, rtol, itol)
    end

    nfcn += 2
    reject = false


    ###### Basic Integration Step
    xph = 0.0 # need to put this here because FORTRAN seems to be looser with scoping
    while true
        if nstep > nmax
            # GOTO 78
            println(" MORE THAN NMAX = ", nmax, " STEPS ARE NEEDED")
            return (-2, x)
        end
        
        if (0.10 * abs(h)) <= abs(x)*uround 
            # GOTO 77
            println("STEP SIZE TOO SMALL, H = ", h)
            return (-3, x)
        end

        if ((x+1.01*h-xend)*posneg) > 0.0
            h = xend-x
            last = true
        end
        
        nstep += 1

        ####### First 6 stages
        # 22
        for i in range(1, n)
            y1[i] = y[i] + h*a21*k1[i]
            fcn(n, x+c2*h, y1, k2)
        end
        # 23
        for i in range(1, n)
            y1[i] = y[i] + h*(a31*k1[i]+a32*k2[i])
            fcn(n, x+c3*h, y1, k3)
        end
        # 24
        for i in range(1, n)
            y1[i] = y[i] + h*(a41*k1[i]+a42*k2[i]+a43*k3[i])
            fcn(n, x+c4*h, y1, k4)
        end
        # 25
        for i in range(1, n)
            y1[i] = y[i] + h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
            fcn(n, x+c5*h, y1, k5)
        end
        # 26
        for i in range(1, n)
            ysti[i] = y[i] + h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
            xph = x + h
            fcn(n, xph, ysti, k6)
        end
        # 27
        for i in range(1, n)
            y1[i] = y[i] + h*(a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
            fcn(n, xph, y1, k2)
        end
        # 28
        for i in range(1, n)
            k4[i] = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k2[i])
        end

        nfcn += 6

        ###### Error Estimation
        err = 0.0
        if itol == 0 
            for i in range(1, n)
                sk = atoli + rtoli*max(abs(y[i]), abs(y1[i]))
                err += (k4[i]/sk)^2
            end
        else
            for i in range(1, n)
                sk = atoli[i] + rtoli[i]*max(abs(y[i]), abs(y1[i]))
                err += (k4[i]/sk)^2
            end
        end
        err = sqrt(err/n)

        ###### Computation of hnew
        fac11 = err^expo1
        ###### Lund-Stabilization
        fac = fac11/(facold^beta)
        ###### we require fac1 <= hnew/h <= fac2
        fac = max(facc2, min(facc1, fac/safe))
        hnew = h/fac
        if err <= 1.0 
            ###### Step is accepted
            facold = max(err, 1e-4)
            naccpt += 1
            ###### Stiffness Detection
            if (mod(naccpt, nstiff) == 0) || (iasti > 0)
                stnum = 0.0
                stden = 0.0
                for i in range(1, n)
                    stnum += (k2[i]-k6[i])^2
                    stden += (y1[i]-ysti[i])^2
                end

                if stden > 0.0
                    hlamb = h*sqrt(stnum/stden)
                end

                if hlamb > 3.25
                    nonsti = 0
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

            # 44
            for i in range(1, n)
                k1[i] = k2[i]
                y[i] = y1[i]
            end

            x = xph

            ###### Normal Exit
            if last 
                h = hnew
                iwork[17] = nfcn 
                iwork[18] = nstep
                iwork[19] = naccpt
                iwork[20] = nrejct
                return (1, x)
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
            hnew = h/min(facc1, fac11/safe)
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
    n, 
    fcn, 
    x, 
    y, 
    xend, 
    posneg, 
    f0, 
    f1, 
    y1, 
    iord,
    hmax, 
    atol, 
    rtol,
    itol
)
    #=
    Compute a first guess for explicit euler as
        h = 0.01 * norm (y0) / norm (f0)
    the increment for explicit euler is small
    compared to the solution
    =#
    dnf = 0.0
    dny = 0.0
    atoli = atol[1]
    rtoli = rtol[1]
    
    if itol == 0
        for i in range(1, n)
            sk = atoli + rtoli*abs(y[i])
            dnf += (f0[i]/sk)^2
            dny += (y[i]/sk)^2
        end
    else 
        for i in range(1, n)
            sk = atol[i] + rtol[i]*abs(y[i])
            dnf += (f0[i]/sk)^2
            dny += (y[i]/sk)^2
        end
    end
    
    if (dnf <= 1.0e-10) || (dny <= 1.0e-10)
        h = 1.0e-6
    else
        h = 0.01*sqrt(dny/dnf)
    end
    h = min(h, hmax)
    h = sign(h, posneg)

    ###### Perform an explicit step
    for i in range(1, n)
        y1[i] = y[i] + h*f0[i]
    end 
    fcn(n, x+h, y1, f1)
    ###### Estimate the second derivative of the solution
    der2 = 0.0
    if itol == 0 
        for i in range(1, n)
            sk = atoli + rtoli*abs(y[i])
            der2 += ((f1[i]-f0[i])/sk)^2
        end
    else
        for i in range(1, n)
            sk = atoli[i] + rtoli[i]*abs(y[i])
            der2 += ((f1[i]-f0[i])/sk)^2
        end
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

function debug_dopcor(
    n,
    fcn, 
    x, 
    y, 
    xend,
    hmax,
    h, 
    rtol, 
    atol, 
    itol, 
    nmax,
    uround, 
    nstiff, 
    safe, 
    beta, 
    fac1, 
    fac2, 
    y1, 
    k1, 
    k2, 
    k3, 
    k4, 
    k5, 
    k6, 
    ysti,
    nrd,  
    rpar,
    ipar
)
    println("n ", n)
    println("fcn ", fcn)
    println("x ", x)
    println("y ", y)
    println("xend ", xend)
    println("hmax ", hmax)
    println("h ", h)
    println("rtol ", rtol)
    println("atol ", atol)
    println("itol ", itol)
    println("nmax ", nmax)
    println("uround ", uround)
    println("nstiff ", nstiff)
    println("safe ",safe)
    println("beta ", beta)
    println("fac1 ", fac1)
    println("fac2 ", fac2)
    println("y1 ", y1)
    println("k1 ", k1)
    println("k2 ", k2)
    println("k3 ", k3)
    println("k4 ", k4)
    println("k5 ", k5)
    println("k6 ", k6)
    println("ysti ", ysti)
    println("nrd ",nrd)
    println("rpar ", rpar)
    println("ipar ", ipar)


end

# FORTRAN sign - returns the value of a with the sign of B
function sign(a,b)
    if b >= 0
        sign = 1.0
    else
        sign = -1.0
    end
    
    return a*sign
end