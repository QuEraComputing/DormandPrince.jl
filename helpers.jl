# FORTRAN sign - returns the value of A with the sign of B
function sign(a,b)
    if b >= 0
        sign = 1.0
    else
        sign = -1.0
    end
    
    return a*sign
end

function do_step!(solver, h)

    # define constants
    c2=0.2
    c3=0.3
    c4=0.8
    c5=8.0/9.0
    a21=0.2
    a31=3.0/40.0
    a32=9.0/40.0
    a41=44.0/45.0
    a42=-56.0/15.0
    a43=32.0/9.0
    a51=19372.0/6561.0
    a52=-25360.0/2187.0
    a53=64448.0/6561.0
    a54=-212.0/729.0
    a61=9017.0/3168.0
    a62=-355.0/33.0
    a63=46732.0/5247.0
    a64=49.0/176.0
    a65=-5103.0/18656.0
    a71=35.0/384.0
    a73=500.0/1113.0
    a74=125.0/192.0
    a75=-2187.0/6784.0
    a76=11.0/84.0
    e1=71.0/57600.0
    e3=-71.0/16695.0
    e4=71.0/1920.0
    e5=-17253.0/339200.0
    e6=22.0/525.0
    e7=-1.0/40.0

    ####### First 6 stages, just set to equality and should work bc everything is vector (no need for loops)
    # 22
    solver.y1 .= solver.y + h * a21 * solver.k1
    solver.f(solver.x + c2 * h, solver.y1, solver.k2)
    # 23
    solver.y1 .= solver.y + h * (a31 * solver.k1 + a32 * solver.k2)
    solver.f(solver.x + c3 * h, solver.y1, solver.k3)
    # 24
    solver.y1 .= solver.y + h * (a41 * solver.k1 + a42 * solver.k2 + a43 * solver.k3)
    solver.f(solver.x + c4 * h, solver.y1, solver.k4)
    # 25
    solver.y1 .= solver.y + h * (a51 * solver.k1 + a52 * solver.k2 + a53 * solver.k3 + a54 * solver.k4)
    solver.f(solver.x + c5*h, solver.y1, solver.k5)
    # 26
    solver.ysti .= solver.y + h * (a61 * solver.k1 + a62 * solver.k2 + a63 * solver.k3 + a64 * solver.k4 + a65 * solver.k5)
    xph = solver.x + h
    solver.f(xph, solver.ysti, solver.k6)
    # 27
    solver.y1 .= solver.y + h * (a71 * solver.k1 + a73 * solver.k3 + a74 * solver.k4 + a75 * solver.k5 + a76 * solver.k6)
    solver.f(xph, solver.y1, solver.k2)
    # 28
    solver.k4 .= h * (e1 * solver.k1 + e3 * solver.k3 + e4 * solver.k4 + e5 * solver.k5 + e6 * solver.k6 + e7 * solver.k2)

end

function error_estimation(solver)

    err = mapreduce(+, solver.consts.atol_iter, solver.consts.rtol_iter, solver.k4, solver.y, solver.ysti) do atoli, rtoli, k4i, yi, ystii
        sk = atoli + rtoli*max(abs(yi), abs(ystii))
        abs(k4i/sk)^2
    end

    err = sqrt(err/length(solver.y))

    return err 
end

function estimate_second_derivative(solver, h)
        
    der2 = mapreduce(+, solver.consts.atol_iter, solver.consts.rtol_iter, solver.k2, solver.k1, solver.y) do atoli, rtoli, f1i, f0i, yi
        sk = atoli + rtoli*abs(yi)
        ((f1i-f0i)/sk)^2
    end

    der2 = sqrt(der2)/h

    return der2

end

function stiffness_detection!(solver, naccpt, h)
    if (mod(naccpt, solver.options.stiffness_test_activation_step) == 0) || (solver.iasti > 0)
        #stnum = 0.0
        #stden = 0.0

        stnum, stden = mapreduce(.+, solver.k2, solver.k6, solver.y1, solver.ysti) do k2i, k6i, y1i, ystii
            #stnum = abs(k2i-k6i)^2 
            #stden = abs(y1i-ystii)^2
            abs(k2i-k6i)^2, abs(y1i-ystii)^2
            # stnum, stden
        end

        if stden > 0.0
            solver.hlamb = h*sqrt(stnum/stden)
        else
            solver.hlamb = Inf
        end

        
        if solver.hlamb > 3.25
            solver.iasti += 1
            if solver.iasti == 15
                # turn this into a debug statement
                @debug "The problem seems to become stiff at $x" 
            end
        else 
            solver.nonsti += 1
            if solver.nonsti == 6
                solver.iasti = 0
            end
        end
    end
end

function euler_first_guess(solver, hmax, posneg)

    dnf, dny = mapreduce(.+, solver.consts.atol_iter, solver.consts.rtol_iter, solver.k1, solver.y) do atoli, rtoli, f0i, yi
        sk = atoli + rtoli*abs(yi)
        abs(f0i/sk)^2, abs(yi/sk)^2 # dnf, dny
    end

   
    # problem with comparing ComplexF64 with Float64
    # take abs of dnf 
    if (dnf <= 1.0e-10) || (dny <= 1.0e-10)
        h = 1.0e-6
    else
        h = 0.01*sqrt(dny/dnf)
    end
    h = min(h, hmax)
    h = h * Base.sign(posneg)

    return h, dnf
end