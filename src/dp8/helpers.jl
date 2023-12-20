function do_step!(solver, h)


    ####### First 6 stages
    # DO 22 I=1,N 
    # 22  Y1(I)=Y(I)+H*A21*K1(I)  
    #     CALL FCN(N,X+C2*H,Y1,K2,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* a21 .* solver.k1
    solver.f(solver.vars.x + c2 * h, solver.y1, solver.k2)
    #     DO 23 I=1,N 
    # 23  Y1(I)=Y(I)+H*(A31*K1(I)+A32*K2(I))  
    #     CALL FCN(N,X+C3*H,Y1,K3,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (a31 .* solver.k1 .+ a32 .* solver.k2)
    solver.f(solver.vars.x + c3 * h, solver.y1, solver.k3)
    #     DO 24 I=1,N 
    # 24  Y1(I)=Y(I)+H*(A41*K1(I)+A43*K3(I))  
    #     CALL FCN(N,X+C4*H,Y1,K4,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (a41 .* solver.k1 .+ a43 .* solver.k3)
    solver.f(solver.vars.x + c4 * h, solver.y1, solver.k4)
    #     DO 25 I=1,N 
    # 25  Y1(I)=Y(I)+H*(A51*K1(I)+A53*K3(I)+A54*K4(I))
    #     CALL FCN(N,X+C5*H,Y1,K5,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (a51 .* solver.k1 .+ a53 .* solver.k3 .+ a54 .* solver.k4)
    solver.f(solver.vars.x + c5*h, solver.y1, solver.k5)
    #     DO 26 I=1,N 
    # 26  Y1(I)=Y(I)+H*(A61*K1(I)+A64*K4(I)+A65*K5(I))
    #     CALL FCN(N,X+C6*H,Y1,K6,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (a61 .* solver.k1 .+ a64 .* solver.k4 .+ a65 .* solver.k5)
    solver.f(solver.vars.x + c6*h, solver.y1, solver.k6)
    #     DO 27 I=1,N 
    # 27  Y1(I)=Y(I)+H*(A71*K1(I)+A74*K4(I)+A75*K5(I)+A76*K6(I))
    #     CALL FCN(N,X+C7*H,Y1,K7,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (
        a71 .* solver.k1 .+ a74 .* solver.k4 .+ a75 .* solver.k5 .+ a76 .* solver.k6
    )
    solver.f(solver.vars.x + c7*h, solver.y1, solver.k7)
    #     DO 28 I=1,N 
    # 28  Y1(I)=Y(I)+H*(A81*K1(I)+A84*K4(I)+A85*K5(I)+A86*K6(I)+A87*K7(I))  
    #     CALL FCN(N,X+C8*H,Y1,K8,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (
        a81 .* solver.k1 .+ a84 .* solver.k4 .+ a85 .* solver.k5 .+ a86 .* solver.k6 .+ a87 .* solver.k7
    )
    solver.f(solver.vars.x + c8*h, solver.y1, solver.k8)
    #     DO 29 I=1,N 
    # 29  Y1(I)=Y(I)+H*(A91*K1(I)+A94*K4(I)+A95*K5(I)+A96*K6(I)+A97*K7(I)
    #    &   +A98*K8(I))
    #     CALL FCN(N,X+C9*H,Y1,K9,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (
        a91 .* solver.k1 .+ a94 .* solver.k4 .+ a95 .* solver.k5 .+ 
        a96 .* solver.k6 .+ a97 .* solver.k7 .+ a98 .* solver.k8
    )
    solver.f(solver.vars.x + c9*h, solver.y1, solver.k9)
    #     DO 30 I=1,N 
    # 30  Y1(I)=Y(I)+H*(A101*K1(I)+A104*K4(I)+A105*K5(I)+A106*K6(I)
    #    &   +A107*K7(I)+A108*K8(I)+A109*K9(I))
    #     CALL FCN(N,X+C10*H,Y1,K10,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (
        a101 .* solver.k1 .+ a104 .* solver.k4 .+ a105 .* solver.k5 .+ 
        a106 .* solver.k6 .+ a107 .* solver.k7 .+ a108 .* solver.k8 .+ a109 .* solver.k9
    )
    solver.f(solver.vars.x + c10*h, solver.y1, solver.k10)
    #     DO 31 I=1,N 
    # 31  Y1(I)=Y(I)+H*(A111*K1(I)+A114*K4(I)+A115*K5(I)+A116*K6(I)
    #    &   +A117*K7(I)+A118*K8(I)+A119*K9(I)+A1110*K10(I))
    #     CALL FCN(N,X+C11*H,Y1,K2,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (
        a111 .* solver.k1 .+ a114 .* solver.k4 .+ a115 .* solver.k5 .+ 
        a116 .* solver.k6 .+ a117 .* solver.k7 .+ a118 .* solver.k8 .+ 
        a119 .* solver.k9 .+ a1110 .* solver.k10
    )
    solver.f(solver.vars.x + c11*h, solver.y1, solver.k2)
    #     XPH=X+H
    xph = solver.vars.x + h
    #     DO 32 I=1,N 
    # 32  Y1(I)=Y(I)+H*(A121*K1(I)+A124*K4(I)+A125*K5(I)+A126*K6(I)
    #    &   +A127*K7(I)+A128*K8(I)+A129*K9(I)+A1210*K10(I)+A1211*K2(I))
    #     CALL FCN(N,XPH,Y1,K3,RPAR,IPAR)
    solver.y1 .= solver.y .+ h .* (
        a121 .* solver.k1 .+ a124 .* solver.k4 .+ a125 .* solver.k5 .+ 
        a126 .* solver.k6 .+ a127 .* solver.k7 .+ a128 .* solver.k8 .+ 
        a129 .* solver.k9 .+ a1210 .* solver.k10 .+ a1211 .* solver.k2
    )
    solver.f(xph, solver.y1, solver.k3)
    #     NFCN=NFCN+11
    #     DO 35 I=1,N 
    #     K4(I)=B1*K1(I)+B6*K6(I)+B7*K7(I)+B8*K8(I)+B9*K9(I)
    #    &   +B10*K10(I)+B11*K2(I)+B12*K3(I)   
    #     K5(I)=Y(I)+H*K4(I)
    solver.k4 .= (
        b1 .* solver.k1 .+ b6 .* solver.k6 .+ b7 .* solver.k7 .+ b8 .* solver.k8 .+ 
        b9 .* solver.k9 .+ b10 .* solver.k10 .+ b11 .* solver.k2 .+ b12 .* solver.k3
    )
    solver.k5 .= solver.y .+ h .* solver.k4

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
    if (mod(naccpt, solver.options.stiffness_test_activation_step) == 0) || (solver.vars.iasti > 0)
        #stnum = 0.0
        #stden = 0.0

        stnum, stden = mapreduce(.+, solver.k2, solver.k6, solver.y1, solver.ysti) do k2i, k6i, y1i, ystii
            #stnum = abs(k2i-k6i)^2 
            #stden = abs(y1i-ystii)^2
            abs(k2i-k6i)^2, abs(y1i-ystii)^2
            # stnum, stden
        end

        if stden > 0.0
            solver.vars.hlamb = h*sqrt(stnum/stden)
        else
            solver.vars.hlamb = Inf
        end

        
        if solver.vars.hlamb > 3.25
            solver.vars.iasti += 1
            if solver.vars.iasti == 15
                @debug "The problem seems to become stiff at $x" 
            end
        else 
            solver.vars.nonsti += 1
            if solver.vars.nonsti == 6
                solver.vars.iasti = 0
            end
        end
    end
end

function euler_first_guess(solver, hmax, posneg)

    dnf, dny = mapreduce(.+, solver.consts.atol_iter, solver.consts.rtol_iter, solver.k1, solver.y) do atoli, rtoli, f0i, yi
        sk = atoli + rtoli*abs(yi)
        abs(f0i/sk)^2, abs(yi/sk)^2 # dnf, dny
    end


    if (dnf <= 1.0e-10) || (dny <= 1.0e-10)
        h = 1.0e-6
    else
        h = 0.01*sqrt(dny/dnf)
    end
    h = min(h, hmax)
    h = h * Base.sign(posneg)

    return h, dnf
end