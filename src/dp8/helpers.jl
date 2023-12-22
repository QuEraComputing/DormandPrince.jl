function do_step!(solver::DP8Solver{T}, h::T) where T

    ####### The 12 stages
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

function error_estimation(solver::DP8Solver{T}, h::T) where T

    err1, err2 = mapreduce(.+, 
        solver.consts.atol_iter, 
        solver.consts.rtol_iter, 
        solver.y,
        solver.k1,
        solver.k2,
        solver.k3,
        solver.k4,
        solver.k5,
        solver.k6,
        solver.k7,
        solver.k8,
        solver.k9,
        solver.k10
        ; init=(zero(T), zero(T))
    ) do atoli, rtoli, yi, k1i, k2i, k3i, k4i, k5i, k6i, k7i, k8i, k9i, k10i

        #     DO 42 I=1,N 
        #     SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(K5(I)))
        #     ERRI=K4(I)-BHH1*K1(I)-BHH2*K9(I)-BHH3*K3(I)
        #     ERR2=ERR2+(ERRI/SK)**2
        #     ERRI=ER1*K1(I)+ER6*K6(I)+ER7*K7(I)+ER8*K8(I)+ER9*K9(I)
        # &      +ER10*K10(I)+ER11*K2(I)+ER12*K3(I)
        # 42    ERR=ERR+(ERRI/SK)**2

        sk = atoli + rtoli*max(abs(yi), abs(k5i))
        erri2 = k4i - bhh1*k1i - bhh2*k9i - bhh3*k3i
        erri1 = er1*k1i + er6*k6i + er7*k7i + er8*k8i + er9*k9i + er10*k10i + er11*k2i + er12*k3i
        
        (abs(erri1)/sk)^2, (abs(erri2)/sk)^2
    end

    # DENO=ERR+0.01D0*ERR2
    # IF (DENO.LE.0.D0) DENO=1.D0
    # ERR=ABS(H)*ERR*SQRT(1.D0/(N*DENO))
    den0 = err1 + 0.01*err2
    den0 = den0 <= 0.0 ? 1.0 : den0
    err = abs(h)*err1*sqrt(1.0/(length(solver.y) * den0))
    return err 
end

function stiffness_detection!(solver::DP8Solver{T}, naccpt::Int, h::T) where T

    if (mod(naccpt, solver.options.stiffness_test_activation_step) == 0) || (solver.vars.iasti > 0)
        #stnum = 0.0
        #stden = 0.0

        stnum, stden = mapreduce(.+, solver.k3, solver.k4, solver.k5, solver.y1) do k3i, k4i, k5i, y1i
            abs(k4i-k3i)^2, abs(k5i-y1i)^2
        end

        if stden > 0.0
            solver.vars.hlamb = abs(h)*sqrt(stnum/stden)
        else
            solver.vars.hlamb = Inf
        end

        
        if solver.vars.hlamb > 6.1
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
