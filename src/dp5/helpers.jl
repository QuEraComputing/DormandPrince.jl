function do_step!(solver::DP5Solver{T}, h::T) where T

    ####### First 6 stages

    solver.y1 .= solver.y .+ h .* a21 .* solver.k1
    solver.f(solver.vars.x + c2 * h, solver.y1, solver.k2)

    solver.y1 .= solver.y .+ h .* (a31 .* solver.k1 .+ a32 .* solver.k2)
    solver.f(solver.vars.x + c3 * h, solver.y1, solver.k3)

    solver.y1 .= solver.y .+ h .* (a41 .* solver.k1 .+ a42 .* solver.k2 .+ a43 .* solver.k3)
    solver.f(solver.vars.x + c4 * h, solver.y1, solver.k4)

    solver.y1 .= solver.y .+ h .* (a51 .* solver.k1 .+ a52 .* solver.k2 .+ a53 .* solver.k3 .+ a54 .* solver.k4)
    solver.f(solver.vars.x + c5*h, solver.y1, solver.k5)
    
    solver.ysti .= solver.y .+ h .* (a61 .* solver.k1 .+ a62 .* solver.k2 .+ a63 .* solver.k3 .+ a64 .* solver.k4 .+ a65 .* solver.k5)
    xph = solver.vars.x + h
    solver.f(xph, solver.ysti, solver.k6)

    solver.y1 .= solver.y .+ h .* (a71 .* solver.k1 .+ a73 .* solver.k3 .+ a74 .* solver.k4 .+ a75 .* solver.k5 .+ a76 .* solver.k6)
    solver.f(xph, solver.y1, solver.k2)

    solver.k4 .= h .* (e1 .* solver.k1 .+ e3 .* solver.k3 .+ e4 .* solver.k4 .+ e5 .* solver.k5 .+ e6 .* solver.k6 .+ e7 .* solver.k2)

end

function error_estimation(solver)

    err = mapreduce(+, solver.consts.atol_iter, solver.consts.rtol_iter, solver.k4, solver.y, solver.ysti) do atoli, rtoli, k4i, yi, ystii
        sk = atoli + rtoli*max(abs(yi), abs(ystii))
        (abs(k4i)/sk)^2
    end

    err = sqrt(err/length(solver.y))

    return err 
end

function stiffness_detection!(solver::DP5Solver{T}, naccpt::Int, h::T) where T
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
            solver.vars.hlamb = abs(h)*sqrt(stnum/stden)
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
