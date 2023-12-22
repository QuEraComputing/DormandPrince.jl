

function euler_first_guess(solver::AbstractDPSolver{T}, hmax::T, posneg::T) where T

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


function estimate_second_derivative(solver::AbstractDPSolver{T}, h::T) where {T}
        
    der2 = mapreduce(+, solver.consts.atol_iter, solver.consts.rtol_iter, solver.k2, solver.k1, solver.y) do atoli, rtoli, f1i, f0i, yi
        sk = atoli + rtoli*abs(yi)
        ((f1i-f0i)/sk)^2
    end

    der2 = sqrt(der2)/h

    return der2

end

function hinit(
    solver::AbstractDPSolver{T}, 
    posneg::T, 
    iord::Int,
    hmax::T
    # f0 arg is k1 from dopcor
    # f1 arg is k2 from dopcor
    # y1 arg is k3 from dopcor
) where T
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
    # copyto!(solver.y1, solver.y + h*solver.k1)
    solver.y1 .= solver.y .+ h .*solver.k1
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