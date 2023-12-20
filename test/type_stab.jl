using DormandPrince: DP5Solver, integrate
using JET: @report_opt


function fcn(x, y, f)
    g(x) = 2.2*2π*sin(2π*x)

    f[1] = -1im * g(x)/2 * y[2]
    f[2] = -1im * g(x)/2 * y[1]
end

solver = DP5Solver(
    fcn,
    0.0,
    ComplexF64[1.0, 0.0]
)

@report_opt integrate(solver, 2π)


