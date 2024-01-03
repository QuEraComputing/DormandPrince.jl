using DormandPrince
using DormandPrince.DP8: dop853, error_estimation
using JET: @report_opt


function fcn(x, y, f)
    g(x) = 2.2*2π*sin(2π*x)

    f[1] = -1im * g(x)/2 * y[2]
    f[2] = -1im * g(x)/2 * y[1]
end

solver = DP8Solver(
    fcn,
    0.0,
    ComplexF64[1.0, 0.0]
)

h = 1e-6
# @report_opt dop853(solver, 1.0, 1.0, h)
# @code_warntype error_estimation(solver, 1e-6)
# @report_opt error_estimation(solver, 1e-6)

@code_warntype integrate!(solver, 2π)
@report_opt integrate!(solver, 2π)


