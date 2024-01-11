# show solver and report

using DormandPrince

function stiff_fcn(x, y, f) 
    f[1] = y[1]^2 - y[1]^3
end

solver = DP5Solver(
    stiff_fcn,
    0.0,
    [0.0001]
    ;
    maximum_allowed_steps=100
)

show(solver)

report = DormandPrince.integrate_core!(solver, 2/0.0001)

show(report)