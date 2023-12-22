using Test
using DormandPrince: DP5Solver, integrate

function stiff_fcn(x, y, f) 
    f[1] = y[1]^2 - y[1]^3
end

@testset "Stiff ODE" begin

    solver = DP5Solver(
        stiff_fcn,
        0.0, # start at 0.0
        [0.0001] # initial value of delta
    )
    
    integrate(solver, 2/0.0001)
    
end
