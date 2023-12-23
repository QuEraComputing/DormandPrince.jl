using Test
using DormandPrince: DP5Solver, DP8Solver, integrate

function stiff_fcn(x, y, f) 
    f[1] = y[1]^2 - y[1]^3
end

@testset "Stiff ODE" begin

    for SolverType in [DP5Solver, DP8Solver]
        solver = SolverType(
            stiff_fcn,
            0.0,
            [0.0001]
        )

        integrate(solver, 2/0.0001)
    end    
end
