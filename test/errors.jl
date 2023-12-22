using Test
using DormandPrince:
    DP5Solver,
    Options,
    LARGER_NMAX_NEEDED,
    STEP_SIZE_BECOMES_TOO_SMALL
    
using DormandPrince. DP5: dopri5

@testset "Larger nmax needed" begin
    solver = DP5Solver(
        stiff_fcn,
        0.0, # start at 0.0
        [0.0001] # delta
        ; maximum_allowed_steps=1
    )

    
    h, report = dopri5(solver, 2/0.0001, 0.1, 0.0)
    @test report.idid == LARGER_NMAX_NEEDED

end

@testset "Step size becomes too small" begin
    solver = DP5Solver(
        stiff_fcn,
        0.0,
        [0.0001]
        ; uround=10000
    )

    
    h, report = dopri5(solver, 2/0.0001, 0.1, 0.0)
    @test report.idid == STEP_SIZE_BECOMES_TOO_SMALL
end
