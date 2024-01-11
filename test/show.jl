# show solver and report

using DormandPrince
using Test

function stiff_fcn(x, y, f) 
    f[1] = y[1]^2 - y[1]^3
end

@testset "DP8 Solver" begin
    solver = DP8Solver(
        stiff_fcn,
        0.0,
        [0.0001]
    )

    show(stdout, MIME"text/plain"(), solver)
end

@testset "Successful Options Check with Successful Integration" begin
    solver = DP5Solver(
        stiff_fcn,
        0.0,
        [0.0001]
    )

    show(stdout, MIME"text/plain"(), solver)

    report = DormandPrince.integrate_core!(solver, 2/0.0001)

    show(stdout, MIME"text/plain"(), report)
    
end

@testset "Successful Options Check with Failed Integration" begin
    solver = DP5Solver(
        stiff_fcn,
        0.0,
        [0.0001];
        maximum_allowed_steps=10
    )

    report = DormandPrince.integrate_core!(solver, 2/0.0001)

    show(stdout, MIME"text/plain"(), report)
end

@testset "Failed Options Check" begin
    solver = DP5Solver(
        stiff_fcn,
        0.0,
        [0.0001];
        uround = 100
    )

    report = DormandPrince.integrate_core!(solver, 2/0.0001)

    show(stdout, MIME"text/plain"(), report)

end