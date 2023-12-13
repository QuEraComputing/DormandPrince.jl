include("../types.jl")
include("../solver.jl")

function fcn(x, y, f)
    f[1] = -1im * y[2]
    f[2] = -1im * y[1]
end

solver = DP5Solver(
    fcn,
    0.0,
    ComplexF64[1.0, 0.0]
)

dopri5(solver, 2π)