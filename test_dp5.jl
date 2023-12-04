include("dp5.jl")
include("dp5_types.jl")

# FCN must follow format FCN(N, X, Y, F, RPAR, IPAR)

function fcn(n, x, y, f)
    f[1] = ℯ^x
    f[2] = ℯ^x
end

y = [1.0, 1.0]
work = zeros(16)

dopri5(
    2,
    fcn,
    0.0,
    y,
    1.0,
    1e-10,
    1e-10,
    0,
    DP5Options(),
    work,# work, should be 8*N
)