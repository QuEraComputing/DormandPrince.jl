include("dp5.jl")

# FCN must follow format FCN(N, X, Y, F, RPAR, IPAR)

function fcn(n, x, y, f)
    f[1] = ℯ^x
end

y = [1.0]
work = zeros(8)

dopri5(
    1,
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