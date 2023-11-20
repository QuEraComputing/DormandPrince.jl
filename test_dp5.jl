include("cleaner_dp5.jl")

# FCN must follow format FCN(N, X, Y, F, RPAR, IPAR)

function fcn(n, x, y, f)
    f[1] = ℯ^x
end

y = [1.0]
work = zeros(29)
iwork = zeros(21)

dopri5(
    1,
    fcn,
    0.0,
    y,
    2.0,
    1e-10,
    1e-10,
    0,
    work,# work, should be 8*N + 5*NRDENS+21 = 8*2 + 5*0 + 21
    iwork, #length should be NRDENS + 21
)