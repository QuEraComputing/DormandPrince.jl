include("dp5.jl")


# differential equation to solve for
function fcn(n, x, y, f)
    f[1] = y[1] - x^2 + 1
end

# initial y
y = [0.5]

# intermediate times to save
times = [1.0, 2.0, 3.0, 4.0]
# place to save values
values = []
# initial x
x = 0.0

for t in times
    xend = t 

    dopri5(
        1,
        fcn,
        x,
        y,
        xend,
        1e-10,
        1e-10,
        0,
        DP5Options(),
        zeros(8),# work, should be 8*N + 5*NRDENS+21 = 8*2 + 5*0 + 21
    )

    push!(values, y[1])
    global x = xend # set next initial x to last x
    
end