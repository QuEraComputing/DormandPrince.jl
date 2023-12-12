using BenchmarkTools

struct Test{StateType}
    x::StateType
end


function test1!(y, a, b, c)
    y = a + 2.0 * ( b + 0.5 * c)
end

function test2!(y, a, b, c)
    y .= a + 2.0 * ( b + 0.5 * c)
end

function test3!(y, a, b, c)
    y .= a .+ 2.0 * ( b .+ 0.5 * c)
end


function test4!(y, a, b, c)
    y .= a .+ 2.0 .* ( b .+ 0.5 .* c)
end


y = zeros(1000)

test = Test(y)


a = rand(1000)
b = rand(1000)
c = rand(1000)

@benchmark test1!(test.x, a, b, c)
@benchmark test2!(test.x, a, b, c)
@benchmark test3!(test.x, a, b, c)
@benchmark test4!(test.x, a, b, c)


