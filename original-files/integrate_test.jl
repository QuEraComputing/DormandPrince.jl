include("integrate.jl")

using Bloqade
using BloqadeExpr: Hamiltonian

nsites = 1;
atoms = generate_sites(ChainLattice(), nsites, scale = 5.74)
h = rydberg_h(atoms; Ω = 11 * 2π, Δ = 0)
reg = zero_state(1)

state = statevec(reg) 
space = YaoSubspaceArrayReg.space(reg)
T = real(eltype(state))
T = isreal(h) ? T : Complex{T}
eq = SchrodingerEquation(h, Hamiltonian(T, h, space))

function fcn(x, y, f)
    eq(f, y, nothing, x)
end

y = statevec(reg)

solver = DP5Solver(
    fcn,
    0.0,
    y
)

iterator = integrate(solver, [1.8, 1.9, 2.0])


for (t, y) in iterator
    println(t)
    println(y)
end


#=
function simple_callback(time, state)
    println("State ", state, " at time ", time)
    return copy(state)
end

integrate(simple_callback, solver, [1.8, 1.9, 2.0])
=#