using Bloqade
using BloqadeExpr: Hamiltonian
include("dp5_solver_interface.jl")
# include("dp5_types.jl")

nsites = 1;
atoms = generate_sites(ChainLattice(), nsites, scale = 5.74)
h = rydberg_h(atoms; Ω = 20 * 2π, Δ = 0)
reg = zero_state(1)
# stop short of creating the SchrodingerProblem, we want the SchrodignerEquation
# which we can than trivially wrap with the fcn(n,x,y,f) that DP5 accepts

# Generate SchrodingerEquation
state = statevec(reg) 
space = YaoSubspaceArrayReg.space(reg)
T = real(eltype(state))
T = isreal(h) ? T : Complex{T}
eq = SchrodingerEquation(h, Hamiltonian(T, h, space))
# invoke eq via eq(dstate, state, p, t::Number) 

function fcn(x, y, f)
    eq(f, y, nothing, x)
end

y = statevec(reg)

solver = DP5Solver(
    fcn,
    0.0,
    y
)

dopri5(solver, 4.1)

