using Bloqade

atoms = generate_sites(ChainLattice(), 1, scale = 5.74)
h = rydberg_h(atoms; Ω = 11 * 2π, Δ = 0)
reg = zero_state(1)

prob = SchrodingerProblem(reg, 1.6, h)
emulate!(prob)

statevec(prob.reg)