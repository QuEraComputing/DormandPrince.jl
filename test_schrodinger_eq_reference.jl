using Bloqade

atoms = generate_sites(ChainLattice(), 1, scale = 5.74)
h = rydberg_h(atoms; Ω = 20 * 2π, Δ = 0)
reg = zero_state(1)

prob = SchrodingerProblem(reg, 4.1, h)
emulate!(prob)

statevec(prob.reg)