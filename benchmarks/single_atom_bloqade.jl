using Bloqade
using BenchmarkTools

function generate_clean_problem()
    atoms = generate_sites(ChainLattice(), 1, scale = 5.74)
    h = rydberg_h(atoms; Ω = 20 * 2π, Δ = 0)
    reg = zero_state(1)
    
    prob = SchrodingerProblem(reg, 1.6, h; algo=DP5())

    return prob
end

emulate!(generate_clean_problem())

# b = @benchmark emulate!(clean_problem) samples=100 evals=1 seconds=172800 setup=(clean_problem=SchrodingerProblem(deepcopy(reg), 1.6, h))
b = @benchmark emulate!(clean_problem) samples=10000 evals=5 seconds=172800 setup=(clean_problem=generate_clean_problem())

# statevec(prob.reg)