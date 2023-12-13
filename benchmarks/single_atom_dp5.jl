include("../solver.jl")
include("../types.jl")
using Bloqade
using BloqadeExpr: Hamiltonian
using BenchmarkTools

function generate_clean_solver()
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

    function fcn(x, y, f)
        eq(f, y, nothing, x)
    end
    
    y = statevec(reg)
    
    solver = DP5Solver(
        fcn,
        0.0,
        y
    )

    return solver
end

# trigger precompilation
dopri5(generate_clean_solver(), 1.6)

# b = @benchmark dopri5(clean_solver, 1.6) samples=100 evals=1 seconds=172800 setup=(clean_solver = DP5Solver(fcn, 0.0, deepcopy(y)))
@benchmark dopri5(clean_solver, 1.6) samples=10000 evals=5 seconds=172800 setup=(clean_solver = generate_clean_solver())