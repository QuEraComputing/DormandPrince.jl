module DormandPrince

using Base.Iterators:repeated, Repeated

# internal imports
include("types.jl")
include("hinit.jl")
include("checks.jl")
include("interface.jl")
include("dp5/mod.jl")
include("dp8/mod.jl")

# export Interface
export AbstractDPSolver, 
    DP5Solver, 
    DP8Solver, 
    integrate!, 
    SolverIterator,
    get_current_state


end # DormandPrince
