module DormandPrince

using Base.Iterators:repeated, Repeated

# internal imports
include("types.jl")
include("hinit.jl")
include("checks.jl")
include("interface.jl")
include("dp5/mod.jl")


using DormandPrince. DP5: core_integrator

# export Interface
export DP5Solver, integrate


end # DormandPrince
