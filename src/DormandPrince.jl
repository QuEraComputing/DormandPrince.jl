module DormandPrince

using Base.Iterators:repeated, Repeated

# internal imports
include("types.jl")
include("interface.jl")
include("dp5_impl/mod.jl")


using DormandPrince.DP5Impl: core_integrator

# export Interface
export DP5Solver, integrate


end # DormandPrince
