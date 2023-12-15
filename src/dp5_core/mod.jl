module DP5Core

# external imports
using Base.Iterators:repeated, Repeated

include("types.jl")
include("helpers.jl")
include("checks.jl")
include("solver.jl")

export DP5Solver, dopri5


end