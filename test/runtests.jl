
# test error paths of dopcor, bypass checks if necessary 

using Test

@testset "Interface" begin
    include("interface.jl")
end

@testset "Checks" begin
    include("checks.jl")
end


@testset "Stiff ODE" begin
    include("stiff.jl")
end

@testset "Error Paths" begin
    include("errors.jl")
end

@testset "Pretty Printing" begin
    include("show.jl")
end