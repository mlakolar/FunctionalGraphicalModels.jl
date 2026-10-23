using FunctionalGraphicalModels
using LinearAlgebra
using Test
using Aqua

@testset "FunctionalGraphicalModels" begin
    include("data.jl")
    include("fpca.jl")
    include("glasso.jl")
    include("neighborhood.jl")
    include("doubly.jl")
    include("separable.jl")
    include("simulation.jl")
    Aqua.test_all(FunctionalGraphicalModels)
end
