using FunctionalGraphicalModels
using LinearAlgebra
using Random
using Test
using Aqua
import FunctionalDataAnalysis as FDA

@testset "FunctionalGraphicalModels" begin
    include("data.jl")
    include("fpca.jl")
    include("tuning.jl")
    include("glasso.jl")
    include("neighborhood.jl")
    include("doubly.jl")
    include("separable.jl")
    include("simulation.jl")
    Aqua.test_all(FunctionalGraphicalModels)
end
