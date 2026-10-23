using FunctionalGraphicalModels
using LinearAlgebra
using Test
using Aqua

@testset "FunctionalGraphicalModels" begin
    include("data.jl")
    include("fpca.jl")
    Aqua.test_all(FunctionalGraphicalModels)
end
