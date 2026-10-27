module FunctionalGraphicalModels
using LinearAlgebra
import SparseArrays,JuMP,Clarabel,HiGHS
module FDA
abstract type MultivariateFunctionalData end
end
include(joinpath(ENV["FGM_SOURCE"],"src/direct.jl"))
include(joinpath(ENV["FGM_SOURCE"],"src/direct_path.jl"))
include(joinpath(ENV["FGM_SOURCE"],"src/spectral_direct.jl"))
include(joinpath(ENV["FGM_SOURCE"],"src/multiscale_path.jl"))
end
using .FunctionalGraphicalModels,LinearAlgebra,Test
include(joinpath(ENV["FGM_SOURCE"],"test/multiscale_path.jl"))
