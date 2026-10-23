module FunctionalGraphicalModels

using LinearAlgebra
import BSpline

export DenseFunctionalData, quadrature_weights, score_covariance, smooth_bspline,
       FPCAResult, fit_fpca, transform_fpca, reconstruct

include("data.jl")
include("fpca.jl")

end
