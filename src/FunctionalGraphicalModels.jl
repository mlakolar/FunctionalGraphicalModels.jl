module FunctionalGraphicalModels

using LinearAlgebra
import BSpline
import ProximalBase

export DenseFunctionalData, quadrature_weights, score_covariance, smooth_bspline,
       FPCAResult, fit_fpca, transform_fpca, reconstruct,
       functional_glasso, FGLassoResult, FGLassoDiagnostic,
       functional_neighborhood, neighborhood_lambda_max, select_neighborhood_scv,
       NeighborhoodResult, NeighborhoodDiagnostic, NeighborhoodSCVResult

include("data.jl")
include("fpca.jl")
include("glasso.jl")
include("neighborhood.jl")

end
