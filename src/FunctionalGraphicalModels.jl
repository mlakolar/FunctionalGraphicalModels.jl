module FunctionalGraphicalModels

using LinearAlgebra
import CovarianceSelection
import Random
import BSplineKit
import ProximalBase
import FunctionalDataAnalysis as FDA
import FunctionalDataAnalysis: quadrature_weights
import KernelSmoothing

export quadrature_weights, score_covariance, smooth_bspline,
       functional_glasso, FGLassoResult, FGLassoDiagnostic,
       functional_neighborhood, neighborhood_lambda_max, select_neighborhood_scv,
       NeighborhoodResult, NeighborhoodDiagnostic, NeighborhoodSCVResult,
       DoublyFunctionalFit, fit_doubly,
       JointGraphicalLassoFit, joint_graphical_lasso, SeparableFit, fit_separable,
       simulate_qiao2019, graph_metrics, FGMComponentSelection, select_fgm_components

include("data.jl")
include("tuning.jl")
include("glasso.jl")
include("neighborhood.jl")
include("doubly.jl")
include("separable.jl")
include("simulation.jl")

end
