# KernelSmoothing renamed these functions on master. Keep the same calls and
# numerical objectives with both current master and the frozen campaign lock.
import KernelSmoothing
function _dfgm_local_polynomial(args...;kwargs...)
    name=isdefined(KernelSmoothing,:local_polynomial) ? :local_polynomial : :localpolynomial
    getfield(KernelSmoothing,name)(args...;kwargs...)
end
function _dfgm_smoother_matrix(args...;kwargs...)
    name=isdefined(KernelSmoothing,:smoother_matrix) ? :smoother_matrix : :smoothermatrix
    getfield(KernelSmoothing,name)(args...;kwargs...)
end
