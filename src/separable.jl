"""
    JointGraphicalLassoFit

The sparse ADMM variable `precision` and the positive definite likelihood
variable `likelihood_precision` for the joint graphical lasso. They agree to
the reported primal residual when `converged`. Sparsity is read from
`precision`; a prematurely stopped sparse variable need not be positive
definite. The objective uses the sparse variable and is `Inf` outside the cone.
"""
struct JointGraphicalLassoFit
    precision::Array{Float64,3}
    likelihood_precision::Array{Float64,3}
    converged::Bool
    iterations::Int
    primal_residual::Float64
    dual_residual::Float64
    objective::Float64
    positive_definite::BitVector
end

function _joint_penalty(A, lambda, alpha)
    p, _, m = size(A)
    value = 0.0
    for j in 1:p, k in 1:p
        j == k && continue
        value += alpha * sum(abs, @view A[j, k, :])
        value += (1 - alpha) * norm(@view A[j, k, :])
    end
    return lambda * value
end

function _joint_objective(covariances, A, lambda, alpha)
    value = _joint_penalty(A, lambda, alpha)
    for l in axes(A, 3)
        C = cholesky(Symmetric(@view A[:, :, l]); check=false)
        isposdef(C) || return Inf
        value += dot(@view(covariances[:, :, l]), @view(A[:, :, l])) - logdet(C)
    end
    return value
end

"""
    joint_graphical_lasso(S; lambda, alpha=0.5, rho=1.0, maxiter=2000,
                         abstol=1e-6, reltol=1e-5)

For symmetric positive semidefinite `p × p × M` input matrices, minimize

    sum_l (tr(S_l * Θ_l) - logdet(Θ_l)) +
    lambda * (alpha * sum_l sum_{j≠k} abs(Θ_l[j,k]) +
              (1-alpha) * sum_{j≠k} norm(Θ_[j,k]))

over positive definite matrices. Diagonals are unpenalized. Both triangular
halves are counted in the penalty, as in equations (10)–(11) of Zapata, Oh,
and Petersen (2022). The likelihood is not multiplied by sample size.

At zero penalty, every input slice must be positive definite; a singular
sample covariance has no finite maximum likelihood precision. For positive
penalties, a negative eigenvalue within the numerical PSD tolerance is
accepted only when a positive definite dual-feasible perturbation can be
certified. Inputs and solver parameters are converted to `Float64`.

ADMM uses `ProximalBase.LogDetLoss` and the exact sparse-group proximal map:
entrywise soft thresholding followed by group shrinkage across components.
The returned residuals and `converged` flag must be inspected after stopping.
"""
function joint_graphical_lasso(S::AbstractArray{<:Real,3}; lambda::Real,
                              alpha::Real=0.5, rho::Real=1.0,
                              maxiter::Integer=2000, abstol::Real=1e-6,
                              reltol::Real=1e-5)
    p, q, m = size(S)
    p == q && p > 0 && m > 0 || throw(DimensionMismatch("S must have size p × p × M with p,M positive"))
    all(isfinite, S) || throw(ArgumentError("S must be finite"))
    isfinite(lambda) && lambda >= 0 || throw(ArgumentError("lambda must be finite and nonnegative"))
    original_lambda = lambda
    lambda, alpha, rho = Float64(lambda), Float64(alpha), Float64(rho)
    abstol, reltol = Float64(abstol), Float64(reltol)
    isfinite(lambda) && lambda >= 0 || throw(ArgumentError("lambda must be finite and nonnegative"))
    original_lambda > 0 && lambda == 0 && throw(ArgumentError("positive lambda must be representable in Float64"))
    isfinite(alpha) && 0 <= alpha <= 1 || throw(ArgumentError("alpha must lie in [0,1]"))
    isfinite(rho) && rho > 0 || throw(ArgumentError("rho must be finite and positive"))
    isfinite(inv(rho)) && inv(rho) > 0 || throw(ArgumentError("inverse rho must be finite and positive in Float64"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    isfinite(abstol) && abstol > 0 || throw(ArgumentError("abstol must be finite and positive"))
    isfinite(reltol) && reltol > 0 || throw(ArgumentError("reltol must be finite and positive"))
    covariances = Array{Float64,3}(S)
    all(isfinite, covariances) || throw(ArgumentError("S entries must be finite in Float64"))
    mineigenvalues = zeros(m)
    diagonal_minima = zeros(m)
    spectral_scales = zeros(m)
    for l in 1:m
        sl = @view covariances[:, :, l]
        isapprox(sl, transpose(sl); atol=1e-12, rtol=1e-10) || throw(ArgumentError("S slices must be symmetric"))
        sl .= sl ./ 2 .+ transpose(sl) ./ 2
        mineigenvalues[l] = minimum(eigvals(Symmetric(sl)))
        spectral_scales[l] = max(1, opnorm(sl))
        mineigenvalues[l] >= -1e-10 * spectral_scales[l] ||
            throw(ArgumentError("S slices must be positive semidefinite"))
        diagonal_minima[l] = minimum(diag(sl))
        diagonal_minima[l] > 0 || throw(ArgumentError("S diagonals must be positive"))
        lambda == 0 && !isposdef(Symmetric(sl)) &&
            throw(ArgumentError("lambda=0 requires positive definite S slices; singular likelihoods have no finite minimizer"))
    end
    if lambda > 0 && any(<(0), mineigenvalues)
        # A sufficient boundedness certificate: mix each covariance with its
        # positive diagonal by the same fraction. The dual perturbation has
        # zero diagonal. Its off-diagonal vectors must lie in the Minkowski
        # sum of the l-infinity and l2 penalty balls. Soft thresholding is the
        # exact distance to the first ball, so the check below is explicit.
        fraction = min(1.0, maximum((max(0, -mineigenvalues[l]) +
            100eps(Float64)*spectral_scales[l]) /
            (diagonal_minima[l]+max(0,-mineigenvalues[l])) for l in 1:m))
        feasible = true
        for l in 1:m
            candidate = (1-fraction) * covariances[:,:,l] +
                        fraction * Diagonal(diag(covariances[:,:,l]))
            feasible &= isposdef(Symmetric(candidate))
        end
        for j in 1:p, k in (j+1):p
            residual = sqrt(sum(max(fraction*abs(covariances[j,k,l])-
                lambda*alpha, 0)^2 for l in 1:m))
            feasible &= residual <= lambda*(1-alpha)
        end
        feasible || throw(ArgumentError("numerically indefinite S is not certified dual feasible at this lambda; use a PSD covariance or a larger penalty"))
    end
    X = zeros(p, p, m)
    Z = zeros(p, p, m)
    for l in 1:m, j in 1:p
        X[j, j, l] = Z[j, j, l] = 1 / covariances[j, j, l]
    end
    U = zeros(p, p, m)
    previous = similar(Z)
    work = zeros(p, p)
    losses = [ProximalBase.LogDetLoss(Matrix(@view covariances[:, :, l])) for l in 1:m]
    converged = false
    primal = dual = Inf
    iterations = 0
    for iteration in 1:maxiter
        iterations = iteration
        copyto!(previous, Z)
        for l in 1:m
            work .= @view(Z[:, :, l]) .- @view(U[:, :, l])
            ProximalBase.prox!(@view(X[:, :, l]), losses[l], work, inv(rho))
        end
        Z .= X .+ U
        for j in 1:p, k in (j + 1):p
            groupnorm2 = 0.0
            for l in 1:m
                # The two entries coincide analytically; averaging removes
                # roundoff asymmetry before applying the symmetric map.
                z = (Z[j, k, l] + Z[k, j, l]) / 2
                z = sign(z) * max(abs(z) - lambda * alpha / rho, 0)
                Z[j, k, l] = Z[k, j, l] = z
                groupnorm2 += z^2
            end
            groupnorm = sqrt(groupnorm2)
            scale = groupnorm == 0 ? 0.0 : max(1 - lambda * (1-alpha) / (rho * groupnorm), 0)
            for l in 1:m
                Z[j, k, l] *= scale
                Z[k, j, l] = Z[j, k, l]
            end
        end
        U .+= X .- Z
        primal = norm(X - Z)
        dual = rho * norm(Z - previous)
        eps_primal = sqrt(length(X)) * abstol + reltol * max(norm(X), norm(Z))
        eps_dual = sqrt(length(X)) * abstol + reltol * rho * norm(U)
        if primal <= eps_primal && dual <= eps_dual
            converged = true
            break
        end
    end
    positive = BitVector([isposdef(Symmetric(@view Z[:, :, l])) for l in 1:m])
    return JointGraphicalLassoFit(Z, X, converged, iterations, primal, dual,
        _joint_objective(covariances, Z, lambda, alpha), positive)
end

"""A shared-basis functional graph fit; `adjacency` is the union of component graphs."""
struct SeparableFit{F}
    precisions::Array{Float64,3}
    inverse_correlations::Array{Float64,3}
    correlations::Array{Float64,3}
    scales::Matrix{Float64}
    component_adjacency::BitArray{3}
    adjacency::BitMatrix
    fpca::F
    diagnostics::JointGraphicalLassoFit
end

"""
    fit_separable(data::DenseFunctionalData; lambda, alpha=0.5,
                  ncomponents=2, threshold=1e-6, kwargs...)

Estimate the partially separable functional graphical model of Zapata, Oh,
and Petersen (2022). FPCA uses the pooled marginal covariance eigenbasis.
Each component is standardized to correlations before fitting the joint
graphical lasso. `precisions` rescales inverse correlations back to the score
units. `threshold` is used only to extract graph support.

The union of component supports describes the truncated functional graph under
Gaussianity and partial separability. If partial separability fails, these
are graphs of the individual projected score vectors; no equivalence to the
full-process conditional independence graph is asserted.
"""
function fit_separable(data::DenseFunctionalData; lambda::Real,
                       alpha::Real=0.5, ncomponents::Integer=2,
                       threshold::Real=1e-6, kwargs...)
    isfinite(threshold) && threshold >= 0 || throw(ArgumentError("threshold must be finite and nonnegative"))
    fpca = fit_fpca(data; ncomponents=ncomponents, basis=:pooled)
    scores = fpca.scores
    n, p, m = size(scores)
    correlations = zeros(p, p, m)
    scales = zeros(p, m)
    for l in 1:m
        A = Matrix(@view scores[:, :, l])
        A .-= sum(A; dims=1) ./ n
        S = transpose(A) * A / n
        scales[:, l] .= sqrt.(diag(S))
        all(>(0), @view(scales[:, l])) || throw(ArgumentError("each retained score component must have positive variance at every node"))
        correlations[:, :, l] .= S ./ (scales[:, l] * transpose(scales[:, l]))
    end
    solver = joint_graphical_lasso(correlations; lambda=lambda, alpha=alpha, kwargs...)
    precision = similar(solver.precision)
    supports = falses(p, p, m)
    adjacency = falses(p, p)
    for l in 1:m
        precision[:, :, l] .= solver.precision[:, :, l] ./ (scales[:, l] * transpose(scales[:, l]))
        for j in 1:p, k in (j + 1):p
            supports[j, k, l] = supports[k, j, l] = abs(solver.precision[j, k, l]) > threshold
            adjacency[j, k] = adjacency[k, j] = adjacency[j, k] || supports[j, k, l]
        end
    end
    return SeparableFit(precision, solver.precision, correlations, scales,
                        supports, adjacency, fpca, solver)
end
