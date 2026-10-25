"""
A pointwise doubly functional graph. Unlike a static functional graph,
`adjacency[:,:,t]` concerns scalar values at the particular time `grid[t]`.
CLIME estimates are symmetric but need not be positive definite; this is
reported explicitly in `positive_definite`.
"""
struct DoublyFunctionalFit{F}
    grid::Vector{Float64}
    covariances::Array{Float64,3}
    precisions::Array{Float64,3}
    adjacency::BitArray{3}
    reconstruction::Array{Float64,3}
    fpca::F
    method::Symbol
    positive_definite::BitVector
    lambda::Vector{Float64}
    threshold::Float64
end

function _fit_doubly(fpca; lambda, threshold, method)
    isfinite(threshold) && threshold >= 0 || throw(ArgumentError("threshold must be finite and nonnegative"))
    scores = FDA.score_array(fpca)
    n, p, _ = size(scores)
    grid = first(fpca.bases).grid
    all(b -> b.grid == grid, fpca.bases) ||
        throw(ArgumentError("doubly functional FPCA requires the same grid for every node"))
    q = length(grid)
    penalties = lambda isa Real ? fill(Float64(lambda),q) : collect(Float64, lambda)
    length(penalties) == q || throw(DimensionMismatch("lambda must be a scalar or have one value per grid point"))
    all(isfinite, penalties) && all(>=(0), penalties) || throw(ArgumentError("lambda must be finite and nonnegative"))
    reconstructed = zeros(n, p, q)
    covariance = zeros(p, p, q)
    precision = similar(covariance)
    adjacency = falses(p, p, q)
    positive = falses(q)
    for j in 1:p
        reconstructed[:,j,:] .= scores[:,j,:] * transpose(fpca.bases[j].values)
    end
    for t in 1:q
        A = Matrix(@view reconstructed[:,:,t])
        # Equation (2), with the estimated marginal means already removed.
        # Conditional scores under irregular schedules need not have exactly
        # zero empirical mean; the paper uses this uncentered second moment.
        covariance[:,:,t] .= transpose(A) * A / n
        precision[:,:,t] .= Matrix(CovarianceSelection.clime(Matrix(@view covariance[:,:,t]), penalties[t]))
        positive[t] = isposdef(Symmetric(@view precision[:,:,t]))
        for j in 1:p, k in (j+1):p
            adjacency[j,k,t] = adjacency[k,j,t] = abs(precision[j,k,t]) > threshold
        end
    end
    DoublyFunctionalFit(Float64.(grid), covariance, precision, adjacency,
                        reconstructed, fpca, method, positive, copy(penalties), Float64(threshold))
end


"""
    fit_doubly(data::FDA.MultivariateFunctionalData; lambda, ncomponents=2,
               threshold=1e-6, grid=nothing, bandwidth=nothing,
               noise_variance=nothing, kernel=:gaussian)

Qiao et al. (2020) doubly functional graph with marginal FDA FPCA and
pointwise columnwise CLIME. `lambda` is required and may be a nonnegative
scalar or one value per grid point. `threshold=1e-6` and `ncomponents=2`
retain FGM's explicit graph threshold and rank defaults. Reconstructed
pointwise covariances use denominator `n`, as in the paper, and exclude means.
CLIME uses the paper's smaller-absolute-entry symmetrization.

Complete grid nodes on one grid use dense quadrature FPCA. Its eigenvalues
use FDA's `n-1` denominator; projected scores and graph covariances retain
their former scaling. FDA resolves signs by mean alignment, followed by the
largest-absolute-entry rule when alignment is negligible.

Irregular nodes require `grid` and a positive `bandwidth`. FDA's marginal
PACE uses the supplied bandwidth for mean, covariance and diagonal smoothing,
`weighting=:subj`, no binning, and all positive covariance modes in the
conditional-score covariance. The relative eigenvalue cutoff `100eps(Float64)`
retains FGM's former spectrum cutoff. This preserves FGM's smoothing settings
rather than FDA's fdapace-derived bandwidth and binning defaults.

`kernel=:gaussian` retains FGM's default and maps to
`KernelSmoothing.GaussianKernel()`; `:epanechnikov` maps to
`KernelSmoothing.EpanechnikovKernel()`. The latter satisfies the paper's
compact-support condition. Gaussian tails do not satisfy that condition, so
the stated theorems do not transfer unchanged. A singular local design raises
an error. `noise_variance` is a nonnegative scalar or one value per node.
When omitted, FDA integrates excess raw diagonal variance over the middle
half of the grid, following Yao, Müller and Wang (2005). FDA's settings store
all resolved PACE choices. Singular observation covariance raises an error
naming the curve; no ridge or pseudoinverse is used.

With bounded observations per curve, the sparse population target is the
conditional-score reconstructed process, as distinguished in Section 3 of
Qiao et al. (2020). It need not equal the original process's precision.
"""
function fit_doubly(data::FDA.MultivariateFunctionalData; lambda,
                    ncomponents::Integer=2, threshold::Real=1e-6,
                    grid=nothing, bandwidth=nothing,
                    noise_variance=nothing, kernel::Symbol=:gaussian)
    complete = all(c -> c isa FDA.GridFunctionalData && c.observed === nothing,
        data.components)
    if complete
        grid === nothing && bandwidth === nothing && noise_variance === nothing &&
            kernel === :gaussian || throw(ArgumentError(
                "grid, bandwidth, noise_variance and nondefault kernel require irregular nodes for fit_doubly"))
        _dense_parts(data)
        fpca = FDA.fit(FDA.FPCA, data; basis=:marginal,
            ncomponents=ncomponents, max_k=ncomponents, method=:dense)
        return _fit_doubly(fpca; lambda, threshold, method=:dense_reconstruction)
    end
    all(c -> c isa FDA.IrregularFunctionalData, data.components) ||
        throw(ArgumentError("fit_doubly requires complete grid nodes or irregular nodes"))
    grid isa AbstractVector{<:Real} || throw(ArgumentError(
        "grid must be supplied as a real vector for irregular fit_doubly"))
    bandwidth isa Real && isfinite(bandwidth) && bandwidth > 0 ||
        throw(ArgumentError("bandwidth must be finite and positive for irregular fit_doubly"))
    kernel in (:gaussian, :epanechnikov) ||
        throw(ArgumentError("kernel must be :gaussian or :epanechnikov"))
    g = collect(Float64, grid)
    length(g) >= 3 || throw(ArgumentError("grid must contain at least three points"))
    FDA.quadrature_weights(g)
    for j in 1:FDA.ncomponents(data)
        all(>=(2), FDA.nobs_per_curve(FDA.component(data, j))) ||
            throw(ArgumentError("every curve at node $j needs at least two observations"))
    end
    noise_variance === nothing || noise_variance isa Real ||
        noise_variance isa AbstractVector || throw(ArgumentError(
            "noise_variance must be a nonnegative scalar or a vector with one value per node"))
    p = FDA.ncomponents(data)
    noise = if noise_variance === nothing
        nothing
    elseif noise_variance isa Real
        converted = Float64(noise_variance)
        isfinite(converted) && converted >= 0 || throw(ArgumentError(
            "noise_variance must be finite and nonnegative"))
        converted
    else
        length(noise_variance) == p || throw(DimensionMismatch(
            "noise_variance must have one value per node"))
        all(x -> x isa Real && isfinite(x) && x >= 0, noise_variance) ||
            throw(ArgumentError("noise_variance must contain finite nonnegative numbers"))
        converted = Float64.(noise_variance)
        all(isfinite, converted) || throw(ArgumentError(
            "noise_variance must remain finite after conversion"))
        converted
    end
    bandwidth = Float64(bandwidth)
    isfinite(bandwidth) && bandwidth > 0 || throw(ArgumentError(
        "bandwidth must remain finite and positive after conversion"))
    smoothing_kernel = kernel === :gaussian ? KernelSmoothing.GaussianKernel() :
        KernelSmoothing.EpanechnikovKernel()
    fpca = FDA.fit(FDA.FPCA, data; basis=:marginal, method=:pace,
        grid=g, ncomponents=ncomponents, max_k=length(g),
        mean_bandwidth=bandwidth, covariance_bandwidth=bandwidth,
        diagonal_bandwidth=bandwidth, kernel=smoothing_kernel,
        weighting=:subj, binning=nothing, noise_variance=noise,
        eigen_tol=100eps(Float64))
    _fit_doubly(fpca; lambda, threshold, method=:sparse_conditional_scores)
end


function Base.show(io::IO, fit::DoublyFunctionalFit)
    print(io, "DoublyFunctionalFit(method=", fit.method, ", lambda=", fit.lambda,
        ", threshold=", fit.threshold, ", fpca=", fit.fpca, ")")
end
