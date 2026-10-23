"""
    SparseFunctionalData(times, values)

Irregular functional observations. `times` and `values` are `n × p` matrices
whose entries are vectors of observation times and observations. Every curve
must have at least two distinct times. The constructor copies and sorts each
curve; observations are not interpolated onto a common grid.
"""
struct SparseFunctionalData
    times::Matrix{Vector{Float64}}
    values::Matrix{Vector{Float64}}
    function SparseFunctionalData(times::AbstractMatrix, values::AbstractMatrix)
        size(times) == size(values) || throw(DimensionMismatch("times and values must have the same n × p size"))
        n, p = size(times)
        n >= 2 && p >= 1 || throw(ArgumentError("at least two subjects and one node are required"))
        t = Matrix{Vector{Float64}}(undef, n, p)
        y = similar(t)
        for j in 1:p, i in 1:n
            ti, yi = Float64.(times[i, j]), Float64.(values[i, j])
            length(ti) == length(yi) || throw(DimensionMismatch("observation times and values differ for curve ($i,$j)"))
            length(ti) >= 2 || throw(ArgumentError("every curve needs at least two observations"))
            all(isfinite, ti) && all(isfinite, yi) || throw(ArgumentError("times and observations must be finite"))
            order = sortperm(ti)
            t[i, j], y[i, j] = ti[order], yi[order]
            all(>(0), diff(t[i, j])) || throw(ArgumentError("observation times must be distinct within each curve"))
        end
        new(t, y)
    end
end

"""
Sparse FPCA with conditional-expectation scores. `covariances[:,:,j]` is the
positive semidefinite projection of the off-diagonal covariance smoother.
`noise_variance[j]` is specified or estimated from the excess diagonal.
"""
struct SparseFPCAFit
    scores::Array{Float64,3}
    eigenfunctions::Array{Float64,3}
    eigenvalues::Matrix{Float64}
    mean::Matrix{Float64}
    weights::Vector{Float64}
    grid::Vector{Float64}
    covariances::Array{Float64,3}
    noise_variance::Vector{Float64}
    kernel::Symbol
    bandwidth::Float64
end

function _dfgm_weights(grid)
    all(isfinite, grid) && length(grid) >= 3 && all(>(0), diff(grid)) ||
        throw(ArgumentError("grid must contain at least three finite increasing points"))
    w = zeros(length(grid))
    w[1], w[end] = (grid[2]-grid[1])/2, (grid[end]-grid[end-1])/2
    for t in 2:(length(grid)-1)
        w[t] = (grid[t+1]-grid[t-1])/2
    end
    w
end

# Kernel normalizing constants cancel in each weighted least-squares fit.
function _dfgm_kernel(x, kernel)
    kernel == :gaussian && return exp(-x*x/2)
    kernel == :epanechnikov && return (abs(x) < 1 ? 1-x*x : 0.0)
    throw(ArgumentError("kernel must be :gaussian or :epanechnikov"))
end

function _dfgm_local_linear(t, y, weights, u, bandwidth, kernel=:gaussian)
    s0 = s1 = s2 = r0 = r1 = 0.0
    for a in eachindex(t)
        x = (t[a] - u) / bandwidth
        k = _dfgm_kernel(x, kernel) * weights[a]
        s0 += k
        s1 += k*x
        s2 += k*x*x
        r0 += k*y[a]
        r1 += k*x*y[a]
    end
    determinant = s0*s2 - s1*s1
    determinant > 100eps(Float64)*max(s0*s2, floatmin(Float64)) ||
        throw(ArgumentError("local-linear smoother is singular; increase bandwidth or observation coverage"))
    (s2*r0 - s1*r1) / determinant
end

function _dfgm_local_surface(t, v, y, weights, u0, v0, bandwidth, kernel=:gaussian)
    G = zeros(3, 3)
    r = zeros(3)
    for a in eachindex(y)
        x1, x2 = (t[a] - u0)/bandwidth, (v[a] - v0)/bandwidth
        k = weights[a] * (kernel == :gaussian ? exp(-(x1*x1+x2*x2)/2) :
            _dfgm_kernel(x1,kernel)*_dfgm_kernel(x2,kernel))
        G[1, 1] += k
        G[1, 2] += k*x1
        G[1, 3] += k*x2
        G[2, 2] += k*x1*x1
        G[2, 3] += k*x1*x2
        G[3, 3] += k*x2*x2
        r[1] += k*y[a]
        r[2] += k*x1*y[a]
        r[3] += k*x2*y[a]
    end
    G[2, 1], G[3, 1], G[3, 2] = G[1, 2], G[1, 3], G[2, 3]
    C = cholesky(Symmetric(G); check=false)
    isposdef(C) || throw(ArgumentError("local-linear covariance smoother is singular; increase bandwidth or coverage"))
    (C \ r)[1]
end

function _dfgm_interpolate(grid, values::AbstractVector, x)
    first(grid) <= x <= last(grid) || throw(ArgumentError("interpolation time is outside the evaluation grid"))
    k = clamp(searchsortedlast(grid, x), 1, length(grid)-1)
    a = (x-grid[k]) / (grid[k+1]-grid[k])
    (1-a)*values[k] + a*values[k+1]
end

"""
    fit_sparse_fpca(data::SparseFunctionalData; grid, ncomponents=2,
                    bandwidth, noise_variance=nothing, kernel=:gaussian)

Estimate marginal FPCA for the sparse-observation DFGM procedure. For each
node, local-linear mean smoothing is followed by a local-linear surface fit
to centered *off-diagonal within-curve* products. Subjects receive equal
total weight. The covariance estimate is symmetrized and projected onto the
positive semidefinite cone using trapezoidal quadrature. Conditional scores
are `Λ Φ' (Φ Λ Φ' + σ² I)⁻¹ (Y-μ)`.

`bandwidth` is an explicit positive bandwidth in the time units of `grid`.
`kernel=:gaussian` retains observations at all distances. The paper's
compact-support kernel condition does not hold for this default; its stated
theorems do not transfer unchanged. `kernel=:epanechnikov` uses a Lipschitz
kernel supported on `[-1,1]`, matching that kernel condition. Compact support
can leave a singular local design, which raises an error rather than silently
changing the bandwidth or fitting a different model.

`noise_variance` may be a nonnegative scalar or one value per node. If omitted,
the average excess of the smoothed raw diagonal over the covariance diagonal
on the middle half of the grid is truncated at zero. This explicit noise
rule is a package default, not an implementation of every fdapace tuning
choice. No bandwidth cross-validation or rank AIC selection is performed.

The evaluation grid must cover every observation time. Linear interpolation
evaluates estimated means and eigenfunctions between grid points; this method
does not extrapolate observations outside the fitted domain.
No diagonal ridge is added. A Moore–Penrose inverse handles a singular
observation covariance when the specified measurement noise is zero.
"""
function fit_sparse_fpca(data::SparseFunctionalData; grid::AbstractVector,
                         ncomponents::Integer=2, bandwidth::Real,
                         noise_variance=nothing, kernel::Symbol=:gaussian)
    bandwidth = Float64(bandwidth)
    isfinite(bandwidth) && bandwidth > 0 || throw(ArgumentError("bandwidth must be finite and positive"))
    kernel in (:gaussian,:epanechnikov) || throw(ArgumentError("kernel must be :gaussian or :epanechnikov"))
    g = collect(Float64, grid)
    w = _dfgm_weights(g)
    n, p = size(data.values)
    all(t -> first(g) <= first(t) && last(t) <= last(g), data.times) ||
        throw(ArgumentError("grid must cover every observation time"))
    q = length(g)
    1 <= ncomponents <= q || throw(ArgumentError("ncomponents must lie between one and the grid length"))
    noise = if noise_variance === nothing
        fill(NaN, p)
    elseif noise_variance isa Real
        fill(Float64(noise_variance), p)
    else
        length(noise_variance) == p || throw(DimensionMismatch("noise_variance must have one value per node"))
        collect(Float64, noise_variance)
    end
    if noise_variance !== nothing
        all(isfinite, noise) && all(>=(0), noise) || throw(ArgumentError("noise_variance must be finite and nonnegative"))
    end
    scores = zeros(n, p, ncomponents)
    eigenfunctions = zeros(q, ncomponents, p)
    eigenvalues = zeros(ncomponents, p)
    means = zeros(p, q)
    covariances = zeros(q, q, p)
    rootw = sqrt.(w)
    for j in 1:p
        times, values, pointweights = Float64[], Float64[], Float64[]
        for i in 1:n
            count = length(data.times[i, j])
            append!(times, data.times[i, j])
            append!(values, data.values[i, j])
            append!(pointweights, fill(1/(n*count), count))
        end
        for t in 1:q
            means[j, t] = _dfgm_local_linear(times, values, pointweights, g[t], bandwidth, kernel)
        end
        centered = [data.values[i,j] .- [_dfgm_interpolate(g, @view(means[j,:]), t) for t in data.times[i,j]] for i in 1:n]
        firsttime, secondtime, products, pairweights = Float64[], Float64[], Float64[], Float64[]
        squaredvalues = Float64[]
        for i in 1:n
            ti, yi = data.times[i,j], centered[i]
            count = length(ti)
            append!(squaredvalues, yi.^2)
            weight = 1/(n*count*(count-1))
            for a in 1:count, b in 1:count
                a == b && continue
                push!(firsttime, ti[a]); push!(secondtime, ti[b])
                push!(products, yi[a]*yi[b]); push!(pairweights, weight)
            end
        end
        K = zeros(q, q)
        for t in 1:q, s in t:q
            K[t,s] = K[s,t] = _dfgm_local_surface(firsttime, secondtime,
                products, pairweights, g[t], g[s], bandwidth, kernel)
        end
        decomposition = eigen(Symmetric((rootw * transpose(rootw)) .* K))
        ordering = sortperm(decomposition.values; rev=true)
        positive = filter(k -> decomposition.values[k] >
            100eps(Float64)*max(maximum(abs, decomposition.values), floatmin(Float64)), ordering)
        length(positive) >= ncomponents || throw(ArgumentError("smoothed covariance at node $j has fewer than ncomponents positive eigenvalues"))
        allvalues = decomposition.values[positive]
        allfunctions = decomposition.vectors[:, positive] ./ rootw
        # Eigenfunction signs are fixed only to make saved results deterministic.
        for l in axes(allfunctions, 2)
            pivot = argmax(abs.(@view allfunctions[:,l]))
            if allfunctions[pivot,l] < 0
                allfunctions[:,l] .*= -1
            end
        end
        covariances[:,:,j] .= (allfunctions .* transpose(allvalues)) * transpose(allfunctions)
        eigenfunctions[:,:,j] .= allfunctions[:,1:ncomponents]
        eigenvalues[:,j] .= allvalues[1:ncomponents]
        if noise_variance === nothing
            interior = max(1, ceil(Int,q/4)):min(q, floor(Int,3q/4))
            diagonal_excess = [_dfgm_local_linear(times, squaredvalues,
                pointweights, g[t], bandwidth, kernel)-covariances[t,t,j] for t in interior]
            noise[j] = max(sum(diagonal_excess)/length(diagonal_excess), 0)
        end
        for i in 1:n
            ti, yi = data.times[i,j], centered[i]
            Phi = [_dfgm_interpolate(g, @view(allfunctions[:,l]), t) for t in ti, l in eachindex(allvalues)]
            observation_covariance = (Phi .* transpose(allvalues)) * transpose(Phi)
            for a in axes(observation_covariance,1)
                observation_covariance[a,a] += noise[j]
            end
            factor = cholesky(Symmetric(observation_covariance); check=false)
            transformed = isposdef(factor) ? factor \ yi : pinv(Symmetric(observation_covariance)) * yi
            scores[i,j,:] .= allvalues[1:ncomponents] .* (transpose(Phi[:,1:ncomponents]) * transformed)
        end
    end
    SparseFPCAFit(scores, eigenfunctions, eigenvalues, means, w, g,
                   covariances, noise, kernel, bandwidth)
end

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
end

function _fit_doubly(fpca; lambda, threshold, method)
    isfinite(threshold) && threshold >= 0 || throw(ArgumentError("threshold must be finite and nonnegative"))
    n, p, m = size(fpca.scores)
    q = length(fpca.grid)
    penalties = lambda isa Real ? fill(Float64(lambda),q) : collect(Float64, lambda)
    length(penalties) == q || throw(DimensionMismatch("lambda must be a scalar or have one value per grid point"))
    all(isfinite, penalties) && all(>=(0), penalties) || throw(ArgumentError("lambda must be finite and nonnegative"))
    reconstructed = zeros(n, p, q)
    covariance = zeros(p, p, q)
    precision = similar(covariance)
    adjacency = falses(p, p, q)
    positive = falses(q)
    for j in 1:p
        reconstructed[:,j,:] .= fpca.scores[:,j,:] * transpose(fpca.eigenfunctions[:,:,j])
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
    DoublyFunctionalFit(Float64.(fpca.grid), covariance, precision, adjacency,
                        reconstructed, fpca, method, positive)
end

"""
    fit_doubly(data::DenseFunctionalData; lambda, ncomponents=2, threshold=1e-6)

The dense, noise-free reconstruction version of the Qiao et al. (2020) doubly
functional graphical model: marginal quadrature FPCA, reconstructed pointwise
covariance with denominator `n`, and exact columnwise CLIME from
`CovarianceSelection.clime`. CLIME applies the paper's smaller-absolute-entry
symmetrization. `lambda` may vary over the grid. Means are excluded from
reconstructions because the target is covariance.

This dense path does not estimate measurement noise or smooth raw covariance
diagonals. Use the irregular-data method for noisy or sparse observations.
"""
function fit_doubly(data::DenseFunctionalData; lambda, ncomponents::Integer=2,
                    threshold::Real=1e-6)
    fpca = fit_fpca(data; ncomponents=ncomponents, basis=:marginal)
    _fit_doubly(fpca; lambda=lambda, threshold=threshold, method=:dense_reconstruction)
end

"""
    fit_doubly(data::SparseFunctionalData; lambda, grid, bandwidth,
               ncomponents=2, noise_variance=nothing, kernel=:gaussian,
               threshold=1e-6)

Sparse/noisy DFGM using `fit_sparse_fpca` followed by pointwise CLIME. With a
bounded number of observations per curve, the population target is the
conditional-score reconstructed process, as distinguished in Section 3 of
Qiao et al. (2020). It need not equal the original process's precision.
"""
function fit_doubly(data::SparseFunctionalData; lambda, grid::AbstractVector,
                    bandwidth::Real, ncomponents::Integer=2,
                    noise_variance=nothing, kernel::Symbol=:gaussian,
                    threshold::Real=1e-6)
    fpca = fit_sparse_fpca(data; grid=grid, ncomponents=ncomponents,
                            bandwidth=bandwidth, noise_variance=noise_variance,
                            kernel=kernel)
    _fit_doubly(fpca; lambda=lambda, threshold=threshold, method=:sparse_conditional_scores)
end
