"""Convergence information for the functional graphical-lasso ADMM solver."""
struct FGLassoDiagnostic
    iterations::Int
    converged::Bool
    primal_residual::Float64
    dual_residual::Float64
    kkt_residual::Float64
    objective::Float64
end

"""
The positive-definite ADMM iterate is `precision`; `sparse_precision` is its
block-sparse consensus copy. Their difference is reported in `diagnostics`.
The graph comes from the sparse copy, with an explicit block-norm threshold.
"""
struct FGLassoResult
    precision::Matrix{Float64}
    sparse_precision::Matrix{Float64}
    covariance::Matrix{Float64}
    block_sizes::Vector{Int}
    adjacency::BitMatrix
    diagnostics::FGLassoDiagnostic
    fpca::Union{Nothing,FDA.MultivariateFPCA}
    settings::NamedTuple
end

function _fgm_blocks(d, sizes)
    bs = sizes isa Integer ? fill(Int(sizes), div(d, sizes > 0 ? sizes : 1)) : Int.(collect(sizes))
    !isempty(bs) && all(>(0), bs) && sum(bs) == d ||
        throw(ArgumentError("positive block sizes must sum to the covariance dimension"))
    stops = cumsum(bs)
    bs, [((j == 1 ? 0 : stops[j-1]) + 1):stops[j] for j in eachindex(bs)]
end

function _fgm_penalty(A, blocks)
    sum((norm(view(A, blocks[j], blocks[k])) for j in eachindex(blocks)
        for k in eachindex(blocks) if j != k); init=0.0)
end

function _fgm_graph(A, blocks, threshold)
    p = length(blocks)
    graph = falses(p, p)
    for j in 1:p, k in j+1:p
        graph[j,k] = graph[k,j] = norm(view(A, blocks[j], blocks[k])) > threshold
    end
    graph
end

function _fgm_kkt(S, X, Z, blocks, lambda)
    G = S - inv(Symmetric(X))
    worst = 0.0
    for j in eachindex(blocks), k in eachindex(blocks)
        g = view(G, blocks[j], blocks[k])
        z = view(Z, blocks[j], blocks[k])
        nz = norm(z)
        v = j == k ? norm(g) : nz > 0 ? norm(g + lambda .* z ./ nz) : max(0, norm(g) - lambda)
        worst = max(worst, v)
    end
    worst
end

"""
    functional_glasso(S, block_sizes, lambda; rho=1, maxiter=10000,
                      abstol=1e-7, reltol=1e-6, threshold=1e-6)
    functional_glasso(scores, lambda; kwargs...)
    functional_glasso(data, lambda; ncomponents, kwargs...)

Qiao, Guo and James's functional graphical lasso (equation 9):
`tr(S*Theta) - logdet(Theta) + lambda * sum(norm(Theta[j,k]) for j != k)`.
The sum is over **ordered** off-diagonal blocks. Entire diagonal blocks are
unpenalized. `block_sizes` is a vector, or a common block size. Scores have
shape `n × p × M`; dense data first use marginal FPCA.

ADMM uses `ProximalBase.LogDetLoss` for the likelihood proximal step. `rho`
is relative to squared average variance. No ridge is added to the data.
Nonconvergence is recorded and warned about. The returned positive-definite
iterate and sparse consensus are separate because finite ADMM tolerance does
not make them identical. This solves the paper's objective with a different
numerical algorithm from its block coordinate descent.

Solver defaults retain FGM's ADMM accuracy and iteration choices; they are
numerical controls, rather than defaults from the paper. Dense inputs are
complete FDA multivariate grid data. FDA's marginal FPCA uses trapezoidal
quadrature, `n-1` eigenvalues, and its mean-alignment sign rule. Explicit
`max_k=ncomponents` preserves the requested FGM truncation. All resolved FPCA
and graph tuning choices are stored on the result and shown by `show`.
"""
function functional_glasso(S::AbstractMatrix{<:Real}, block_sizes, lambda::Real;
                           rho::Real=1, maxiter::Integer=10000,
                           abstol::Real=1e-7, reltol::Real=1e-6,
                           threshold::Real=1e-6)
    Base.require_one_based_indexing(S)
    d = LinearAlgebra.checksquare(S)
    d > 0 || throw(ArgumentError("covariance cannot be empty"))
    all(isfinite, S) && isapprox(S, transpose(S); atol=1e-12, rtol=1e-10) ||
        throw(ArgumentError("covariance must be finite and symmetric"))
    all(x -> isfinite(x) && x >= 0, (lambda, threshold)) ||
        throw(ArgumentError("lambda and threshold must be finite and nonnegative"))
    isfinite(rho) && rho > 0 && maxiter > 0 || throw(ArgumentError("rho and maxiter must be positive"))
    all(x -> isfinite(x) && x > 0, (abstol, reltol)) || throw(ArgumentError("tolerances must be positive"))
    bs, blocks = _fgm_blocks(d, block_sizes)
    C = Matrix{Float64}(Symmetric(S))
    mineigen = eigmin(Symmetric(C))
    mineigen >= -1e-10 * max(1, opnorm(C)) || throw(ArgumentError("covariance must be positive semidefinite"))
    for b in blocks
        isposdef(Symmetric(C[b,b])) || throw(ArgumentError("each diagonal covariance block must be positive definite; reduce truncation"))
    end
    lambda == 0 && !isposdef(Symmetric(C)) && throw(ArgumentError("lambda=0 requires positive-definite covariance"))
    if mineigen < 0 && lambda > 0
        # A feasible dual covariance lies on the segment from C to its block
        # diagonal. Its minimum eigenvalue has this concavity lower bound.
        offmax = maximum((norm(C[a,b]) for a in blocks for b in blocks if a != b); init=0.0)
        fraction = offmax == 0 ? 1.0 : min(1.0, lambda/offmax)
        mindiag = minimum(eigmin(Symmetric(C[b,b])) for b in blocks)
        (1-fraction)*mineigen + fraction*mindiag > 0 ||
            throw(ArgumentError("covariance is numerically indefinite at this penalty; rescale or increase lambda"))
    end
    scale = tr(C) / d
    r = Float64(rho) * scale^2
    isfinite(scale) && scale > 0 && isfinite(r) && r > 0 && isfinite(inv(r)) ||
        throw(ArgumentError("covariance scale and rho exceed Float64 range; rescale the covariance"))
    X = Matrix{Float64}(I, d, d) / scale
    Z = copy(X)
    U = zeros(d,d)
    old = similar(Z)
    input = similar(Z)
    loss = ProximalBase.LogDetLoss(C)
    converged = false
    primal = dual = Inf
    iter = 0
    for iteration in 1:maxiter
        iter = iteration
        input .= Z .- U
        ProximalBase.prox!(X, loss, input, inv(r))
        old .= Z
        Z .= X .+ U
        for j in eachindex(blocks), k in j+1:length(blocks)
            a, b = blocks[j], blocks[k]
            z = view(Z, a, b)
            nz = norm(z)
            z .*= nz == 0 ? 0 : max(0, 1 - lambda / (r * nz))
            Z[b,a] .= transpose(z)
        end
        U .+= X .- Z
        primal = norm(X-Z)
        dual = r * norm(Z-old)
        epspri = d*abstol/scale + reltol*max(norm(X),norm(Z))
        epsdual = d*abstol*scale + reltol*r*norm(U)
        if primal <= epspri && dual <= epsdual
            converged = true
            break
        end
    end
    kkt = _fgm_kkt(C, X, Z, blocks, lambda)
    obj = dot(C,X) - logdet(Symmetric(X)) + lambda*_fgm_penalty(X,blocks)
    info = FGLassoDiagnostic(iter, converged, primal, dual, kkt, obj)
    converged || @warn "functional_glasso did not converge" iterations=iter primal dual
    settings = (; lambda=Float64(lambda), rho=Float64(rho), maxiter=Int(maxiter),
        abstol=Float64(abstol), reltol=Float64(reltol), threshold=Float64(threshold),
        block_sizes=copy(bs))
    FGLassoResult(X,Z,C,bs,_fgm_graph(Z,blocks,threshold),info,nothing,settings)
end

function functional_glasso(scores::AbstractArray{<:Real,3}, lambda::Real; kwargs...)
    functional_glasso(score_covariance(scores), size(scores,3), lambda; kwargs...)
end

function functional_glasso(data::FDA.MultivariateFunctionalData, lambda::Real;
                           ncomponents::Integer, kwargs...)
    _dense_parts(data)
    fpca = FDA.fit(FDA.FPCA, data; ncomponents, max_k=ncomponents, basis=:marginal)
    fit = functional_glasso(FDA.score_array(fpca), lambda; kwargs...)
    FGLassoResult(fit.precision,fit.sparse_precision,fit.covariance,fit.block_sizes,
                 fit.adjacency,fit.diagnostics,fpca,fit.settings)
end

function Base.show(io::IO, fit::FGLassoResult)
    print(io, "FGLassoResult(settings=", fit.settings, ", diagnostics=",
        fit.diagnostics, ", fpca=", fit.fpca, ")")
end
