"""
    FPCAResult

A fitted common-grid FPCA projection. `scores` is `n × p × M`,
`eigenfunctions` is `T × M × p`, `eigenvalues` is `M × p`, and `mean` is
`p × T`. Eigenfunctions are orthonormal under the stored quadrature weights.
"""
struct FPCAResult{T<:AbstractFloat}
    scores::Array{T,3}
    eigenfunctions::Array{T,3}
    eigenvalues::Matrix{T}
    mean::Matrix{T}
    grid::Vector{T}
    weights::Vector{T}
end

"""
    fit_fpca(data; ncomponents, basis=:marginal)

Weighted FPCA of densely observed curves, using trapezoidal integration.
`:marginal` estimates a separate basis for each node. `:pooled` estimates
the common basis from the average of the node covariance kernels, as in
partial separability estimation. In both cases `eigenvalues[m,j]` is the
empirical variance of component `m` at node `j`, with denominator `n`.
The fitted mean and basis are reused by [`transform_fpca`](@ref).
"""
function fit_fpca(data::DenseFunctionalData; ncomponents::Integer,
                   basis::Symbol=:marginal)
    n, p, nt = size(data.values)
    n >= 2 || throw(ArgumentError("FPCA requires at least two subjects"))
    basis in (:marginal, :pooled) || throw(ArgumentError("basis must be :marginal or :pooled"))
    limit = basis == :marginal ? min(nt, n - 1) : min(nt, p * (n - 1))
    1 <= ncomponents <= limit || throw(ArgumentError("ncomponents must be between 1 and $limit"))
    T = eltype(data.values)
    w = quadrature_weights(data.grid)
    sw = sqrt.(w)
    mu = dropdims(sum(data.values; dims=1); dims=1) ./ n
    centered = data.values .- reshape(mu, 1, p, nt)
    phi = Array{T,3}(undef, nt, ncomponents, p)
    function eigenbasis(C)
        e = eigen(Symmetric(C))
        inds = length(e.values):-1:(length(e.values) - ncomponents + 1)
        e.vectors[:, inds] ./ sw
    end
    if basis == :pooled
        C = zeros(T, nt, nt)
        for j in 1:p
            A = centered[:, j, :] .* transpose(sw)
            C .+= transpose(A) * A / (n * p)
        end
        common = eigenbasis(C)
        for j in 1:p
            phi[:, :, j] .= common
        end
    else
        for j in 1:p
            A = centered[:, j, :] .* transpose(sw)
            phi[:, :, j] .= eigenbasis(transpose(A) * A / n)
        end
    end
    scores = Array{T,3}(undef, n, p, ncomponents)
    eigvals = Matrix{T}(undef, ncomponents, p)
    for j in 1:p
        scores[:, j, :] .= centered[:, j, :] * (w .* phi[:, :, j])
        eigvals[:, j] .= vec(sum(abs2, scores[:, j, :]; dims=1)) ./ n
    end
    FPCAResult{T}(scores, phi, eigvals, mu, copy(data.grid), w)
end

"""
    transform_fpca(fit, data)

Project new subjects using the fitted means and eigenfunctions. The node
count and grid must match the training data. Returns an `n × p × M` array.
"""
function transform_fpca(fit::FPCAResult, data::DenseFunctionalData)
    data.grid == fit.grid || throw(ArgumentError("new data must use the fitted grid"))
    n, p, _ = size(data.values)
    p == size(fit.mean, 1) || throw(DimensionMismatch("new data has the wrong number of nodes"))
    m = size(fit.eigenfunctions, 2)
    T = promote_type(eltype(data.values), eltype(fit.scores))
    scores = Array{T,3}(undef, n, p, m)
    for j in 1:p
        C = data.values[:, j, :] .- transpose(fit.mean[j, :])
        scores[:, j, :] .= C * (fit.weights .* fit.eigenfunctions[:, :, j])
    end
    scores
end

"""
    reconstruct(fit, scores=fit.scores)

Reconstruct common-grid curves from an `n × p × M` score array and the
fitted mean. Returns [`DenseFunctionalData`](@ref).
"""
function reconstruct(fit::FPCAResult, scores::AbstractArray{<:Real,3}=fit.scores)
    n, p, m = size(scores)
    (p, m) == (size(fit.mean, 1), size(fit.eigenfunctions, 2)) ||
        throw(DimensionMismatch("scores have the wrong node or component count"))
    all(isfinite, scores) || throw(ArgumentError("scores must be finite"))
    T = promote_type(float(eltype(scores)), eltype(fit.scores))
    values = Array{T,3}(undef, n, p, length(fit.grid))
    for j in 1:p
        values[:, j, :] .= scores[:, j, :] * transpose(fit.eigenfunctions[:, :, j]) .+
                            transpose(fit.mean[j, :])
    end
    DenseFunctionalData(values, fit.grid)
end
