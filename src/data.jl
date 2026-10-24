# FGM's graph and P-spline routines require complete common-grid observations.
# FDA owns the input representation, validation, copying and quadrature rule.
function _dense_parts(data::FDA.MultivariateFunctionalData)
    values = FDA.to_array(data)
    grid = copy((FDA.component(data, 1)::FDA.GridFunctionalData).grid)
    quadrature_weights(grid)
    all(isfinite, values) || throw(ArgumentError("curve values must be finite"))
    values, grid
end

function _dense_data(values::AbstractArray{<:Real,3}, grid::AbstractVector{<:Real};
                     ids=nothing, metadata::NamedTuple=(;))
    Base.require_one_based_indexing(values, grid)
    size(values, 2) >= 1 || throw(ArgumentError("at least one node is required"))
    FDA.MultivariateFunctionalData([
        FDA.GridFunctionalData(values[:, j, :], grid; ids, metadata)
        for j in axes(values, 2)])
end

"""
    score_covariance(scores; center=true, corrected=false)

Covariance of an `n × p × M` score array, with coordinates ordered by node:
`(node 1, components 1:M), (node 2, components 1:M), …`. The default uses
the likelihood normalization `1/n`; `corrected=true` uses `1/(n-1)`.
"""
function score_covariance(scores::AbstractArray{<:Real,3}; center::Bool=true,
                          corrected::Bool=false)
    Base.require_one_based_indexing(scores)
    n, p, m = size(scores)
    n >= 1 && p >= 1 && m >= 1 || throw(ArgumentError("all score dimensions must be positive"))
    corrected && n < 2 && throw(ArgumentError("corrected covariance requires two subjects"))
    all(isfinite, scores) || throw(ArgumentError("scores must be finite"))
    X = float.(reshape(permutedims(scores, (1, 3, 2)), n, p * m))
    if center
        X .-= sum(X; dims=1) ./ n
    end
    Matrix(Symmetric(transpose(X) * X / (n - Int(corrected))))
end

# Values on `grid` of the `nbasis` B-splines of the given degree whose knots are
# equally spaced on [first(grid), last(grid)] and extended by `degree` knots with
# the same spacing on each side. BSplineKit takes the left limit at the last
# knot, so with degree 0 the final grid point keeps its basis function.
function _bspline_design(grid::AbstractVector{T}, nbasis::Integer,
                         degree::Integer) where {T<:AbstractFloat}
    a, b = first(grid), last(grid)
    h = (b - a) / (nbasis - degree)
    knots = vcat(a .- h .* (degree:-1:1), range(a, b; length=nbasis - degree + 1),
                 b .+ h .* (1:degree))
    basis = BSplineKit.BSplineBasis(BSplineKit.BSplineOrder(degree + 1), knots;
                                    augment=Val(false))
    BSplineKit.collocation_matrix(basis, grid, BSplineKit.Derivative(0), Matrix{T})
end

"""
    smooth_bspline(data; nbasis, degree=3, penalty=0)

Smooth each common-grid curve in a basis of `nbasis` B-splines of the given
`degree`, with equally spaced knots on the grid interval (BSplineKit.jl). Fits
minimize weighted squared error plus `penalty` times the squared second
differences of the basis coefficients. This is dense-curve smoothing, not
sparse-data PACE. The input and output are FDA `MultivariateFunctionalData`
with complete `GridFunctionalData` nodes on one common grid. Smoothing settings
are stored in each output node's `metadata`; subject IDs are retained.

The defaults `degree=3` and `penalty=0` retain FGM's cubic, unpenalized
P-spline convention. The second-difference penalty is distinct from FDA's
derivative-penalized spline smoother.
"""
function smooth_bspline(data::FDA.MultivariateFunctionalData; nbasis::Integer,
                        degree::Integer=3, penalty::Real=0)
    values, grid = _dense_parts(data)
    degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
    nbasis >= degree + 1 || throw(ArgumentError("nbasis must be at least degree + 1"))
    nbasis <= length(grid) || throw(ArgumentError("nbasis cannot exceed the number of grid points"))
    isfinite(penalty) && penalty >= 0 || throw(ArgumentError("penalty must be finite and nonnegative"))
    B = _bspline_design(grid, nbasis, degree)
    w = quadrature_weights(grid)
    WB = w .* B
    A = transpose(B) * WB
    if penalty > 0 && nbasis > 2
        D = diff(diff(Matrix{eltype(B)}(I, nbasis, nbasis); dims=1); dims=1)
        A .+= penalty .* (transpose(D) * D)
    end
    F = cholesky(Symmetric(A))
    n, p, nt = size(values)
    out = similar(values)
    for j in 1:p, i in 1:n
        c = F \ (transpose(WB) * view(values, i, j, :))
        out[i, j, :] .= B * c
    end
    metadata = (; method=:pspline, nbasis=Int(nbasis), degree=Int(degree),
        penalty, penalty_order=2, weighting=:trapezoidal, grid=copy(grid))
    _dense_data(out, grid; ids=FDA.ids(data), metadata)
end
