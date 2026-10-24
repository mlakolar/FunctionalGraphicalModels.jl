"""
    quadrature_weights(grid)

Trapezoidal integration weights on a finite, strictly increasing grid.
The weights sum to `last(grid) - first(grid)`.
"""
function quadrature_weights(grid::AbstractVector{<:Real})
    length(grid) >= 2 || throw(ArgumentError("the grid needs at least two points"))
    all(isfinite, grid) || throw(ArgumentError("the grid must be finite"))
    all(>(0), diff(grid)) || throw(ArgumentError("the grid must be strictly increasing"))
    x = float.(collect(grid))
    w = similar(x)
    w[1] = (x[2] - x[1]) / 2
    w[end] = (x[end] - x[end - 1]) / 2
    for t in 2:length(x)-1
        w[t] = (x[t + 1] - x[t - 1]) / 2
    end
    w
end

"""
    DenseFunctionalData(values, grid)

Multivariate curves observed on a common grid. `values[i,j,t]` is the
measurement for subject `i`, node `j`, and grid point `t`. Inputs are copied.
Missing values and irregular subject-specific grids require preprocessing.
"""
struct DenseFunctionalData{T<:AbstractFloat}
    values::Array{T,3}
    grid::Vector{T}
    DenseFunctionalData{T}(values::Array{T,3}, grid::Vector{T}) where {T<:AbstractFloat} =
        new{T}(values, grid)
end

function DenseFunctionalData(values::AbstractArray{<:Real,3}, grid::AbstractVector{<:Real})
    Base.require_one_based_indexing(values, grid)
    n, p, nt = size(values)
    n >= 1 && p >= 1 || throw(ArgumentError("at least one subject and node are required"))
    nt == length(grid) || throw(DimensionMismatch("the last data dimension must match the grid"))
    quadrature_weights(grid)
    all(isfinite, values) || throw(ArgumentError("curve values must be finite"))
    T = float(promote_type(eltype(values), eltype(grid)))
    DenseFunctionalData{T}(Array{T,3}(values), Vector{T}(grid))
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
sparse-data PACE.
"""
function smooth_bspline(data::DenseFunctionalData; nbasis::Integer,
                        degree::Integer=3, penalty::Real=0)
    degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
    nbasis >= degree + 1 || throw(ArgumentError("nbasis must be at least degree + 1"))
    nbasis <= length(data.grid) || throw(ArgumentError("nbasis cannot exceed the number of grid points"))
    isfinite(penalty) && penalty >= 0 || throw(ArgumentError("penalty must be finite and nonnegative"))
    B = _bspline_design(data.grid, nbasis, degree)
    w = quadrature_weights(data.grid)
    WB = w .* B
    A = transpose(B) * WB
    if penalty > 0 && nbasis > 2
        D = diff(diff(Matrix{eltype(B)}(I, nbasis, nbasis); dims=1); dims=1)
        A .+= penalty .* (transpose(D) * D)
    end
    F = cholesky(Symmetric(A))
    n, p, nt = size(data.values)
    out = similar(data.values)
    for j in 1:p, i in 1:n
        c = F \ (transpose(WB) * view(data.values, i, j, :))
        out[i, j, :] .= B * c
    end
    DenseFunctionalData(out, data.grid)
end
