"""
    FGMComponentSelection

Time-fold selection of spline dimension `nbasis` and marginal FPCA truncation
`ncomponents`. Rows of `losses` index `nbasis_grid`; columns index
`ncomponents_grid`. `fold_losses` additionally indexes validation folds.
`reasons` explains invalid candidates, whose losses are `Inf`. `smoothed` and
`fpca` are refitted to all observed times at the selected dimensions.
"""
struct FGMComponentSelection{D<:DenseFunctionalData,F<:FPCAResult}
    nbasis::Int
    ncomponents::Int
    nbasis_grid::Vector{Int}
    ncomponents_grid::Vector{Int}
    losses::Matrix{Float64}
    fold_losses::Array{Float64,3}
    fold_ids::Vector{Int}
    reasons::Matrix{String}
    smoothed::D
    fpca::F
end

function _fgm_dimension_grid(values, name)
    raw = collect(values)
    !isempty(raw) || throw(ArgumentError("$name cannot be empty"))
    all(x -> x isa Integer && x > 0, raw) ||
        throw(ArgumentError("$name must contain positive integers"))
    grid = Int.(raw)
    length(unique(grid)) == length(grid) || throw(ArgumentError("$name must not contain duplicates"))
    grid
end

function _fgm_time_folds(nt, fold_ids, folds, rng)
    if fold_ids === nothing
        folds isa Integer && 2 <= folds <= nt ||
            throw(ArgumentError("folds must be between 2 and the number of time points"))
        order = Random.randperm(rng, nt)
        ids = zeros(Int, nt)
        for (position, t) in enumerate(order)
            ids[t] = mod1(position, folds)
        end
        return ids
    end
    length(fold_ids) == nt || throw(DimensionMismatch("fold_ids must have one entry per time point"))
    all(x -> x isa Integer && x > 0, fold_ids) ||
        throw(ArgumentError("fold_ids must contain positive integer labels"))
    ids = Int.(collect(fold_ids))
    labels = sort(unique(ids))
    length(labels) >= 2 && labels == collect(1:length(labels)) ||
        throw(ArgumentError("fold_ids must use contiguous labels 1:J for at least two folds"))
    ids
end

# The spline support is the original complete interval, even if endpoints are
# held out. Only values at `training_times` enter the coefficient fit.
function _fgm_training_smoother(data::DenseFunctionalData, training_times;
                                nbasis::Integer, degree::Integer=3,
                                penalty::Real=0)
    n, p, nt = size(data.values)
    degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
    nbasis >= degree + 1 || throw(ArgumentError("nbasis must be at least degree + 1"))
    isfinite(penalty) && penalty >= 0 || throw(ArgumentError("penalty must be finite and nonnegative"))
    raw = collect(training_times)
    all(t -> t isa Integer && 1 <= t <= nt, raw) ||
        throw(ArgumentError("training_times must index the observed grid"))
    train = sort(Int.(raw))
    length(unique(train)) == length(train) || throw(ArgumentError("training_times must not repeat"))
    length(train) >= max(2, nbasis) ||
        throw(ArgumentError("the training fold has too few time points for the spline dimension"))
    B = _bspline_design(data.grid, nbasis, degree)
    Bt = B[train, :]
    rank(Bt) == nbasis || throw(ArgumentError("the training spline design is rank deficient"))
    weights = quadrature_weights(data.grid[train])
    WB = weights .* Bt
    gram = transpose(Bt) * WB
    if penalty > 0 && nbasis > 2
        D = diff(diff(Matrix{eltype(B)}(I, nbasis, nbasis); dims=1); dims=1)
        gram .+= penalty .* (transpose(D) * D)
    end
    factor = cholesky(Symmetric(gram))
    fitted = similar(data.values)
    for j in 1:p, i in 1:n
        coefficients = factor \ (transpose(WB) * data.values[i, j, train])
        fitted[i, j, :] .= B * coefficients
    end
    DenseFunctionalData(fitted, data.grid)
end

function _fgm_retained_rank_reason(fpca::FPCAResult, m)
    for j in axes(fpca.eigenvalues, 2)
        largest = fpca.eigenvalues[1, j]
        tolerance = max(size(fpca.scores, 1), length(fpca.grid)) * eps(eltype(fpca.scores)) * largest
        smallest = minimum(@view fpca.eigenvalues[1:m, j])
        if !(largest > 0 && smallest > tolerance)
            return "retained FPCA scores are rank deficient at node $j"
        end
    end
    ""
end

function _fgm_truncate_fpca(fpca::FPCAResult{T}, m) where {T}
    1 <= m <= size(fpca.scores, 3) || throw(ArgumentError("invalid FPCA truncation"))
    FPCAResult{T}(fpca.scores[:, :, 1:m], fpca.eigenfunctions[:, 1:m, :],
                  fpca.eigenvalues[1:m, :], copy(fpca.mean), copy(fpca.grid),
                  copy(fpca.weights))
end

"""
    FunctionalGraphicalModels._fgm_fold_prediction(data, training_times;
        nbasis, ncomponents, degree=3, penalty=0)

Internal cross-fitting primitive used by time-fold selection. Fits curves,
means and FPCA using only outcomes at `training_times`, then returns the full
grid `prediction`, `smoothed` training-curve fits, and their truncated `fpca`.
No outcomes at held-out times are used to estimate any of these objects.
"""
function _fgm_fold_prediction(data::DenseFunctionalData, training_times;
                              nbasis::Integer, ncomponents::Integer,
                              degree::Integer=3, penalty::Real=0)
    ncomponents <= nbasis || throw(ArgumentError("ncomponents must not exceed nbasis"))
    smoothed = _fgm_training_smoother(data, training_times; nbasis, degree, penalty)
    fpca = fit_fpca(smoothed; ncomponents)
    reason = _fgm_retained_rank_reason(fpca, ncomponents)
    isempty(reason) || throw(ArgumentError(reason))
    (; prediction=reconstruct(fpca), smoothed, fpca)
end

"""
    select_fgm_components(data; nbasis_grid=5:10, ncomponents_grid=1:6,
        fold_ids=nothing, folds=5, rng=Random.default_rng(), degree=3, penalty=0)

Select the spline dimension `L` and FPCA dimension `M` using the observed-time
cross-validation of Qiao, Guo and James, Section 3.3. Each fold holds out a
group of time points on every curve. Spline coefficients, the mean functions
and marginal FPCA are estimated from the remaining outcomes. The criterion
is mean squared prediction error over all held-out observations, after
`M`-component reconstruction. The chosen model is refitted at all times.

This implementation uses one common randomized partition of the shared grid
for all curves. Fold sizes differ by at most one. Pass `fold_ids` to specify
the partition; when supplied, `folds` is ignored. Spline knots retain the full
original interval even when its endpoints are held out. Spline fits use the
same trapezoidal weighted least-squares convention as `smooth_bspline`.
Integration on the full grid uses fitted predictions, never held-out outcomes.

Candidates with `M > L`, `M > n-1`, insufficient training points, or numerical
rank deficiency receive `Inf` and an explanation in `reasons`. At least one
candidate must be valid in every fold. Exact loss ties follow input-grid order.
This chooses `L` and `M`; it does not select the graph regularization penalty.

Source: <https://personal.lse.ac.uk/qiaox/qiao.links/research_papers/2017_FGM.pdf>,
Section 3.3.
"""
function select_fgm_components(data::DenseFunctionalData; nbasis_grid=5:10,
                               ncomponents_grid=1:6, fold_ids=nothing,
                               folds=5, rng=Random.default_rng(),
                               degree::Integer=3, penalty::Real=0)
    n, p, nt = size(data.values)
    n >= 2 || throw(ArgumentError("component selection requires two subjects"))
    degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
    isfinite(penalty) && penalty >= 0 || throw(ArgumentError("penalty must be finite and nonnegative"))
    Ls = _fgm_dimension_grid(nbasis_grid, "nbasis_grid")
    Ms = _fgm_dimension_grid(ncomponents_grid, "ncomponents_grid")
    ids = _fgm_time_folds(nt, fold_ids, folds, rng)
    nfolds = maximum(ids)
    fold_losses = fill(Inf, length(Ls), length(Ms), nfolds)
    reasons = fill("", length(Ls), length(Ms))
    for (li, L) in enumerate(Ls)
        valid = Int[]
        for (mi, M) in enumerate(Ms)
            reason = L < degree + 1 ? "nbasis must be at least degree + 1" :
                     M > L ? "ncomponents must not exceed nbasis" :
                     M > min(n-1, nt) ? "ncomponents exceeds the available subject or grid rank" : ""
            if isempty(reason)
                push!(valid, mi)
            else
                reasons[li, mi] = reason
            end
        end
        isempty(valid) && continue
        for fold in 1:nfolds
            train = findall(!=(fold), ids)
            test = findall(==(fold), ids)
            local fpca
            try
                smoothed = _fgm_training_smoother(data, train; nbasis=L, degree, penalty)
                fpca = fit_fpca(smoothed; ncomponents=maximum(Ms[valid]))
            catch err
                err isa InterruptException && rethrow()
                message = "fold $fold: " * sprint(showerror, err)
                for mi in valid
                    reasons[li, mi] *= isempty(reasons[li, mi]) ? message : "; " * message
                end
                continue
            end
            for mi in valid
                M = Ms[mi]
                reason = _fgm_retained_rank_reason(fpca, M)
                if !isempty(reason)
                    message = "fold $fold: " * reason
                    reasons[li, mi] *= isempty(reasons[li, mi]) ? message : "; " * message
                    continue
                end
                prediction = reconstruct(_fgm_truncate_fpca(fpca, M))
                error = prediction.values[:, :, test] .- data.values[:, :, test]
                fold_losses[li, mi, fold] = sum(abs2, error) / (n * p * length(test))
            end
        end
    end
    losses = fill(Inf, length(Ls), length(Ms))
    for li in eachindex(Ls), mi in eachindex(Ms)
        if all(isfinite, @view fold_losses[li, mi, :])
            losses[li, mi] = sum(fold_losses[li, mi, fold] * count(==(fold), ids)
                                 for fold in 1:nfolds) / nt
        end
    end
    best = argmin(losses)
    isfinite(losses[best]) ||
        throw(ArgumentError("all (nbasis,ncomponents) candidates failed: " *
                            join(unique(vec(reasons)), "; ")))
    L, M = Ls[best[1]], Ms[best[2]]
    refit = _fgm_fold_prediction(data, collect(1:nt); nbasis=L, ncomponents=M,
                                degree, penalty)
    FGMComponentSelection(L, M, Ls, Ms, losses, fold_losses, ids, reasons,
                           refit.smoothed, refit.fpca)
end
