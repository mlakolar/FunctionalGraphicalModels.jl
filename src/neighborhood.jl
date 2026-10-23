"""
    NeighborhoodDiagnostic

Optimization diagnostics for one target node. `kkt_residual` is the largest
Frobenius norm of a group KKT violation for the objective in EJS (2024), (13).
`converged` requires this residual and the ADMM feasibility residual to meet
the requested tolerances. Iteration limits are reported without claiming convergence.
"""
struct NeighborhoodDiagnostic
    iterations::Int
    converged::Bool
    objective::Float64
    kkt_residual::Float64
    primal_residual::Float64
    dual_residual::Float64
end

"""
    NeighborhoodResult

Functional neighborhood selection fit. `directed[j,k]` means node `k` was
selected in the regression for target `j`. `adjacency` combines these decisions
with the requested AND/OR rule. `coefficients[j]` is a `((p-1)M) × M` matrix:
successive `M × M` row blocks correspond to source nodes in increasing order,
omitting `j`, and are the transposes of the paper's regression matrices.

`target_bases` stores the sampled target FPCA bases for a dense-data fit, and is
`nothing` for scores supplied by the caller. All quantities use `Float64`.
"""
struct NeighborhoodResult
    adjacency::BitMatrix
    directed::BitMatrix
    coefficients::Vector{Matrix{Float64}}
    lambdas::Vector{Float64}
    thresholds::Vector{Float64}
    diagnostics::Vector{NeighborhoodDiagnostic}
    score_means::Vector{Matrix{Float64}}
    target_bases::Union{Nothing,Array{Float64,3}}
end

function _neighborhood_parameters(value, p, name)
    values = value isa Real ? fill(Float64(value), p) : Float64.(collect(value))
    length(values) == p || throw(DimensionMismatch("$name must have one value per node"))
    all(x -> isfinite(x) && x >= 0, values) ||
        throw(ArgumentError("$name must be finite and nonnegative"))
    values
end

function _neighborhood_validate_scores(scores)
    Base.require_one_based_indexing(scores)
    ndims(scores) == 3 || throw(DimensionMismatch("scores must have dimensions n × p × M"))
    n, p, m = size(scores)
    n >= 2 || throw(ArgumentError("at least two observations are required"))
    p >= 2 || throw(ArgumentError("at least two nodes are required"))
    m >= 1 || throw(ArgumentError("at least one component is required"))
    all(isfinite, scores) || throw(ArgumentError("scores must be finite"))
    n, p, m
end

function _neighborhood_design(scores, target; center=true)
    n, p, m = size(scores)
    means = zeros(p, m)
    if center
        for k in 1:p, l in 1:m
            means[k,l] = sum(view(scores, :, k, l)) / n
        end
    end
    X = Matrix{Float64}(undef, n, (p-1)*m)
    Y = Matrix{Float64}(undef, n, m)
    q = 0
    for k in 1:p
        if k == target
            for l in 1:m, i in 1:n
                Y[i,l] = scores[i,k,l] - means[k,l]
            end
        else
            for l in 1:m, i in 1:n
                X[i,q+l] = scores[i,k,l] - means[k,l]
            end
            q += m
        end
    end
    X, Y, means
end

function _neighborhood_groups(m, number_of_groups)
    q = m * number_of_groups
    [[row + (column-1)*q for column in 1:m for row in ((k-1)*m+1):(k*m)]
     for k in 1:number_of_groups]
end

function _neighborhood_kkt(B, gradient, groups, lambda)
    violation = 0.0
    for group in groups
        b = view(B, group)
        g = view(gradient, group)
        bnorm = norm(b)
        residual = iszero(bnorm) ? max(norm(g)-lambda, 0.0) :
                   norm(g .+ (lambda / bnorm) .* b)
        violation = max(violation, residual)
    end
    violation
end

# Solve (13) without materializing the Kronecker-expanded multiresponse design.
# Group weights are one: the default sqrt(group size) in GroupNormL2 would
# otherwise multiply the paper's lambda by M.
function _neighborhood_group_lasso(X, Y, m, lambda;
                                   maxiter=10000, atol=1e-7, rtol=1e-6, rho=1.0)
    maxiter isa Integer && maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    isfinite(atol) && atol > 0 || throw(ArgumentError("atol must be finite and positive"))
    isfinite(rtol) && rtol >= 0 || throw(ArgumentError("rtol must be finite and nonnegative"))
    rho = Float64(rho)
    isfinite(rho) && rho > 0 || throw(ArgumentError("rho must be finite and positive in Float64"))
    all(isfinite,X) && all(isfinite,Y) ||
        throw(ArgumentError("centered scores overflow Float64; rescale the supplied scores"))
    n, q = size(X)
    groups = _neighborhood_groups(m, div(q,m))
    penalty = ProximalBase.GroupNormL2(lambda, groups, ones(length(groups)))
    gram = transpose(X) * X / n
    cross = transpose(X) * Y / n
    all(isfinite,gram) && all(isfinite,cross) ||
        throw(ArgumentError("score cross-products overflow Float64; rescale the supplied scores"))
    scale = maximum(group -> norm(view(cross, group)), groups)
    tolerance = atol + rtol * max(scale, lambda)
    B = zeros(q, m)
    if scale <= lambda
        diagnostic = NeighborhoodDiagnostic(0, true, sum(abs2,Y)/(2n),
                                            0.0, 0.0, 0.0)
        return B, diagnostic
    end
    system = gram + rho*I
    all(isfinite,system) || throw(ArgumentError("rho plus score curvature overflows Float64"))
    factor = cholesky(Symmetric(system))
    smooth = copy(B)
    dual = copy(B)
    previous = copy(B)
    tmp = copy(B)
    gradient = copy(B)
    primal_residual = dual_residual = kkt = Inf
    converged = false
    iterations = 0
    for iteration in 1:maxiter
        iterations = iteration
        previous .= B
        @. tmp = cross + rho*(B-dual)
        smooth .= factor \ tmp
        @. tmp = smooth + dual
        ProximalBase.prox!(B, penalty, tmp, inv(rho))
        @. dual += smooth - B
        all(isfinite,B) && all(isfinite,dual) ||
            error("neighborhood iterations overflowed Float64; rescale the supplied scores")
        primal_residual = norm(smooth-B)
        dual_residual = rho * norm(B-previous)
        mul!(gradient, gram, B)
        gradient .-= cross
        kkt = _neighborhood_kkt(B, gradient, groups, lambda)
        feasible_tolerance = sqrt(length(B))*atol + rtol*max(norm(smooth),norm(B))
        if primal_residual <= feasible_tolerance && kkt <= tolerance
            converged = true
            break
        end
    end
    objective = sum(abs2, Y-X*B)/(2n) + penalty(B)
    B, NeighborhoodDiagnostic(iterations, converged, objective, kkt,
                              primal_residual, dual_residual)
end

function _neighborhood_assemble(coefficients, diagnostics, means, lambdas,
                                thresholds, p, m, rule, target_bases)
    rule in (:and, :or) || throw(ArgumentError("rule must be :and or :or"))
    directed = falses(p,p)
    for j in 1:p
        group = 0
        for k in 1:p
            k == j && continue
            group += 1
            block = view(coefficients[j], ((group-1)*m+1):(group*m), :)
            directed[j,k] = norm(block) > thresholds[j]
        end
    end
    adjacency = rule == :and ? directed .& transpose(directed) :
                              directed .| transpose(directed)
    NeighborhoodResult(adjacency, directed, coefficients, lambdas, thresholds,
                       diagnostics, means, target_bases)
end

function _functional_neighborhood(target_scores, lambdas, thresholds;
                                   rule=:or, center=true, target_bases=nothing,
                                   kwargs...)
    _, p, m = _neighborhood_validate_scores(target_scores[1])
    rule in (:and, :or) || throw(ArgumentError("rule must be :and or :or"))
    coefficients = Matrix{Float64}[]
    diagnostics = NeighborhoodDiagnostic[]
    means = Matrix{Float64}[]
    for j in 1:p
        size(target_scores[j]) == size(target_scores[1]) ||
            throw(DimensionMismatch("all target score arrays must have the same size"))
        _neighborhood_validate_scores(target_scores[j])
        X,Y,mu = _neighborhood_design(target_scores[j], j; center)
        B, diagnostic = _neighborhood_group_lasso(X,Y,m,lambdas[j]; kwargs...)
        push!(coefficients, B)
        push!(diagnostics, diagnostic)
        push!(means, mu)
    end
    _neighborhood_assemble(coefficients, diagnostics, means, lambdas, thresholds,
                           p, m, rule, target_bases)
end

"""
    functional_neighborhood(scores, lambda; threshold=0, rule=:or, center=true,
                            maxiter=10000, atol=1e-7, rtol=1e-6, rho=1)
    functional_neighborhood(data::DenseFunctionalData, lambda; ncomponents, ...)

Estimate the graph with Zhao, Zhai, Wang and Kolar's EJS (2024) neighborhood
selection. Scores have dimensions `n × p × M` and must use a common orthonormal
basis for all nodes. This low-level method cannot check the basis provenance.
For each target, solve the exact objective (13):

`||Y - X*B||_F^2/(2n) + lambda * sum(norm(B_k))`.

`threshold` is the Frobenius-norm threshold in (14), and `rule` is `:and` or
`:or`. A scalar `lambda`/`threshold` applies to all targets; a vector specifies
one value per target. Scores are centered by default. Inspect every
`result.diagnostics[j].converged` before using a fit.

For dense curves, each target's marginal FPCA basis is estimated, and **all**
nodes are projected onto that target's basis (Algorithm 1 and Section 2.4).
This differs from regressing the separate marginal FPCA scores of each node.
Tuning is user-supplied; see [`select_neighborhood_scv`](@ref) for Algorithm 2.

Source: <https://arxiv.org/html/2105.02487v3#S2.SS3>.
"""
function functional_neighborhood(scores::AbstractArray{<:Real,3}, lambda;
                                  threshold=0.0, kwargs...)
    _, p, _ = _neighborhood_validate_scores(scores)
    lambdas = _neighborhood_parameters(lambda,p,"lambda")
    thresholds = _neighborhood_parameters(threshold,p,"threshold")
    _functional_neighborhood(fill(scores,p),lambdas,thresholds; kwargs...)
end

function _neighborhood_target_scores(data::DenseFunctionalData, ncomponents)
    fpca = fit_fpca(data; ncomponents, basis=:marginal)
    n,p,t = size(data.values)
    m = size(fpca.eigenfunctions,2)
    target_scores = [zeros(n,p,m) for _ in 1:p]
    for j in 1:p, k in 1:p, l in 1:m, i in 1:n
        target_scores[j][i,k,l] = sum((data.values[i,k,s]-fpca.mean[k,s]) *
            fpca.weights[s] * fpca.eigenfunctions[s,l,j] for s in 1:t)
    end
    target_scores, Array{Float64,3}(fpca.eigenfunctions)
end

function functional_neighborhood(data::DenseFunctionalData, lambda;
                                  ncomponents, threshold=0.0, kwargs...)
    target_scores, bases = _neighborhood_target_scores(data,ncomponents)
    p = size(data.values,2)
    lambdas = _neighborhood_parameters(lambda,p,"lambda")
    thresholds = _neighborhood_parameters(threshold,p,"threshold")
    _functional_neighborhood(target_scores,lambdas,thresholds; target_bases=bases,kwargs...)
end

"""
    neighborhood_lambda_max(scores, target; center=true)

Smallest group-lasso penalty for which the zero regression is optimal for
`target` (EJS Proposition 1): `maximum_k norm(X_k'Y)/n`. Scores must use the
target's common basis and have dimensions `n × p × M`.
"""
function neighborhood_lambda_max(scores::AbstractArray{<:Real,3}, target::Integer;
                                   center=true)
    n,p,m = _neighborhood_validate_scores(scores)
    1 <= target <= p || throw(ArgumentError("target must index a node"))
    X,Y,_ = _neighborhood_design(scores,target; center)
    cross = transpose(X)*Y/n
    maximum(group -> norm(view(cross,group)), _neighborhood_groups(m,p-1))
end

"""
    NeighborhoodSCVResult

Selective cross-validation result. `fit` contains the selected full-data
group-lasso fits. `criteria[a,b,j]` is the mean SCV-RSS for penalty candidate
`lambdas[a]`, threshold candidate `thresholds[b]`, and target `j`.
`candidate_converged[a,j]` records optimization success. Unconverged candidates
receive infinite criteria and cannot be selected.
"""
struct NeighborhoodSCVResult
    fit::NeighborhoodResult
    criteria::Array{Float64,3}
    lambdas::Vector{Float64}
    thresholds::Vector{Float64}
    fold_ids::Vector{Int}
    candidate_converged::BitMatrix
end

function _neighborhood_scv_rss(X, Y, selected_columns, fold_ids)
    criterion = 0.0
    folds = unique(fold_ids)
    for fold in folds
        train = findall(!=(fold),fold_ids)
        test = findall(==(fold),fold_ids)
        Ytrain = Y[train,:]
        Ytest = Y[test,:]
        if isempty(selected_columns)
            criterion += sum(abs2,Ytest)
        else
            Xtrain = X[train,selected_columns]
            # A rank-deficient selected model has multiple minimizers; pinv
            # chooses the minimum-norm unpenalized least-squares refit in (18).
            refit = pinv(Xtrain)*Ytrain
            criterion += sum(abs2,Ytest-X[test,selected_columns]*refit)
        end
        # Count selected NODE groups, not scalar coefficients (equation 19).
        criterion += log(length(test)) * length(selected_columns)/size(Y,2)
    end
    criterion/length(folds)
end

function _select_neighborhood_scv(target_scores, lambda_grid, threshold_grid;
                                   folds=5, fold_ids=nothing, rule=:or,
                                   center=true, target_bases=nothing, kwargs...)
    n,p,m = _neighborhood_validate_scores(target_scores[1])
    lambdas = Float64.(collect(lambda_grid))
    thresholds = Float64.(collect(threshold_grid))
    !isempty(lambdas) && !isempty(thresholds) ||
        throw(ArgumentError("candidate grids must be nonempty"))
    all(x -> isfinite(x) && x >= 0,lambdas) ||
        throw(ArgumentError("lambda candidates must be finite and nonnegative"))
    all(x -> isfinite(x) && x >= 0,thresholds) ||
        throw(ArgumentError("threshold candidates must be finite and nonnegative"))
    rule in (:and,:or) || throw(ArgumentError("rule must be :and or :or"))
    if fold_ids === nothing
        folds isa Integer && 2 <= folds <= n ||
            throw(ArgumentError("folds must be an integer between 2 and n"))
        fold_ids = mod1.(1:n,folds)
    else
        fold_ids = Int.(collect(fold_ids))
        length(fold_ids) == n || throw(DimensionMismatch("fold_ids must have length n"))
        length(unique(fold_ids)) >= 2 || throw(ArgumentError("at least two folds are required"))
    end
    criteria = fill(Inf,length(lambdas),length(thresholds),p)
    candidate_converged = falses(length(lambdas),p)
    coefficients = Matrix{Float64}[]
    diagnostics = NeighborhoodDiagnostic[]
    means = Matrix{Float64}[]
    selected_lambdas = zeros(p)
    selected_thresholds = zeros(p)
    for j in 1:p
        _neighborhood_validate_scores(target_scores[j])
        X,Y,mu = _neighborhood_design(target_scores[j],j; center)
        best_value = Inf
        best_B = nothing
        best_diagnostic = nothing
        for (a,lambda) in pairs(lambdas)
            B,diagnostic = _neighborhood_group_lasso(X,Y,m,lambda; kwargs...)
            candidate_converged[a,j] = diagnostic.converged
            diagnostic.converged || continue
            block_norms = [norm(view(B,((k-1)*m+1):(k*m),:)) for k in 1:(p-1)]
            for (b,threshold) in pairs(thresholds)
                selected_groups = findall(>(threshold),block_norms)
                columns = reduce(vcat,[collect(((k-1)*m+1):(k*m)) for k in selected_groups]; init=Int[])
                score = _neighborhood_scv_rss(X,Y,columns,fold_ids)
                criteria[a,b,j] = score
                if score < best_value
                    best_value = score
                    best_B = B
                    best_diagnostic = diagnostic
                    selected_lambdas[j] = lambda
                    selected_thresholds[j] = threshold
                end
            end
        end
        best_B === nothing && error("no converged SCV candidate for target $j; increase maxiter or relax tolerances")
        push!(coefficients,best_B)
        push!(diagnostics,best_diagnostic)
        push!(means,mu)
    end
    fit = _neighborhood_assemble(coefficients,diagnostics,means,selected_lambdas,
                                 selected_thresholds,p,m,rule,target_bases)
    NeighborhoodSCVResult(fit,criteria,lambdas,thresholds,fold_ids,candidate_converged)
end

"""
    select_neighborhood_scv(scores, lambda_grid, threshold_grid; folds=5,
                            fold_ids=nothing, rule=:or, kwargs...)
    select_neighborhood_scv(data::DenseFunctionalData, lambda_grid, threshold_grid;
                            ncomponents, kwargs...)

Select target-specific `(lambda, threshold)` pairs by EJS (2024), Algorithm 2.
For each candidate penalty, select node groups using all observations. Refit
unpenalized least squares on each fold's training observations with the support
held fixed. Average `RSS + log(n_test)*number_of_selected_nodes` across folds
(equations 18–19). This intentionally uses full-data support selection, as in
the paper; it is a model-selection criterion, not an independent estimate of
out-of-sample risk. Scores are centered once on the full data when `center=true`;
the fold refits omit intercepts, as specified in equation (18).

`fold_ids` specifies a reproducible partition. By default observations are
assigned cyclically to `folds` folds without consuming random numbers. Ties
choose the first candidate in the supplied grids. The dense overload estimates
the target FPCA bases once before SCV; basis dimension is user-supplied and is
not selected by this function.
"""
function select_neighborhood_scv(scores::AbstractArray{<:Real,3},lambda_grid,
                                  threshold_grid; kwargs...)
    _,p,_ = _neighborhood_validate_scores(scores)
    _select_neighborhood_scv(fill(scores,p),lambda_grid,threshold_grid; kwargs...)
end

function select_neighborhood_scv(data::DenseFunctionalData,lambda_grid,
                                  threshold_grid; ncomponents,kwargs...)
    target_scores,bases = _neighborhood_target_scores(data,ncomponents)
    _select_neighborhood_scv(target_scores,lambda_grid,threshold_grid;
                             target_bases=bases,kwargs...)
end

