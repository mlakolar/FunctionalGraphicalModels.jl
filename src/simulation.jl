"""
    simulate_qiao2019(rng; n=100, p=50, ncomponents=5, ntime=100,
                       model=:ar2, noise_sd=0.5, random_diagonal_margin=0.1)

Generate the finite Fourier models in Section 5.1 of Qiao, Guo and James,
*Functional Graphical Models*. The returned named tuple contains noisy
`data`, latent Fourier `scores`, their node-major `precision`, the true
undirected `adjacency`, and the simulation `settings`.

Models are `:ar2` (diagonal 1, first two precision lags 0.4 and 0.2),
`:alternating` (AR(2) on alternating groups of ten nodes, with the other
groups isolated), and `:random` (independent undirected edges with probability
0.1 and precision weight 0.5). The paper leaves the random model diagonal
unspecified except for positive definiteness. Here it is
`0.5 * maximum_degree + random_diagonal_margin`, which guarantees strict
diagonal dominance. This choice is recorded in `settings`.

The package Fourier convention on `[0,1]` is `1, √2sin(2πt), √2cos(2πt),
√2sin(4πt), √2cos(4πt), …`. The paper specifies five Fourier functions,
without this ordering convention. `ncomponents != 5` is a generalization of
its generating design. A generated dataset alone is not a reproduction of
the paper's tuning, competitor comparisons, or repeated simulation results.
"""
function simulate_qiao2019(rng::Random.AbstractRNG; n::Integer=100, p::Integer=50,
                           ncomponents::Integer=5, ntime::Integer=100,
                           model::Symbol=:ar2, noise_sd::Real=0.5,
                           random_diagonal_margin::Real=0.1)
    n >= 1 && p >= 1 && ncomponents >= 1 ||
        throw(ArgumentError("n, p and ncomponents must be positive"))
    ntime >= 2 || throw(ArgumentError("ntime must be at least two"))
    model in (:ar2, :alternating, :random) ||
        throw(ArgumentError("model must be :ar2, :alternating or :random"))
    isfinite(noise_sd) && noise_sd >= 0 ||
        throw(ArgumentError("noise_sd must be finite and nonnegative"))
    isfinite(random_diagonal_margin) && random_diagonal_margin > 0 ||
        throw(ArgumentError("random_diagonal_margin must be finite and positive"))
    Q = Matrix{Float64}(I, p, p)
    if model == :random
        for j in 2:p, k in 1:j-1
            if Random.rand(rng) < 0.1
                Q[j, k] = Q[k, j] = 0.5
            end
        end
        maxdegree = maximum(sum(!iszero, view(Q, j, :)) - 1 for j in 1:p)
        diagonal = 0.5 * maxdegree + random_diagonal_margin
        for j in 1:p
            Q[j, j] = diagonal
        end
    else
        for j in 1:p, lag in 1:2
            k = j - lag
            k >= 1 || continue
            same_active_block = fld(j - 1, 10) == fld(k - 1, 10) &&
                                iseven(fld(j - 1, 10))
            if model == :ar2 || same_active_block
                Q[j, k] = Q[k, j] = lag == 1 ? 0.4 : 0.2
            end
        end
        diagonal = 1.0
    end
    adjacency = .!iszero.(Q)
    for j in 1:p
        adjacency[j, j] = false
    end
    precision = kron(Q, Matrix{Float64}(I, ncomponents, ncomponents))
    # Θ = U'U; rows Z/U' therefore have covariance Θ^(-1).
    U = cholesky(Symmetric(precision)).U
    X = Random.randn(rng, n, p * ncomponents) / transpose(U)
    scores = permutedims(reshape(X, n, ncomponents, p), (1, 3, 2))
    grid = collect(range(0.0, 1.0; length=ntime))
    fourier = ones(ntime, ncomponents)
    for m in 2:ncomponents
        frequency = fld(m, 2)
        fourier[:, m] .= sqrt(2) .* (iseven(m) ? sin.(2pi * frequency .* grid) :
                                                cos.(2pi * frequency .* grid))
    end
    values = Array{Float64,3}(undef, n, p, ntime)
    for j in 1:p
        values[:, j, :] .= scores[:, j, :] * transpose(fourier)
    end
    if noise_sd > 0
        values .+= noise_sd .* Random.randn(rng, n, p, ntime)
    end
    settings = (; n=Int(n), p=Int(p), ncomponents=Int(ncomponents), ntime=Int(ntime),
                 model, noise_sd=Float64(noise_sd), diagonal=Float64(diagonal),
                 random_diagonal_margin=Float64(random_diagonal_margin),
                 random_diagonal_rule="0.5 * maximum_degree + random_diagonal_margin",
                 fourier_convention="1, sqrt(2)sin(2pi*t), sqrt(2)cos(2pi*t), ...",
                 paper_coefficient_dimension=ncomponents == 5,
                 source="https://personal.lse.ac.uk/qiaox/qiao.links/research_papers/2017_FGM.pdf",
                 section="5.1", reproduction_status=:generating_design)
    (; data=DenseFunctionalData(values, grid), scores, precision, adjacency, settings)
end

"""
    graph_metrics(estimated, truth)

Compare symmetric binary adjacency matrices. Each undirected pair is counted
once and diagonals are ignored. Returns `tp`, `fp`, `tn`, `fn`, true and false
positive rates, precision, recall, F1, specificity, and accuracy. Any rate with
an empty denominator is `NaN` (for example, recall for an empty true graph).
"""
function graph_metrics(estimated::AbstractMatrix, truth::AbstractMatrix)
    size(estimated) == size(truth) || throw(DimensionMismatch("adjacency matrices must have equal size"))
    p = size(truth, 1)
    size(truth, 2) == p || throw(DimensionMismatch("adjacency matrices must be square"))
    issymmetric(estimated) && issymmetric(truth) ||
        throw(ArgumentError("adjacency matrices must be symmetric"))
    all(x -> x == 0 || x == 1, estimated) && all(x -> x == 0 || x == 1, truth) ||
        throw(ArgumentError("adjacency matrices must contain only zeros and ones"))
    tp = fp = tn = fn = 0
    for j in 2:p, k in 1:j-1
        if !iszero(truth[j, k])
            !iszero(estimated[j, k]) ? (tp += 1) : (fn += 1)
        else
            !iszero(estimated[j, k]) ? (fp += 1) : (tn += 1)
        end
    end
    ratio(a, b) = iszero(b) ? NaN : a / b
    tpr = ratio(tp, tp + fn)
    fpr = ratio(fp, fp + tn)
    precision = ratio(tp, tp + fp)
    f1 = ratio(2tp, 2tp + fp + fn)
    specificity = ratio(tn, tn + fp)
    accuracy = ratio(tp + tn, tp + fp + tn + fn)
    (; tp, fp, tn, fn, tpr, fpr, precision, recall=tpr, f1, specificity, accuracy)
end
