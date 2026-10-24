# Run only through m4_compare.jl, which gives each package revision its own process.
using FunctionalGraphicalModels
using LinearAlgebra
using Random
using Serialization
using Pkg

const FGM = FunctionalGraphicalModels
const MIGRATED = !isdefined(FGM, :DenseFunctionalData)
if MIGRATED
    @eval import FunctionalDataAnalysis as FDA
end

function dense_data(values, grid)
    MIGRATED || return FGM.DenseFunctionalData(values, grid)
    FDA.MultivariateFunctionalData([
        FDA.GridFunctionalData(values[:, j, :], grid) for j in axes(values, 2)])
end

function sparse_data(times, values)
    MIGRATED || return FGM.SparseFunctionalData(times, values)
    FDA.MultivariateFunctionalData([
        FDA.IrregularFunctionalData(times[:, j], values[:, j]) for j in axes(times, 2)])
end

function fpca_parts(f)
    if MIGRATED
        return (; scores=FDA.score_array(f),
            bases=cat((b.values for b in f.bases)...; dims=3),
            eigenvalues=hcat(f.eigenvalues...),
            means=reduce(vcat, (mu.values for mu in f.means)),
            weights=copy(first(f.bases).weights), settings=repr(f.settings))
    end
    (; scores=copy(f.scores), bases=copy(f.eigenfunctions),
        eigenvalues=copy(f.eigenvalues), means=copy(f.mean),
        weights=copy(f.weights), settings="FGM 6f0eac3: trapezoidal quadrature; covariance denominator n; unmodified LAPACK eigenvector signs")
end

function dense_reconstruction_fixture()
    n, p, q = 60, 3, 17
    grid = collect(range(0, 1; length=q))
    values = zeros(n, p, q)
    for i in 1:n, t in 1:q
        phi = sqrt(2) * cos(2pi * grid[t])
        a, b, c, d, e, f = [sqrt(2) * sin(2pi * l * i / n) for l in 1:6]
        values[i, 1, t] = 5 + a + 0.25b * phi
        values[i, 2, t] = -2 + 0.4a + c + 0.25d * phi
        values[i, 3, t] = 1 + e + 0.25f * phi
    end
    dense_data(values, grid)
end

function snapshot(fit, parts; fixture, settings, coordinates, precision,
        covariance, adjacency=fit.adjacency, scores=parts.scores,
        coefficients=nothing, noise_variance=nothing)
    (; fixture, settings=repr(settings), fpca_settings=parts.settings,
        coordinates, scores, covariance=Array(covariance),
        precision=precision === nothing ? nothing : Array(precision),
        coefficients=coefficients === nothing ? nothing : Array(coefficients),
        adjacency=Array(adjacency), bases=parts.bases, weights=parts.weights,
        eigenvalues=parts.eigenvalues, nsubjects=size(scores, 1), noise_variance)
end

function comparison_fixtures()
    rows = Dict{String,Any}()
    # Exact covariance fixture from test/glasso.jl. This checks the solver alone.
    C = [2.0 0 0.3 0; 0 1.5 0 0.2; 0.3 0 1.2 0; 0 0.2 0 1.0]
    settings = (; lambda=0.07, block_sizes=[2, 2], rho=1.0, maxiter=10000,
        abstol=1e-7, reltol=1e-6, threshold=1e-6)
    fit = functional_glasso(C, settings.block_sizes, settings.lambda;
        rho=settings.rho, maxiter=settings.maxiter, abstol=settings.abstol,
        reltol=settings.reltol, threshold=settings.threshold)
    fit.diagnostics.converged || error("glasso covariance fixture did not converge")
    rows["glasso_covariance"] = (; fixture="test/glasso.jl: node-wise orthogonal rotations",
        settings=repr(settings), fpca_settings="not applicable", coordinates=:direct,
        scores=nothing, covariance=copy(fit.covariance), precision=copy(fit.precision),
        coefficients=nothing,
        adjacency=Array(fit.adjacency), bases=nothing, weights=nothing,
        eigenvalues=nothing, nsubjects=0, noise_variance=nothing)

    # The dense reconstruction fixture has full-rank score covariance. The
    # test/fpca.jl fixture shares two latent scores across nodes and is singular.
    data = dense_reconstruction_fixture()
    settings = (; lambda=0.07, ncomponents=2, basis=:marginal, rho=1.0,
        maxiter=10000, abstol=1e-7, reltol=1e-6, threshold=1e-6)
    fit = functional_glasso(data, settings.lambda; ncomponents=settings.ncomponents,
        rho=settings.rho, maxiter=settings.maxiter, abstol=settings.abstol,
        reltol=settings.reltol, threshold=settings.threshold)
    fit.diagnostics.converged || error("glasso dense fixture did not converge")
    rows["glasso_dense"] = snapshot(fit, fpca_parts(fit.fpca);
        fixture="test/doubly.jl: full-rank dense reconstruction, through glasso entry point", settings, coordinates=:node_blocks,
        covariance=fit.covariance, precision=fit.precision)

    # Exact target-specific projection fixture from test/neighborhood.jl.
    n, p, q = 30, 3, 41
    grid = collect(range(0, 1; length=q))
    rng = MersenneTwister(301)
    a, b = randn(rng, n, p), randn(rng, n, p)
    values = zeros(n, p, q)
    for i in 1:n, j in 1:p, s in 1:q
        values[i, j, s] = 2a[i, j] * sinpi(j * grid[s]) + 0.3b[i, j] * cospi(grid[s])
    end
    data = dense_data(values, grid)
    settings = (; lambda=0.08, ncomponents=2, basis=:marginal, threshold=0.0,
        rule=:or, center=true, maxiter=10000, atol=1e-7, rtol=1e-6, rho=1.0,
        seed=301, rng="MersenneTwister")
    fit = functional_neighborhood(data, settings.lambda; ncomponents=settings.ncomponents,
        threshold=settings.threshold, rule=settings.rule, center=settings.center,
        maxiter=settings.maxiter, atol=settings.atol, rtol=settings.rtol, rho=settings.rho)
    all(d -> d.converged, fit.diagnostics) || error("neighborhood fixture did not converge")
    # The migrated estimator retains its own FDA fit. The base does not retain
    # FPCA means/eigenvalues, so only there recover them with its identical fit.
    parts = fpca_parts(MIGRATED ? fit.fpca :
        FGM.fit_fpca(data; ncomponents=settings.ncomponents, basis=settings.basis))
    parts = merge(parts, (; bases=copy(fit.target_bases)))
    projected = zeros(n, p, 2, p)
    covariance = zeros(2p, 2p, p)
    for target in 1:p
        for node in 1:p
            projected[:, node, :, target] .=
                (values[:, node, :] .- transpose(parts.means[node, :])) *
                (parts.weights .* parts.bases[:, :, target])
        end
        covariance[:, :, target] .= score_covariance(projected[:, :, :, target])
    end
    rows["neighborhood_dense"] = snapshot(fit, parts;
        fixture="test/neighborhood.jl: target-specific FPCA projection", settings,
        coordinates=:target_blocks, covariance, precision=nothing, scores=projected,
        coefficients=cat(fit.coefficients...; dims=3))

    # Exact dense reconstruction fixture from test/doubly.jl.
    settings = (; lambda=0.0, ncomponents=2, basis=:marginal, threshold=1e-6)
    fit = fit_doubly(dense_reconstruction_fixture(); lambda=settings.lambda,
        ncomponents=settings.ncomponents, threshold=settings.threshold)
    rows["doubly_dense"] = snapshot(fit, fpca_parts(fit.fpca);
        fixture="test/doubly.jl: dense reconstruction and CLIME", settings,
        coordinates=:pointwise, covariance=fit.covariances, precision=fit.precisions)

    # Exact positive-noise conditional-score fixture from test/doubly.jl.
    n = 60
    times = Matrix{Vector{Float64}}(undef, n, 1)
    values = similar(times)
    scores = [sqrt(2) * sin(2pi * i / n) for i in 1:n]
    for i in 1:n
        times[i, 1] = collect(range(0.02mod(i, 3), 1 - 0.01mod(i, 4);
            length=5 + mod(i, 3)))
        values[i, 1] = fill(scores[i], length(times[i, 1]))
    end
    grid = collect(range(0, 1; length=9))
    settings = (; lambda=0.1, grid, bandwidth=0.3, ncomponents=1,
        noise_variance=0.04, kernel=:gaussian, weighting=:subj, threshold=1e-6)
    fit = fit_doubly(sparse_data(times, values); lambda=settings.lambda, grid,
        bandwidth=settings.bandwidth, ncomponents=settings.ncomponents,
        noise_variance=settings.noise_variance, kernel=settings.kernel,
        threshold=settings.threshold)
    rows["doubly_sparse"] = snapshot(fit, fpca_parts(fit.fpca);
        fixture="test/doubly.jl: analytic conditional scores with supplied noise", settings,
        coordinates=:pointwise, covariance=fit.covariances, precision=fit.precisions)

    # Exact two-mode fixture from test/doubly.jl, with estimated noise.
    # Its strictly positive estimate avoids a singular zero-noise CE system.
    n = 40
    schedule = [0.02, 0.06, 0.14, 0.27, 0.44, 0.73, 0.98]
    times = Matrix{Vector{Float64}}(undef, n, 1)
    values = similar(times)
    for i in 1:n
        a, b = sqrt(2) * sin(2pi * i / n), 0.7sqrt(2) * cos(2pi * i / n)
        times[i, 1] = copy(schedule)
        values[i, 1] = a .+ b .* (2schedule .- 1)
    end
    settings = (; lambda=0.1, grid=collect(range(0, 1; length=17)),
        bandwidth=0.2, ncomponents=1, noise_variance=nothing,
        kernel=:gaussian, weighting=:subj, threshold=1e-6)
    fit = fit_doubly(sparse_data(times, values); lambda=settings.lambda,
        grid=settings.grid, bandwidth=settings.bandwidth,
        ncomponents=settings.ncomponents, kernel=settings.kernel,
        threshold=settings.threshold)
    noise = MIGRATED ? [FDA.component(fit.fpca, 1).noise_variance] :
        copy(fit.fpca.noise_variance)
    rows["doubly_sparse_estimated_noise"] = snapshot(fit, fpca_parts(fit.fpca);
        fixture="test/doubly.jl: two covariance modes, estimated noise", settings,
        coordinates=:pointwise, covariance=fit.covariances, precision=fit.precisions,
        noise_variance=noise)

    # Exact partially separable fixture from test/separable.jl.
    n, p, q = 60, 3, 21
    grid = collect(range(0, 1; length=q))
    values = zeros(n, p, q)
    for i in 1:n, t in 1:q
        a, b, c, d, e, f = [sqrt(2) * sin(2pi * l * i / n) for l in 1:6]
        phi = sqrt(2) * cos(2pi * grid[t])
        values[i, 1, t] = a + 0.3b * phi
        values[i, 2, t] = 0.5a + c + 0.3d * phi
        values[i, 3, t] = e + 0.3f * phi
    end
    settings = (; lambda=0.1, ncomponents=2, basis=:pooled, alpha=0.5,
        threshold=1e-6, rho=1.0, maxiter=10000, abstol=1e-8, reltol=1e-7)
    fit = fit_separable(dense_data(values, grid); lambda=settings.lambda,
        ncomponents=settings.ncomponents, alpha=settings.alpha, threshold=settings.threshold,
        rho=settings.rho, maxiter=settings.maxiter, abstol=settings.abstol, reltol=settings.reltol)
    fit.diagnostics.converged || error("separable fixture did not converge")
    parts = fpca_parts(fit.fpca)
    covariance = zeros(p, p, settings.ncomponents)
    for l in 1:settings.ncomponents
        A = parts.scores[:, :, l]
        A .-= sum(A; dims=1) / n
        covariance[:, :, l] .= transpose(A) * A / n
    end
    rows["separable_dense"] = snapshot(fit, parts;
        fixture="test/separable.jl: partially separable dense data", settings,
        coordinates=:components, covariance, precision=fit.precisions)

    # FDA PR #16 preserves the base's pooled rank range above n-1.
    n, p, q = 3, 2, 9
    grid = collect(range(0, 1; length=q))
    values = randn(MersenneTwister(416), n, p, q)
    settings = (; lambda=2.0, ncomponents=p*(n-1), basis=:pooled, alpha=0.5,
        threshold=1e-6, rho=1.0, maxiter=2000, abstol=1e-6, reltol=1e-5,
        seed=416, rng="MersenneTwister")
    fit = fit_separable(dense_data(values, grid); lambda=settings.lambda,
        ncomponents=settings.ncomponents, alpha=settings.alpha,
        threshold=settings.threshold, rho=settings.rho, maxiter=settings.maxiter,
        abstol=settings.abstol, reltol=settings.reltol)
    fit.diagnostics.converged || error("separable pooled-rank fixture did not converge")
    parts = fpca_parts(fit.fpca)
    covariance = zeros(p, p, settings.ncomponents)
    for l in 1:settings.ncomponents
        A = parts.scores[:, :, l]
        A .-= sum(A; dims=1) / n
        covariance[:, :, l] .= transpose(A) * A / n
    end
    rows["separable_pooled_rank"] = snapshot(fit, parts;
        fixture="test/separable.jl: pooled rank exceeds n-1", settings,
        coordinates=:components, covariance, precision=fit.precisions)
    rows
end

length(ARGS) == 1 || error("m4_compare_worker.jl requires one snapshot output path")
dependencies = Dict(string(info.name) => (; version=string(info.version),
    revision=info.git_revision, tree_sha1=string(info.tree_hash), source=info.git_source)
    for info in values(Pkg.dependencies())
    if info.name in ("FunctionalGraphicalModels", "FunctionalDataAnalysis", "KernelSmoothing",
        "CovarianceSelection", "ProximalBase", "ProximalSolvers"))
serialize(only(ARGS), (; julia_version=string(VERSION), dependencies, rows=comparison_fixtures()))
