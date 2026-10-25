@testset "Dense doubly functional reconstruction and CLIME" begin
    n,p,q = 60,3,17
    grid = collect(range(0,1; length=q))
    values = zeros(n,p,q)
    for i in 1:n, t in 1:q
        phi = sqrt(2)*cos(2pi*grid[t])
        a,b,c,d,e,f = [sqrt(2)*sin(2pi*l*i/n) for l in 1:6]
        values[i,1,t] = 5 + a+0.25*b*phi
        values[i,2,t] = -2 + 0.4*a+c+0.25*d*phi
        values[i,3,t] = 1 + e+0.25*f*phi
    end
    data = FunctionalGraphicalModels._dense_data(values,grid)
    fit = fit_doubly(data;lambda=0.0,ncomponents=2)
    @test fit.method == :dense_reconstruction
    @test fit.lambda == zeros(q)
    @test fit.threshold == 1e-6
    @test size(fit.adjacency) == (p,p,q)
    @test all(fit.positive_definite)
    for t in 1:q
        reconstructed = fit.reconstruction[:,:,t]
        @test fit.covariances[:,:,t] ≈ reconstructed'*reconstructed/n atol=1e-10
        @test fit.precisions[:,:,t] ≈ inv(fit.covariances[:,:,t]) atol=1e-7
        @test all(!fit.adjacency[j,j,t] for j in 1:p)
    end
    @test maximum(abs, sum(fit.reconstruction;dims=1)) < 1e-9
    @test_throws DimensionMismatch fit_doubly(FunctionalGraphicalModels._dense_data(values,grid);
        lambda=[0.1],ncomponents=2)
    @test_throws ArgumentError fit_doubly(FunctionalGraphicalModels._dense_data(values,grid);
        lambda=-0.1,ncomponents=2)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,ncomponents=2,bandwidth=0.2)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,ncomponents=2,noise_variance=0.04)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,ncomponents=2,grid=grid)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,ncomponents=2,kernel=:epanechnikov)
    mixed = FDA.MultivariateFunctionalData([
        FDA.component(data, 1),
        FDA.IrregularFunctionalData([copy(grid) for _ in 1:n],
            [values[i,2,:] for i in 1:n])])
    @test_throws ArgumentError fit_doubly(mixed;lambda=0.1,ncomponents=2)
    incomplete_values = Matrix{Union{Missing,Float64}}(values[:,1,:])
    incomplete_values[1,1] = missing
    incomplete = FDA.MultivariateFunctionalData([
        FDA.GridFunctionalData(incomplete_values, grid)])
    @test_throws ArgumentError fit_doubly(incomplete;lambda=0.1,ncomponents=2)
end

@testset "Sparse FDA covariance smoothing and analytic conditional scores" begin
    # Constant random curves have covariance C(s,t)=1. The schedules repeat
    # every 12 subjects, and every schedule group has variance-one scores.
    n = 60
    times = Matrix{Vector{Float64}}(undef,n,1)
    values = similar(times)
    scores = [sqrt(2)*sin(2pi*i/n) for i in 1:n]
    for i in 1:n
        times[i,1] = collect(range(0.02*mod(i,3),1-0.01*mod(i,4);
                                  length=5+mod(i,3)))
        values[i,1] = fill(scores[i],length(times[i,1]))
    end
    data = FDA.MultivariateFunctionalData([
        FDA.IrregularFunctionalData(times[:,1], values[:,1])])
    grid = collect(range(0,1; length=9))
    fit = fit_doubly(data;lambda=0.1,grid=grid,bandwidth=0.3,
        ncomponents=1,noise_variance=0.04)
    fpca = fit.fpca
    node = FDA.component(fpca, 1)
    @test fpca isa FDA.MultivariateFPCA
    @test node.settings.kernel isa FunctionalGraphicalModels.KernelSmoothing.GaussianKernel
    @test node.settings.mean_bandwidth == 0.3
    @test node.settings.covariance_bandwidth == 0.3
    @test node.settings.diagonal_bandwidth == 0.3
    @test node.settings.weighting == :subj
    @test node.settings.binning === nothing
    @test node.settings.max_k == length(grid)
    @test node.settings.eigen_tol == 100eps(Float64)
    @test node.eigenvalues[1] ≈ 1 atol=1e-8
    @test node.basis.values[:,1] ≈ ones(9) atol=1e-8
    @test node.settings.full_projected_covariance ≈ ones(9,9) atol=1e-8
    @test maximum(abs,node.mean.values) < 1e-10
    for i in 1:n
        count = length(times[i,1])
        @test node.scores[i,1] ≈ count/(count+0.04)*scores[i] atol=1e-7
    end
    @test fit.method == :sparse_conditional_scores
    @test fit.lambda == fill(0.1,length(grid))
    @test fit.threshold == 1e-6
    @test all(isfinite,fit.precisions)
    @test !any(fit.adjacency)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,grid=grid,
        ncomponents=1,bandwidth=0)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,grid=grid,
        ncomponents=1,bandwidth=0.3,noise_variance=-1)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,grid=grid,
        ncomponents=2,bandwidth=0.3,noise_variance=0.04)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,grid=grid,
        ncomponents=1,bandwidth=0.3,kernel=:unknown)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,
        grid=collect(range(0.1,0.9;length=9)),ncomponents=1,bandwidth=0.3)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,bandwidth=0.3)
    @test_throws ArgumentError fit_doubly(data;lambda=0.1,grid=grid)
    # A sufficiently wide compact kernel retains a nonsingular local design.
    compact_graph = fit_doubly(data;lambda=0.1,grid=grid,bandwidth=1.1,
        ncomponents=1,noise_variance=0.04,kernel=:epanechnikov)
    compact = FDA.component(compact_graph.fpca, 1)
    @test compact.settings.kernel isa FunctionalGraphicalModels.KernelSmoothing.EpanechnikovKernel
    @test compact.eigenvalues[1] ≈ 1 atol=1e-8
    @test compact.scores ≈ node.scores atol=1e-7

    # Each noise value must reach the corresponding FDA node without changing
    # subject order. Independent errors shrink its conditional scores more.
    two_nodes = FDA.MultivariateFunctionalData([
        FDA.component(data, 1), FDA.component(data, 1)])
    multinoise = fit_doubly(two_nodes;lambda=0.1,grid=grid,bandwidth=0.3,
        ncomponents=1,noise_variance=[0.04,0.08])
    @test FDA.component(multinoise.fpca, 1).noise_variance == 0.04
    @test FDA.component(multinoise.fpca, 2).noise_variance == 0.08
    for i in 1:n
        count = length(times[i,1])
        @test FDA.scores(multinoise.fpca, 2)[i,1] ≈ count/(count+0.08)*scores[i] atol=1e-7
    end
    for invalid_noise in ([-0.1,0.04], [NaN,0.04], [nothing,0.04], [:bad,0.04], "bad")
        @test_throws ArgumentError fit_doubly(two_nodes;lambda=0.1,
            grid=grid,bandwidth=0.3,ncomponents=1,noise_variance=invalid_noise)
    end
    @test_throws DimensionMismatch fit_doubly(two_nodes;lambda=0.1,
        grid=grid,bandwidth=0.3,ncomponents=1,noise_variance=[0.04])
end

@testset "Sparse BLUP uses covariance modes beyond retained scores" begin
    n = 40
    schedule = [0.02,0.06,0.14,0.27,0.44,0.73,0.98]
    times = Matrix{Vector{Float64}}(undef,n,1)
    values = similar(times)
    for i in 1:n
        a,b = sqrt(2)*sin(2pi*i/n),0.7*sqrt(2)*cos(2pi*i/n)
        times[i,1] = copy(schedule)
        values[i,1] = a .+ b.*(2schedule.-1)
    end
    grid = collect(range(0,1;length=17))
    noise = 0.05
    data = FDA.MultivariateFunctionalData([
        FDA.IrregularFunctionalData(times[:,1], values[:,1])])
    graph = fit_doubly(data;lambda=0.1,grid=grid,ncomponents=1,
        bandwidth=0.2,noise_variance=noise)
    fit = FDA.component(graph.fpca, 1)
    # H independently represents linear interpolation on the evaluation grid.
    H = zeros(length(schedule),length(grid))
    for (a,t) in enumerate(schedule)
        k = searchsortedlast(grid,t)
        fraction = (t-grid[k])/(grid[k+1]-grid[k])
        H[a,k],H[a,k+1] = 1-fraction,fraction
    end
    covariance = fit.settings.full_projected_covariance
    w = fit.weights
    @test count(>(1e-8),eigvals(Symmetric(sqrt.(w).*covariance.*transpose(sqrt.(w))))) >= 2
    full_observation_covariance = H*covariance*H' + noise*I
    retained_function = H*fit.basis.values[:,1]
    retained_variance = fit.eigenvalues[1]
    score_observation_covariance = retained_variance.*retained_function
    truncated_observation_covariance = retained_variance.*
        (retained_function*retained_function') + noise*I
    wrong_scores = zeros(n)
    for i in 1:n
        centered = values[i,1]-H*vec(fit.mean.values)
        expected = dot(score_observation_covariance,full_observation_covariance\centered)
        wrong_scores[i] = dot(score_observation_covariance,truncated_observation_covariance\centered)
        @test fit.scores[i,1] ≈ expected atol=1e-8
    end
    @test norm(fit.scores[:,1]-wrong_scores) > 1e-3

    # Estimated noise follows FDA's default estimator: the rotated diagonal
    # (Yao, Mueller and Wang 2005, which the DFGM paper cites) integrated over
    # the full interval, chosen by the FDA issue #17 benchmark (2026-10-25).
    # Earlier values on this fixture: PSD-diagonal average 0.09556629502099738,
    # surface diagonal over the middle half 0.11159325235704108. The curves are
    # noise-free, so the new estimate (0.031) is also the closest to zero.
    estimated_graph = fit_doubly(data; lambda=0.1, grid=grid,
        ncomponents=1, bandwidth=0.2)
    estimated = FDA.component(estimated_graph.fpca, 1)
    @test estimated.settings.noise_variance_source == :estimated
    raw = estimated.settings.covariance_input
    @test raw.noise_details.noise_method == :rotated
    @test raw.noise_details.noise_interval == :full
    @test raw.noise_variance ≈ 0.030864993605722504 rtol=1e-10
    # FDA now floors estimated PACE noise at 0.025 times its maximum fitted
    # signal variance. Check both the unchanged raw estimate and the value
    # actually used by the conditional solve; supplied noise remains unchanged.
    expected_floor = 0.025 * maximum(diag(estimated.settings.full_projected_covariance))
    @test estimated.settings.noise_floor ≈ expected_floor rtol=1e-12
    @test estimated.noise_variance ≈ max(raw.noise_variance,expected_floor) rtol=1e-12
    @test estimated.settings.noise_floor_applied == (raw.noise_variance < expected_floor)
    # The adopted estimate must also reach the conditional-score solve.
    observation_covariance = H*estimated.settings.full_projected_covariance*H' +
        estimated.noise_variance*I
    crosscovariance = estimated.eigenvalues[1] .* (H*estimated.basis.values[:,1])
    for i in 1:n
        centered = values[i,1] - H*vec(estimated.mean.values)
        expected = dot(crosscovariance, observation_covariance\centered)
        @test estimated.scores[i,1] ≈ expected rtol=1e-10 atol=1e-12
    end
end
