@testset "Direct covariance and precision" begin
    # Exhaustive visit-count expectation for constant trajectories. This tests
    # both inclusion correction and omission of repeated diagonal visits.
    for m in (2,7,40),a in (.0001,.2,.9,1.)
        pi1,pi2 = FunctionalGraphicalModels._inclusion(m,a)
        expected1 = sum(binomial(big(m),r)*a^r*(1-a)^(m-r) for r in 1:m)
        expected2 = sum(binomial(big(m),r)*a^r*(1-a)^(m-r) for r in 2:m)
        @test pi1 ≈ Float64(expected1) rtol=1e-10
        @test pi2 ≈ Float64(expected2) rtol=1e-10
    end
    # Zero stochastic allowance and negligible positive bias recover an inverse.
    S = [1.0 .3 .0;.3 1.5 .2;.0 .2 .8]
    for normkind in (:linf,:l2), lp_solver in (:highs,:clarabel)
        result = bias_aclime(reshape(S,3,3,1);u=[0.],v=[0.],bias=[1e-7],bias_norm=normkind,lp_solver)
        @test result.columns ≈ inv(S) atol=2e-5
        @test maximum(result.violations)<1e-6
        @test all(result.budgets .>= vec(sum(abs.(result.columns);dims=1)).-1e-6)
    end
    # Adversarial star: interval LP controls column l1 and squared Frobenius.
    O = Matrix{Float64}(I,8,8); O[1,2:end].=.1; O[2:end,1].=.1
    W = copy(O); W[1,2:end].-=.04
    Z = interval_symmetrize(W)
    @test issymmetric(Z)
    @test opnorm(Z-O,1)<=2opnorm(W-O,1)+1e-6
    @test sum(abs2,Z-O)<=2sum(abs2,W-O)+1e-6
    @test all(min.(W,transpose(W)).-1e-7 .<= Z .<= max.(W,transpose(W)).+1e-7)
    # Explicit raw products on a full-domain window: no noise-diagonal use.
    fd = FDA.MultivariateFunctionalData([
        FDA.IrregularFunctionalData([[.1,.9],[.2,.8]],[[1.,3.],[2.,4.]];domain=(0.,1.)),
        FDA.IrregularFunctionalData([[.15,.85],[.25,.75]],[[2.,4.],[3.,5.]];domain=(0.,1.))])
    C = direct_covariance(fd,[.5];bandwidths=[.5])
    w = [2*.75*(1-((t-.5)/.5)^2) for t in (.1,.9)]
    w2 = [2*.75*(1-((t-.5)/.5)^2) for t in (.2,.8)]
    @test C.matrices[1,1,1,1] ≈ (prod(w)*3+prod(w2)*8)/2
    @test C.pair_subjects[1,1,1]==2
    @test_throws ArgumentError direct_covariance(fd,[0.];bandwidths=[.2])
    B = direct_covariance(fd,[0.];bandwidths=[.5],boundary=:local_linear)
    @test all(isfinite,B.matrices)
    @test_throws ArgumentError direct_covariance(fd,[.5];bandwidths=[.5],design=:synchronous)
    @test_throws ArgumentError bias_aclime(reshape(S,3,3,1);u=[0.],v=[0.],bias=[0.])
    # Scaled-lasso reduction on a diagonal positive Gram matrix and PSD receipt.
    D = DirectCovariance([.5],[.2],reshape(Matrix(Diagonal([1.,2.,3.])),3,3,1,1),
        200,20,fill(200,3,1,1),(;design=:independent))
    f = fit_direct(D;method=:scaled_lasso,bias_constant=1.)
    @test f.precisions[:,:,1] ≈ Matrix(Diagonal([1.,.5,1/3])) atol=1e-6
    @test f.diagnostics[1].covariance_psd_displacement==0
    @test !any(f.adjacency)
    # Multiscale l2 intersection remains feasible even for zero statistics.
    r = bias_aclime(zeros(3,3,2);u=[.1,.2],v=[.1,.2],bias=[.02,.005],bias_norm=:l2)
    @test all(isfinite,r.columns)
    @test maximum(r.violations)<1e-5
end

@testset "Stratified synchronous pair normalization" begin
    # Two-point Gauss quadrature in each of five strata integrates every
    # independent quadratic Epanechnikov factor exactly. Exhausting all
    # 2^10 choices gives an independent expectation check for the U-statistic.
    k,strata=2,5
    labels=repeat(1:strata;inner=k)
    acc=0.
    for mask in 0:1023
        times=[(s-1+.5+(((mask>>(r-1))&1)==1 ? 1 : -1)/(2sqrt(3)))/strata
            for (r,s) in enumerate(labels)]
        R=reshape([2*.75*(1-((t-.5)/.5)^2) for t in times],1,:)
        acc+=only(FunctionalGraphicalModels._stratified_pairs(R,labels,k,strata,1.))
    end
    @test acc/1024≈1 atol=1e-12
    times=sort!([(s-1+q)/5 for s in 1:5 for q in (.25,.75)])
    fd=FDA.MultivariateFunctionalData([
        FDA.IrregularFunctionalData([times,times],[ones(10),2ones(10)];domain=(0.,1.))])
    C=direct_covariance(fd,[.5];bandwidths=[.5],design=:stratified_synchronous)
    R=reshape([2*.75*(1-((t-.5)/.5)^2) for t in times],1,:)
    expected=2.5only(FunctionalGraphicalModels._stratified_pairs(R,labels,2,5,1.))
    @test C.matrices[1,1,1,1]≈expected
    @test C.settings.strata==5
    @test_throws ArgumentError direct_covariance(fd,[.5];bandwidths=[.5],design=:stratified_synchronous,strata=4)
end

include(joinpath(@__DIR__,"..","experiments","dfgm_factorized.jl"))
# The campaign helpers intentionally preserve FDA's earlier surface-diagonal,
# middle-half noise estimator. New FDA versions default to rotated/full noise.
# Compare the same explicit estimator on both sides, retaining old-FDA support.
function _dfgm_raw_reference_covariance(fd,meanfit;kwargs...)
    signature=Tuple{Type{FDA.SmoothedCovariance},typeof(fd),typeof(meanfit)}
    if hasmethod(FDA.fit,signature,(:noise_method,:noise_interval))
        FDA.fit(FDA.SmoothedCovariance,fd,meanfit;kwargs...,
            noise_method=:surface_diagonal,noise_interval=:middle_half)
    else
        FDA.fit(FDA.SmoothedCovariance,fd,meanfit;kwargs...)
    end
end

@testset "Exact irregular factorization against raw FDA fits" begin
    rng=MersenneTwister(710)
    times=[sort!(rand(rng,m)) for m in (4,6,5,7,4,6)]
    values=[sinpi.(t).+.3randn(rng,length(t)) for t in times]
    fd=FDA.IrregularFunctionalData(times,values;domain=(0.,1.))
    grid=collect(range(0,1;length=9))
    kernel=KernelSmoothing.GaussianKernel()
    meanfit=FDA.fit(FDA.SmoothedMean,fd;grid,bandwidth=.3,kernel,weighting=:subj)
    for bandwidth in (.3,FDA.CurveCV([.2,.4];folds=6))
        factored=_dfgm_factorized_covariance(fd,meanfit,grid;bandwidth)
        raw=_dfgm_raw_reference_covariance(fd,meanfit;grid,bandwidth,kernel,weighting=:subj,binning=nothing)
        @test factored.surface≈raw.surface atol=1e-10
        @test factored.noise_variance≈raw.noise_variance atol=1e-10
        @test factored.bandwidth==raw.bandwidth
        if bandwidth isa FDA.CurveCV
            @test factored.selection.errors≈raw.selection.errors atol=1e-10
        end
    end
end

@testset "DFGM design and projected target" begin
    settings = dfgm2020_settings()
    @test length(settings)==12
    @test all(s.n==200 && s.replicates==100 for s in settings)
    sim = simulate_dfgm2020(MersenneTwister(72);p=10,n=8,design=:irregular,
        schedule=:T6_9,integration_points=21,evaluation_grid=[.25,.5,.75])
    @test all(all(6 .<= FDA.nobs_per_curve(node) .<= 9) for node in sim.data.components)
    @test all(isposdef(Symmetric(sim.covariance[:,:,t])) for t in 1:3)
    @test all(sim.covariance[:,:,t]*sim.precision[:,:,t]≈I for t in 1:3)
    again = simulate_dfgm2020(MersenneTwister(72);p=10,n=8,design=:irregular,
        schedule=:T6_9,integration_points=21,evaluation_grid=[.25,.5,.75])
    @test sim.coefficient_covariance==again.coefficient_covariance
    @test FDA.component(sim.data,1).values==FDA.component(again.data,1).values
    @test length(sim.graph.events)==20sim.graph.change_count
end

include(joinpath(@__DIR__,"..","experiments","dfgm_common.jl"))
@testset "Exact common-grid aggregation against raw FDA fits" begin
    rng = MersenneTwister(9007)
    times = collect(range(0,1;length=6))
    values = [randn(rng).*sinpi.(times)+.3randn(rng,6) for _ in 1:8]
    fd = FDA.IrregularFunctionalData([copy(times) for _ in 1:8],values;domain=(0.,1.))
    grid = collect(range(0,1;length=9))
    kernel = KernelSmoothing.GaussianKernel()
    meanfit = FDA.fit(FDA.SmoothedMean,fd;grid,bandwidth=.3,kernel,weighting=:subj)
    for bandwidth in (.3,FDA.CurveCV([.2,.4];folds=8))
        compact = _dfgm_common_covariance(fd,meanfit,grid;bandwidth)
        raw = _dfgm_raw_reference_covariance(fd,meanfit;grid,bandwidth,kernel,weighting=:subj,binning=nothing)
        @test compact.surface ≈ raw.surface atol=1e-10
        @test compact.noise_variance ≈ raw.noise_variance atol=1e-10
        @test compact.bandwidth==raw.bandwidth
        if bandwidth isa FDA.CurveCV
            @test compact.selection.errors ≈ raw.selection.errors atol=1e-10
        end
    end
end
