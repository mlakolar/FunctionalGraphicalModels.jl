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
    for normkind in (:linf,:l2)
        result = bias_aclime(reshape(S,3,3,1);u=[0.],v=[0.],bias=[1e-7],bias_norm=normkind)
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
