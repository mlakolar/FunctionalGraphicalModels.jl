@testset "Joint graphical lasso: analytic solutions and diagnostics" begin
    S = reshape([1.0, 0.45, 0.45, 1.0], 2, 2, 1)
    # With one component, the elementwise and group penalties sum to the
    # ordinary graphical-lasso penalty, independently of alpha.
    expected = inv([1.0 0.35; 0.35 1.0])
    for alpha in (0.0, 0.35, 1.0)
        fit = joint_graphical_lasso(S; lambda=0.1, alpha=alpha,
            abstol=1e-9, reltol=1e-8, maxiter=5000)
        @test fit.converged
        @test all(fit.positive_definite)
        @test fit.precision[:,:,1] ≈ expected atol=2e-6
        @test norm(fit.precision-fit.likelihood_precision) == fit.primal_residual
        @test isfinite(fit.objective)
    end
    S2 = cat(S, reshape([1.0, -0.25, -0.25, 1.0], 2, 2, 1); dims=3)
    zero_penalty = joint_graphical_lasso(S2; lambda=0.0,
        abstol=1e-9, reltol=1e-8, maxiter=5000)
    @test zero_penalty.converged
    for l in 1:2
        @test zero_penalty.precision[:,:,l] ≈ inv(S2[:,:,l]) atol=2e-6
    end
    diagonal = joint_graphical_lasso(S2; lambda=2.0, alpha=0.0)
    @test diagonal.converged
    @test diagonal.precision[1,2,:] == zeros(2)
    @test diagonal.precision[1,1,:] ≈ ones(2) atol=1e-4
    stopped = joint_graphical_lasso(S; lambda=0.1, maxiter=1,
        abstol=1e-14, reltol=1e-14)
    @test !stopped.converged
    @test stopped.iterations == 1
    @test_throws ArgumentError joint_graphical_lasso(S; lambda=-0.1)
    @test_throws ArgumentError joint_graphical_lasso(S; lambda=0.1, alpha=1.1)
    @test_throws ArgumentError joint_graphical_lasso(S; lambda=0.1, rho=0)
    @test_throws ArgumentError joint_graphical_lasso(S; lambda=0.1, maxiter=0)
    @test_throws DimensionMismatch joint_graphical_lasso(zeros(2,3,1); lambda=0.1)
    @test_throws ArgumentError joint_graphical_lasso(reshape(ones(2,2),2,2,1); lambda=0.0)
    almost_psd = reshape([1.0,1+1e-12,1+1e-12,1.0],2,2,1)
    @test_throws ArgumentError joint_graphical_lasso(almost_psd; lambda=1e-14)
    @test joint_graphical_lasso(almost_psd; lambda=0.1, maxiter=1) isa JointGraphicalLassoFit
    @test_throws ArgumentError joint_graphical_lasso(S; lambda=0.1, rho=nextfloat(0.0))
end

@testset "Pooled rank can exceed each node's rank" begin
    # Two rank-(n-1) node covariances can pool to rank p*(n-1).
    # FDA PR #16 restores this range; no local pooled eigenstep is needed.
    n, p, q = 3, 2, 9
    grid = collect(range(0, 1; length=q))
    values = randn(MersenneTwister(416), n, p, q)
    data = FunctionalGraphicalModels._dense_data(values, grid)
    m = p * (n - 1)
    fit = fit_separable(data; lambda=2.0, ncomponents=m)
    @test fit.diagnostics.converged
    @test fit.fpca isa FDA.MultivariateFPCA
    @test fit.fpca.basis == :pooled
    @test fit.fpca.settings.max_k == m
    @test size(FDA.score_array(fit.fpca)) == (n, p, m)
    @test length(fit.fpca.pooled.eigenvalues) == m
    @test all(>(0), fit.fpca.pooled.eigenvalues)
    @test !any(fit.adjacency)
    @test FDA.to_array(FDA.reconstruct(fit.fpca)) ≈ values atol=1e-12
    @test_throws ArgumentError fit_separable(data; lambda=2.0, ncomponents=m+1)
end

@testset "Cross-component penalty normalization" begin
    # alpha=1 decouples the component likelihoods and penalties exactly.
    correlations = [0.55,-0.35,0.05]
    S = zeros(2,2,length(correlations))
    for l in eachindex(correlations)
        S[:,:,l] .= [1.0 correlations[l]; correlations[l] 1.0]
    end
    independent = joint_graphical_lasso(S;lambda=0.15,alpha=1.0,
        abstol=1e-9,reltol=1e-8,maxiter=5000)
    @test independent.converged
    for l in eachindex(correlations)
        r = sign(correlations[l])*max(abs(correlations[l])-0.15,0)
        @test independent.precision[:,:,l] ≈ inv([1.0 r; r 1.0]) atol=2e-6
    end
    # For M identical components at alpha=0, symmetry gives M equal
    # precision estimates. Dividing the objective by M leaves a scalar
    # graphical lasso with effective penalty lambda/sqrt(M).
    m,lambda,r = 3,0.15,0.55
    duplicated = repeat(reshape([1.0,r,r,1.0],2,2,1);outer=(1,1,m))
    grouped = joint_graphical_lasso(duplicated;lambda=lambda,alpha=0.0,
        abstol=1e-9,reltol=1e-8,maxiter=5000)
    expected_r = r-lambda/sqrt(m)
    expected = inv([1.0 expected_r; expected_r 1.0])
    @test grouped.converged
    for l in 1:m
        @test grouped.precision[:,:,l] ≈ expected atol=2e-6
    end
end

@testset "Partially separable dense data" begin
    n, p, q = 60, 3, 21
    grid = collect(range(0,1; length=q))
    values = zeros(n,p,q)
    for i in 1:n, t in 1:q
        a = sqrt(2)*sin(2pi*i/n)
        b = sqrt(2)*sin(4pi*i/n)
        c = sqrt(2)*sin(6pi*i/n)
        d = sqrt(2)*sin(8pi*i/n)
        e = sqrt(2)*sin(10pi*i/n)
        f = sqrt(2)*sin(12pi*i/n)
        phi = sqrt(2)*cos(2pi*grid[t])
        values[i,1,t] = a + 0.3*b*phi
        values[i,2,t] = 0.5*a+c + 0.3*d*phi
        values[i,3,t] = e + 0.3*f*phi
    end
    fit = fit_separable(FunctionalGraphicalModels._dense_data(values,grid);
        lambda=0.1, ncomponents=2, alpha=0.5, abstol=1e-8, reltol=1e-7)
    @test fit.diagnostics.converged
    @test fit.settings.threshold == 1e-6
    @test fit.diagnostics.settings.alpha == 0.5
    @test fit.fpca.settings.max_k == 2
    @test occursin("settings=", sprint(show, fit))
    @test size(fit.precisions) == (p,p,2)
    @test all(fit.adjacency .== dropdims(any(fit.component_adjacency; dims=3); dims=3))
    @test fit.adjacency[1,2]
    @test !fit.adjacency[1,3]
    @test !fit.adjacency[2,3]
    @test all(!fit.adjacency[j,j] for j in 1:p)
    for l in 1:2
        @test diag(fit.correlations[:,:,l]) ≈ ones(p)
        @test fit.precisions[:,:,l] ≈ fit.inverse_correlations[:,:,l] ./
            (fit.scales[:,l]*fit.scales[:,l]')
    end
end
