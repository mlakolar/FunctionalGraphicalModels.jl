@testset "Functional graphical lasso objective and optimality" begin
    # With one node, the entire diagonal block is unpenalized.
    S = [2.0 0.4; 0.4 1.0]
    fit = functional_glasso(S, [2], 100.0)
    @test fit.diagnostics.converged
    @test fit.precision ≈ inv(S) atol=3e-6
    @test !any(fit.adjacency)
    # For scalar blocks the exact dual solution clips the covariance.
    lambda = 0.1
    exact = inv([2.0 0.3;0.3 1.0])
    fit = functional_glasso(S, 1, lambda)
    @test fit.diagnostics.converged
    @test fit.precision ≈ exact atol=3e-6
    @test fit.diagnostics.kkt_residual < 1e-5
    @test fit.adjacency == Bool[0 1;1 0]
    emptyfit = functional_glasso(S, 1, 0.4)
    @test !any(emptyfit.adjacency)
    @test diag(emptyfit.precision) ≈ [0.5,1.0] atol=3e-6
    @test isposdef(Symmetric(fit.precision))
    @test norm(fit.precision-fit.sparse_precision) < 1e-5
    # A non-scalar off-diagonal block is thresholded by Frobenius norm.
    C = [2.0 0 0.3 0;0 1.5 0 0.2;0.3 0 1.2 0;0 0.2 0 1.0]
    f = functional_glasso(C, [2,2], hypot(0.3,0.2)+0.01)
    @test f.diagnostics.converged
    @test !any(f.adjacency)
    @test f.precision[1:2,1:2] ≈ inv(C[1:2,1:2]) atol=4e-6
    @test functional_glasso(C,[2,2],0).precision ≈ inv(C) atol=4e-6
    @test_throws ArgumentError functional_glasso(S,[1],0.1)
    @test_throws ArgumentError functional_glasso(S,1,-1)
    @test_throws ArgumentError functional_glasso([1.0 1+1e-11;1+1e-11 1],1,1e-12)
    @test_throws ArgumentError functional_glasso(1e200 .* S,1,0.1)
    @test_throws ArgumentError functional_glasso([1.0 2;2 1],1,0.1)
    @test_throws ArgumentError functional_glasso(ones(2,2),[2],0.1)
    @test_throws ArgumentError functional_glasso(S,1,0.1;threshold=-1)
    @test_logs (:warn,r"did not converge") begin
        failed = functional_glasso(S,1,lambda;maxiter=1)
        @test !failed.diagnostics.converged
    end
    # Node-wise orthogonal rotations preserve block graph and transform precision.
    Q = [0.0 1;-1 0]
    R = zeros(4,4); R[1:2,1:2]=Q; R[3:4,3:4]=Matrix{Float64}(I,2,2)
    f1 = functional_glasso(C,[2,2],0.07)
    f2 = functional_glasso(R*C*R',[2,2],0.07)
    @test f1.adjacency == f2.adjacency
    @test f2.precision ≈ R*f1.precision*R' atol=1e-5
end

@testset "FDA dense graphical lasso integration" begin
    sim = simulate_qiao2019(MersenneTwister(46); n=40, p=3, ncomponents=2,
                            ntime=15, noise_sd=0)
    fit = functional_glasso(sim.data, 0.1; ncomponents=2)
    scores = FDA.score_array(fit.fpca)
    direct = functional_glasso(scores, 0.1)
    @test fit.precision ≈ direct.precision atol=1e-10
    @test fit.covariance ≈ score_covariance(scores) atol=1e-12
    @test fit.adjacency == direct.adjacency
    @test fit.settings.lambda == 0.1
    @test fit.settings.threshold == 1e-6
    @test all(s -> s.maxK == 2, fit.fpca.settings.node_settings)
    @test occursin("settings=", sprint(show, fit))
end
