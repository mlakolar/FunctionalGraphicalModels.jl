using LinearAlgebra
using Random

@testset "Functional neighborhood selection" begin
    rng = MersenneTwister(241024)
    n,m = 40,2
    raw = randn(rng,n,2m)
    raw .-= sum(raw; dims=1)/n
    X = sqrt(n)*Matrix(qr(raw).Q)[:,1:2m]
    truth = [1.0 0.4; -0.2 0.6; 0.1 0.2; 0.0 -0.1]
    Y = X*truth
    scores = zeros(n,3,m)
    scores[:,1,:] = Y
    scores[:,2,:] = X[:,1:m]
    scores[:,3,:] = X[:,m+1:2m]
    lambda = 0.35
    result = functional_neighborhood(scores,lambda; atol=1e-9,rtol=1e-8)
    expected = copy(truth)
    for block in (1:m,m+1:2m)
        expected[block,:] .*= max(1-lambda/norm(truth[block,:]),0)
    end
    @test result.coefficients[1] ≈ expected atol=2e-7
    @test all(d -> d.converged,result.diagnostics)
    @test maximum(d.kkt_residual for d in result.diagnostics) < 2e-7
    @test result.diagnostics[1].objective ≈
        sum(abs2,Y-X*expected)/(2n) + lambda*(norm(expected[1:m,:])+norm(expected[m+1:2m,:])) atol=1e-9
    @test result.directed[1,:] == [false,true,false]
    @test issymmetric(result.adjacency)
    @test all(!result.adjacency[j,j] for j in 1:3)

    lambda_max = neighborhood_lambda_max(scores,1)
    @test lambda_max ≈ norm(truth[1:m,:])
    zero_fit = functional_neighborhood(scores,fill(2lambda_max,3))
    @test iszero(zero_fit.coefficients[1])
    @test zero_fit.diagnostics[1].iterations == 0

    shifted = copy(scores)
    shifted .+= reshape([5.0,-2.0,1.0],1,3,1)
    shifted_fit = functional_neighborhood(shifted,lambda; atol=1e-9,rtol=1e-8)
    @test shifted_fit.coefficients[1] ≈ result.coefficients[1] atol=1e-8

    and_fit = functional_neighborhood(scores,lambda; rule=:and)
    or_fit = functional_neighborhood(scores,lambda; rule=:or)
    @test and_fit.adjacency == (and_fit.directed .& transpose(and_fit.directed))
    @test or_fit.adjacency == (or_fit.directed .| transpose(or_fit.directed))
    limited = functional_neighborhood(scores,lambda; maxiter=1,atol=1e-14,rtol=0)
    @test !limited.diagnostics[1].converged
    @test limited.diagnostics[1].iterations == 1
    @test_throws ArgumentError functional_neighborhood(scores,-0.1)
    @test_throws ArgumentError functional_neighborhood(scores,NaN)
    @test_throws ArgumentError functional_neighborhood(scores,lambda; threshold=-0.1)
    @test_throws ArgumentError functional_neighborhood(scores,lambda; rule=:xor)
    @test_throws ArgumentError functional_neighborhood(scores,lambda; rho=0)
    @test_throws DimensionMismatch functional_neighborhood(scores,[0.1,0.2])
    @test_throws ArgumentError functional_neighborhood(fill(NaN,2,2,1),0.1)
end

@testset "Target-specific FPCA projection for EJS" begin
    n,p,t = 30,3,41
    grid = collect(range(0,1; length=t))
    rng = MersenneTwister(301)
    values = zeros(n,p,t)
    # Distinct target covariances require distinct bases. Regressing each
    # node's marginal FPCA scores would fail this explicit-projection identity.
    a = randn(rng,n,p)
    b = randn(rng,n,p)
    for i in 1:n,k in 1:p,s in 1:t
        values[i,k,s] = 2a[i,k]*sinpi(k*grid[s]) + 0.3b[i,k]*cospi(grid[s])
    end
    data = DenseFunctionalData(values,grid)
    fpca = fit_fpca(data; ncomponents=2,basis=:marginal)
    fit = functional_neighborhood(data,0.08; ncomponents=2)
    for j in 1:p
        projected = zeros(n,p,2)
        for i in 1:n,k in 1:p,l in 1:2
            projected[i,k,l] = sum((values[i,k,s]-fpca.mean[k,s])*
                fpca.weights[s]*fpca.eigenfunctions[s,l,j] for s in 1:t)
        end
        explicit = functional_neighborhood(projected,0.08)
        @test fit.coefficients[j] ≈ explicit.coefficients[j] atol=1e-9
    end
    @test size(fit.target_bases) == (t,2,p)
end

@testset "Selective cross-validation (EJS Algorithm 2)" begin
    rng = MersenneTwister(42)
    scores = randn(rng,30,3,1)
    scores[:,1,1] = 1.5scores[:,2,1] + 0.05randn(rng,30)
    selected = select_neighborhood_scv(scores,[0.01,10.0],[0.0,100.0]; folds=3)
    @test selected.fit.directed[1,2]
    @test !selected.fit.directed[1,3]
    @test size(selected.criteria) == (2,2,3)
    @test all(selected.candidate_converged)
    @test all(d -> d.converged,selected.fit.diagnostics)
    @test selected.fold_ids == mod1.(1:30,3)
    # Full thresholding removes every source. SCV then reduces to the mean
    # held-out RSS from training-fold intercepts, with zero complexity penalty.
    mu = sum(scores; dims=1)/30
    Y = reshape(scores[:,1,1].-mu[1,1,1],:,1)
    null_rss = sum(sum(abs2,Y[findall(==(f),selected.fold_ids),:].-
        sum(Y[findall(!=(f),selected.fold_ids),:];dims=1)/20) for f in 1:3)/3
    @test selected.criteria[1,2,1] ≈ null_rss
    @test_throws ArgumentError select_neighborhood_scv(scores,Float64[],[0.0])
    @test_throws ArgumentError select_neighborhood_scv(scores,[0.1],[0.0]; folds=1)
    @test_throws DimensionMismatch select_neighborhood_scv(scores,[0.1],[0.0]; fold_ids=[1,2])
end
