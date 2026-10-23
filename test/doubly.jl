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
    fit = fit_doubly(DenseFunctionalData(values,grid);lambda=0.0,ncomponents=2)
    @test fit.method == :dense_reconstruction
    @test size(fit.adjacency) == (p,p,q)
    @test all(fit.positive_definite)
    for t in 1:q
        reconstructed = fit.reconstruction[:,:,t]
        @test fit.covariances[:,:,t] ≈ reconstructed'*reconstructed/n atol=1e-10
        @test fit.precisions[:,:,t] ≈ inv(fit.covariances[:,:,t]) atol=1e-7
        @test all(!fit.adjacency[j,j,t] for j in 1:p)
    end
    @test maximum(abs, sum(fit.reconstruction;dims=1)) < 1e-9
    @test_throws DimensionMismatch fit_doubly(DenseFunctionalData(values,grid);
        lambda=[0.1],ncomponents=2)
    @test_throws ArgumentError fit_doubly(DenseFunctionalData(values,grid);
        lambda=-0.1,ncomponents=2)
end

@testset "Sparse covariance smoothing and analytic conditional scores" begin
    # Constant random curves have covariance C(s,t)=1. Observation schedules
    # repeat every 12 subjects. The 60 sinusoidal scores have mean zero and
    # variance one within every schedule group, so covariance smoothing is exact.
    n = 60
    times = Matrix{Vector{Float64}}(undef,n,1)
    values = similar(times)
    scores = [sqrt(2)*sin(2pi*i/n) for i in 1:n]
    for i in 1:n
        times[i,1] = collect(range(0.02*mod(i,3),1-0.01*mod(i,4);
                                  length=5+mod(i,3)))
        values[i,1] = fill(scores[i],length(times[i,1]))
    end
    data = SparseFunctionalData(times,values)
    grid = collect(range(0,1; length=9))
    fpca = fit_sparse_fpca(data;grid=grid,ncomponents=1,
        bandwidth=0.3,noise_variance=0.04)
    @test fpca.eigenvalues[1,1] ≈ 1 atol=1e-8
    @test fpca.eigenfunctions[:,1,1] ≈ ones(9) atol=1e-8
    @test fpca.covariances[:,:,1] ≈ ones(9,9) atol=1e-8
    @test maximum(abs,fpca.mean) < 1e-10
    for i in 1:n
        count = length(times[i,1])
        @test fpca.scores[i,1,1] ≈ count/(count+0.04)*scores[i] atol=1e-7
    end
    fit = fit_doubly(data;lambda=0.1,grid=grid,bandwidth=0.3,
        ncomponents=1,noise_variance=0.04)
    @test fit.method == :sparse_conditional_scores
    @test all(isfinite,fit.precisions)
    @test !any(fit.adjacency)
    @test_throws ArgumentError fit_sparse_fpca(data;grid=grid,
        ncomponents=1,bandwidth=0)
    @test_throws ArgumentError fit_sparse_fpca(data;grid=grid,
        ncomponents=1,bandwidth=0.3,noise_variance=-1)
    @test_throws ArgumentError fit_sparse_fpca(data;grid=grid,
        ncomponents=2,bandwidth=0.3,noise_variance=0.04)
    duplicate_times = deepcopy(times)
    duplicate_times[1,1][2] = duplicate_times[1,1][1]
    @test_throws ArgumentError SparseFunctionalData(duplicate_times,values)
    @test_throws DimensionMismatch SparseFunctionalData(times,values[1:end-1,:])
end
