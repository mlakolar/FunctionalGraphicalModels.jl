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
    @test fpca.kernel == :gaussian
    @test fpca.bandwidth == 0.3
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
    @test_throws ArgumentError fit_sparse_fpca(data;grid=grid,
        ncomponents=1,bandwidth=0.3,kernel=:unknown)
    @test_throws ArgumentError fit_sparse_fpca(data;grid=collect(range(0.1,0.9;length=9)),
        ncomponents=1,bandwidth=0.3)
    # A sufficiently wide compact kernel retains a nonsingular local design.
    compact = fit_sparse_fpca(data;grid=grid,ncomponents=1,
        bandwidth=1.1,noise_variance=0.04,kernel=:epanechnikov)
    @test compact.kernel == :epanechnikov
    @test compact.eigenvalues[1,1] ≈ 1 atol=1e-8
    @test compact.scores ≈ fpca.scores atol=1e-7
    compact_graph = fit_doubly(data;lambda=0.1,grid=grid,bandwidth=1.1,
        ncomponents=1,noise_variance=0.04,kernel=:epanechnikov)
    @test compact_graph.fpca.kernel == :epanechnikov
end

@testset "Sparse local-plane covariance weighted least-squares oracle" begin
    # Each subject contributes all ordered off-diagonal products. Different
    # observation counts make equal-subject weighting differ from equal-pair
    # weighting. The oracle uses a rectangular QR least-squares solve instead
    # of the implementation's accumulated 3-by-3 normal equations.
    schedules = [[0.0,0.2],[0.15,0.55,0.9],[0.05,0.35,0.7,1.0]]
    centered = [[2.0,-1.0],[0.3,1.7,-0.4],[-1.2,0.5,2.0,0.7]]
    firsttime,secondtime,products,pairweights = Float64[],Float64[],Float64[],Float64[]
    for i in eachindex(schedules)
        count = length(schedules[i])
        for a in 1:count, b in 1:count
            a == b && continue
            push!(firsttime,schedules[i][a]); push!(secondtime,schedules[i][b])
            push!(products,centered[i][a]*centered[i][b])
            push!(pairweights,1/(length(schedules)*count*(count-1)))
        end
    end
    @test length(products) == 2+6+12
    @test sum(pairweights) ≈ 1
    u,v,h = 0.45,0.6,1.1
    x,y = (firsttime.-u)./h,(secondtime.-v)./h
    design = hcat(ones(length(products)),x,y)
    for kernel in (:gaussian,:epanechnikov)
        # Normalized kernels include their constants and 1/h factors here.
        # Those constants cancel, as required by the local fit.
        kernelvalue(z) = kernel == :gaussian ? exp(-z*z/2)/(sqrt(2pi)*h) :
            0.75*max(1-z*z,0)/h
        kernelweights = kernelvalue.(x).*kernelvalue.(y)
        weights = pairweights.*kernelweights
        oracle = (sqrt.(weights).*design) \ (sqrt.(weights).*products)
        estimate = FunctionalGraphicalModels._dfgm_local_surface(firsttime,
            secondtime,products,pairweights,u,v,h,kernel)
        @test estimate ≈ oracle[1] atol=1e-11
        equal_pair_oracle = (sqrt.(kernelweights).*design) \
            (sqrt.(kernelweights).*products)
        @test abs(estimate-equal_pair_oracle[1]) > 1e-3
    end
    @test FunctionalGraphicalModels._dfgm_kernel(1.1,:epanechnikov) == 0
    @test FunctionalGraphicalModels._dfgm_kernel(1.1,:gaussian) > 0
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
    fit = fit_sparse_fpca(SparseFunctionalData(times,values);grid=grid,
        ncomponents=1,bandwidth=0.2,noise_variance=noise)
    # H independently represents linear interpolation on the evaluation grid.
    H = zeros(length(schedule),length(grid))
    for (a,t) in enumerate(schedule)
        k = searchsortedlast(grid,t)
        fraction = (t-grid[k])/(grid[k+1]-grid[k])
        H[a,k],H[a,k+1] = 1-fraction,fraction
    end
    covariance = fit.covariances[:,:,1]
    w = fit.weights
    @test count(>(1e-8),eigvals(Symmetric(sqrt.(w).*covariance.*transpose(sqrt.(w))))) >= 2
    full_observation_covariance = H*covariance*H' + noise*I
    retained_function = H*fit.eigenfunctions[:,1,1]
    retained_variance = fit.eigenvalues[1,1]
    score_observation_covariance = retained_variance.*retained_function
    truncated_observation_covariance = retained_variance.*
        (retained_function*retained_function') + noise*I
    wrong_scores = zeros(n)
    for i in 1:n
        centered = values[i,1]-H*fit.mean[1,:]
        expected = dot(score_observation_covariance,full_observation_covariance\centered)
        wrong_scores[i] = dot(score_observation_covariance,truncated_observation_covariance\centered)
        @test fit.scores[i,1,1] ≈ expected atol=1e-8
    end
    @test norm(fit.scores[:,1,1]-wrong_scores) > 1e-3
end
