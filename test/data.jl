@testset "Dense functional data and integration" begin
    grid = [0.0, 0.1, 0.4, 1.0]
    w = quadrature_weights(grid)
    @test w ≈ [0.05, 0.2, 0.45, 0.3]
    @test sum(w) ≈ 1
    @test dot(w, grid) ≈ 0.5
    @test_throws ArgumentError quadrature_weights([0.0])
    @test_throws ArgumentError quadrature_weights([0.0, 0.0, 1.0])
    @test_throws ArgumentError quadrature_weights([1.0, 0.0])
    @test_throws ArgumentError quadrature_weights([0.0, Inf])
    values = reshape(collect(1.0:24.0), 3, 2, 4)
    d = DenseFunctionalData(values, grid)
    values[1] = -10
    grid[1] = -1
    @test d.values[1] == 1
    @test d.grid[1] == 0
    @test_throws DimensionMismatch DenseFunctionalData(ones(3, 2, 3), d.grid)
    @test_throws ArgumentError DenseFunctionalData(fill(NaN, 3, 2, 4), d.grid)
    @test_throws ArgumentError DenseFunctionalData(zeros(0, 2, 4), d.grid)
end

@testset "Node-major score covariance" begin
    X = [-1.0 -2.0 2.0 1.0; -1.0 2.0 -2.0 1.0;
          1.0 -2.0 -2.0 -1.0; 1.0 2.0 2.0 -1.0]
    scores = permutedims(reshape(X, 4, 2, 2), (1, 3, 2))
    @test score_covariance(scores) ≈ transpose(X) * X / 4
    @test score_covariance(scores; corrected=true) ≈ transpose(X) * X / 3
    shifted = copy(scores)
    shifted[:, 1, 1] .+= 3
    @test score_covariance(shifted) ≈ score_covariance(scores)
    @test score_covariance(shifted; center=false)[1, 1] ≈ 10
    @test_throws ArgumentError score_covariance(zeros(1, 2, 2); corrected=true)
    @test_throws ArgumentError score_covariance(fill(NaN, 3, 2, 2))
end

@testset "BSpline integration" begin
    grid = collect(range(0.0, 1.0; length=21))
    values = zeros(3, 2, length(grid))
    for j in 1:2, i in 1:3
        values[i, j, :] .= i .+ j .* grid
    end
    d = DenseFunctionalData(values, grid)
    smoothed = smooth_bspline(d; nbasis=2, degree=1)
    @test smoothed.values ≈ d.values atol=1e-12
    @test smoothed.grid == d.grid
    @test all(isfinite, smooth_bspline(d; nbasis=5, penalty=0.01).values)
    @test_throws ArgumentError smooth_bspline(d; nbasis=2)
    @test_throws ArgumentError smooth_bspline(d; nbasis=5, penalty=-1)
    @test_throws ArgumentError smooth_bspline(d; nbasis=22)
end
