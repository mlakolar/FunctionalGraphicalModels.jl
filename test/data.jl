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

@testset "B-spline smoothing" begin
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
    # One constant basis function fits each curve by its weighted mean, at the
    # last grid point too.
    w = quadrature_weights(grid)
    level = [sum(w .* d.values[i, j, :]) / sum(w) for i in 1:3, j in 1:2]
    @test smooth_bspline(d; nbasis=1, degree=0).values ≈ repeat(level, 1, 1, 21)
    @test_throws ArgumentError smooth_bspline(d; nbasis=2)
    @test_throws ArgumentError smooth_bspline(d; nbasis=5, penalty=-1)
    @test_throws ArgumentError smooth_bspline(d; nbasis=22)

    # Knots are equally spaced on the grid interval and extended by `degree`
    # knots on each side. With knots on the grid, hat functions give the
    # identity, cubic B-splines take the values 1/6, 2/3, 1/6, and piecewise
    # constants put the last grid point in the last interval.
    design = FunctionalGraphicalModels._bspline_design
    @test design(grid, 21, 1) ≈ Matrix(I, 21, 21) atol=1e-14
    @test design(grid, 23, 3)[11, 11:13] ≈ [1, 4, 1] ./ 6
    @test design(grid, 20, 0) == [Matrix(1.0I, 20, 20); [zeros(1, 19) 1.0]]
    for degree in 0:3
        B = design(grid, 7, degree)
        @test size(B) == (21, 7)
        @test all(>=(0), B)
        @test vec(sum(B; dims=2)) ≈ ones(21)
    end
end
