import FunctionalDataAnalysis as FDA

@testset "FDA functional data and integration" begin
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
    d = FDA.MultivariateFunctionalData([FDA.GridFunctionalData(values[:, j, :], grid) for j in 1:2])
    values[1] = -10
    grid[1] = -1
    @test FDA.to_array(d)[1] == 1
    @test FDA.component(d, 1).grid[1] == 0
    @test_throws DimensionMismatch FunctionalGraphicalModels._dense_data(ones(3, 2, 3), FDA.component(d, 1).grid)
    @test_throws ArgumentError FunctionalGraphicalModels._dense_data(fill(NaN, 3, 2, 4), FDA.component(d, 1).grid)
    @test_throws ArgumentError FunctionalGraphicalModels._dense_data(zeros(0, 2, 4), FDA.component(d, 1).grid)
    @test_throws ArgumentError FunctionalGraphicalModels._dense_data(zeros(3, 0, 4), FDA.component(d, 1).grid)
    copied_values, copied_grid = FunctionalGraphicalModels._dense_parts(d)
    copied_values[1] = -3
    copied_grid[1] = -3
    @test FDA.to_array(d)[1] == 1
    @test FDA.component(d, 1).grid[1] == 0
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
    d = FDA.MultivariateFunctionalData([FDA.GridFunctionalData(values[:, j, :], grid) for j in 1:2])
    smoothed = smooth_bspline(d; nbasis=2, degree=1)
    @test FDA.to_array(smoothed) ≈ FDA.to_array(d) atol=1e-12
    @test FDA.component(smoothed, 1).grid == FDA.component(d, 1).grid
    @test all(isfinite, FDA.to_array(smooth_bspline(d; nbasis=5, penalty=0.01)))
    @test smoothed isa FDA.MultivariateFunctionalData
    @test FDA.component(smoothed, 1).metadata.method === :pspline
    @test FDA.component(smoothed, 1).metadata.nbasis == 2
    @test FDA.component(smoothed, 1).metadata.degree == 1
    @test FDA.component(smoothed, 1).metadata.penalty == 0
    labelled = FunctionalGraphicalModels._dense_data(values, grid; ids=[:a, :b, :c])
    @test FDA.ids(smooth_bspline(labelled; nbasis=2, degree=1)) == [:a, :b, :c]
    # One constant basis function fits each curve by its weighted mean, at the
    # last grid point too.
    w = quadrature_weights(grid)
    level = [sum(w .* FDA.to_array(d)[i, j, :]) / sum(w) for i in 1:3, j in 1:2]
    @test FDA.to_array(smooth_bspline(d; nbasis=1, degree=0)) ≈ repeat(level, 1, 1, 21)
    @test_throws ArgumentError smooth_bspline(d; nbasis=2)
    @test_throws ArgumentError smooth_bspline(d; nbasis=5, penalty=-1)
    @test_throws ArgumentError smooth_bspline(d; nbasis=22)

    # Dense-only FGM smoothing rejects irregular, missing, and mismatched grids.
    irregular = FDA.MultivariateFunctionalData([FDA.IrregularFunctionalData(
        [grid, grid, grid], [values[i, 1, :] for i in 1:3])])
    @test_throws ArgumentError smooth_bspline(irregular; nbasis=2, degree=1)
    missing_values = Matrix{Union{Missing,Float64}}(values[:, 1, :])
    missing_values[1, 1] = missing
    incomplete = FDA.MultivariateFunctionalData([FDA.GridFunctionalData(missing_values, grid)])
    @test_throws ArgumentError smooth_bspline(incomplete; nbasis=2, degree=1)
    different_grid = FDA.MultivariateFunctionalData([
        FDA.GridFunctionalData(values[:, 1, :], grid),
        FDA.GridFunctionalData(values[:, 2, :], grid .+ 0.01)])
    @test_throws ArgumentError smooth_bspline(different_grid; nbasis=2, degree=1)

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
