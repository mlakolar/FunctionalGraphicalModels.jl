@testset "Weighted FPCA and held-out projection" begin
    grid = [0.0, 0.1, 0.25, 0.5, 0.8, 1.0]
    w = quadrature_weights(grid)
    f1 = ones(length(grid))
    f1 ./= sqrt(dot(w, f1 .^ 2))
    f2 = grid .- dot(w, grid) / sum(w)
    f2 ./= sqrt(dot(w, f2 .^ 2))
    a = [-1.0, -1.0, 1.0, 1.0]
    b = [-1.0, 1.0, -1.0, 1.0]
    values = zeros(4, 2, length(grid))
    for i in 1:4
        values[i, 1, :] .= 2 .+ sqrt(2) .* a[i] .* f1 .+ 0.5 .* b[i] .* f2
        values[i, 2, :] .= -1 .+ 0.3 .* a[i] .* f1 .+ 2 .* b[i] .* f2
    end
    d = FunctionalGraphicalModels._dense_data(values, grid)
    for basis in (:marginal, :pooled)
        fit = FDA.fit(FDA.FPCA, d; ncomponents=2, maxK=2, basis)
        scores = FDA.score_array(fit)
        @test size(scores) == (4, 2, 2)
        @test all(b -> size(b.values) == (6, 2), fit.bases)
        @test all(v -> length(v) == 2, fit.eigenvalues)
        @test cat((reshape(x, 4, 1, 2) for x in FDA.predict(fit, d))...; dims=2) ≈ scores atol=1e-12
        @test FDA.to_array(FDA.reconstruct(fit)) ≈ values atol=1e-12
        @test vec(fit.means[1].values) ≈ fill(2, length(grid))
        @test vec(fit.means[2].values) ≈ fill(-1, length(grid))
        for j in 1:2
            phi = fit.bases[j].values
            @test transpose(phi) * (w .* phi) ≈ Matrix{Float64}(I, 2, 2) atol=1e-12
            @test fit.eigenvalues[j] ≈ vec(sum(abs2, scores[:, j, :]; dims=1)) ./ 3
        end
        if basis == :pooled
            @test fit.bases[1].values == fit.bases[2].values
        else
            @test fit.eigenvalues[1] ≈ (4 / 3) .* [2, 0.25] atol=1e-12
            @test fit.eigenvalues[2] ≈ (4 / 3) .* [4, 0.09] atol=1e-12
        end
        # A new-subject shift must be projected against the training mean.
        heldout = FunctionalGraphicalModels._dense_data(values[1:1, :, :] .+ 3, grid)
        new_scores = FDA.predict(fit, heldout)
        @test FDA.to_array(FDA.reconstruct(fit, new_scores)) ≈ FDA.to_array(heldout) atol=1e-12
        @test_throws ArgumentError FDA.predict(fit, FunctionalGraphicalModels._dense_data(values, grid .+ 0.1))
        @test_throws DimensionMismatch FDA.predict(fit, FunctionalGraphicalModels._dense_data(values[:, 1:1, :], grid))
        @test_throws DimensionMismatch FDA.reconstruct(fit, [zeros(4, 1), zeros(4, 1)])
    end
    @test_throws ArgumentError FDA.fit(FDA.FPCA, d; ncomponents=0)
    @test_throws ArgumentError FDA.fit(FDA.FPCA, d; ncomponents=4)
    @test_throws ArgumentError FDA.fit(FDA.FPCA, d; ncomponents=2, basis=:invalid)
    @test_throws ArgumentError FDA.fit(FDA.FPCA, FunctionalGraphicalModels._dense_data(values[1:1, :, :], grid); ncomponents=1)
end
