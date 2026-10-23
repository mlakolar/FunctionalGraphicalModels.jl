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
    d = DenseFunctionalData(values, grid)
    for basis in (:marginal, :pooled)
        fit = fit_fpca(d; ncomponents=2, basis)
        @test size(fit.scores) == (4, 2, 2)
        @test size(fit.eigenfunctions) == (6, 2, 2)
        @test size(fit.eigenvalues) == (2, 2)
        @test transform_fpca(fit, d) ≈ fit.scores atol=1e-12
        @test reconstruct(fit).values ≈ values atol=1e-12
        @test fit.mean[1, :] ≈ fill(2, length(grid))
        @test fit.mean[2, :] ≈ fill(-1, length(grid))
        for j in 1:2
            phi = fit.eigenfunctions[:, :, j]
            @test transpose(phi) * (w .* phi) ≈ Matrix{Float64}(I, 2, 2) atol=1e-12
            @test fit.eigenvalues[:, j] ≈ vec(sum(abs2, fit.scores[:, j, :]; dims=1)) ./ 4
        end
        if basis == :pooled
            @test fit.eigenfunctions[:, :, 1] == fit.eigenfunctions[:, :, 2]
        else
            @test fit.eigenvalues[:, 1] ≈ [2, 0.25] atol=1e-12
            @test fit.eigenvalues[:, 2] ≈ [4, 0.09] atol=1e-12
        end
        # A new-subject shift must be projected against the training mean.
        heldout = DenseFunctionalData(values[1:1, :, :] .+ 3, grid)
        new_scores = transform_fpca(fit, heldout)
        @test reconstruct(fit, new_scores).values ≈ heldout.values atol=1e-12
        @test_throws ArgumentError transform_fpca(fit, DenseFunctionalData(values, grid .+ 0.1))
        @test_throws DimensionMismatch transform_fpca(fit, DenseFunctionalData(values[:, 1:1, :], grid))
        @test_throws DimensionMismatch reconstruct(fit, zeros(4, 2, 1))
    end
    @test_throws ArgumentError fit_fpca(d; ncomponents=0)
    @test_throws ArgumentError fit_fpca(d; ncomponents=4)
    @test_throws ArgumentError fit_fpca(d; ncomponents=2, basis=:invalid)
    @test_throws ArgumentError fit_fpca(DenseFunctionalData(values[1:1, :, :], grid); ncomponents=1)
end
