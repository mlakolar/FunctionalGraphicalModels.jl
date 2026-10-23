@testset "Observed-time spline and FPCA selection" begin
    grid = collect(range(0.0, 1.0; length=17))
    amplitudes = [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0]
    values = zeros(6, 2, length(grid))
    for i in 1:6
        values[i, 1, :] .= 1 .+ 2 .* grid .+ amplitudes[i] .* (0.5 .+ grid)
        values[i, 2, :] .= -1 .+ grid .+ amplitudes[i] .* (2 .- grid)
    end
    data = DenseFunctionalData(values, grid)
    fold_ids = [mod1(t, 4) for t in eachindex(grid)]
    selection = select_fgm_components(data; nbasis_grid=[2, 3, 20],
        ncomponents_grid=[1, 2, 6], degree=1, fold_ids)
    @test selection.nbasis in (2, 3)
    @test selection.ncomponents == 1
    @test selection.fold_ids == fold_ids
    @test size(selection.losses) == (3, 3)
    @test size(selection.fold_losses) == (3, 3, 4)
    @test maximum(selection.losses[1:2, 1]) < 1e-20
    @test all(isinf, selection.losses[:, 2:3])
    @test all(isinf, selection.losses[3, :])
    @test all(!isempty, selection.reasons[:, 2:3])
    @test !isempty(selection.reasons[3, 1])
    @test reconstruct(selection.fpca).values ≈ values atol=1e-11
    @test selection.smoothed.values ≈ values atol=1e-11

    # A manually computed, weighted aggregate of held-out errors agrees with
    # the reported CV criterion for one fixed (L,M) candidate.
    single = select_fgm_components(data; nbasis_grid=[2], ncomponents_grid=[1],
                                    degree=1, fold_ids)
    squared_error = 0.0
    for fold in 1:4
        train = findall(!=(fold), fold_ids)
        test = findall(==(fold), fold_ids)
        fit = FunctionalGraphicalModels._fgm_fold_prediction(data, train;
            nbasis=2, ncomponents=1, degree=1)
        squared_error += sum(abs2, fit.prediction.values[:, :, test] .- values[:, :, test])
    end
    @test single.losses[1, 1] ≈ squared_error / length(values) atol=1e-25

    # Changing only held-out outcomes cannot alter curve fits, FPCA means or
    # predictions in the training-only primitive.
    train = findall(!=(1), fold_ids)
    heldout = findall(==(1), fold_ids)
    original = FunctionalGraphicalModels._fgm_fold_prediction(data, train;
        nbasis=2, ncomponents=1, degree=1)
    perturbed_values = copy(values)
    perturbed_values[:, :, heldout] .+= 100
    perturbed = FunctionalGraphicalModels._fgm_fold_prediction(
        DenseFunctionalData(perturbed_values, grid), train;
        nbasis=2, ncomponents=1, degree=1)
    @test original.prediction.values == perturbed.prediction.values
    @test original.smoothed.values == perturbed.smoothed.values
    @test original.fpca.mean == perturbed.fpca.mean
    @test original.fpca.scores == perturbed.fpca.scores

    randomized1 = select_fgm_components(data; nbasis_grid=[2], ncomponents_grid=[1],
        degree=1, folds=4, rng=MersenneTwister(81))
    randomized2 = select_fgm_components(data; nbasis_grid=[2], ncomponents_grid=[1],
        degree=1, folds=4, rng=MersenneTwister(81))
    @test randomized1.fold_ids == randomized2.fold_ids
    @test randomized1.losses == randomized2.losses
    fold_sizes = [count(==(fold), randomized1.fold_ids) for fold in 1:4]
    @test maximum(fold_sizes) - minimum(fold_sizes) <= 1
    @test_throws ArgumentError select_fgm_components(data; nbasis_grid=[20], ncomponents_grid=[1], degree=1, fold_ids)
    @test_throws ArgumentError select_fgm_components(data; nbasis_grid=Int[], fold_ids)
    @test_throws ArgumentError select_fgm_components(data; nbasis_grid=[2, 2], degree=1, fold_ids)
    @test_throws DimensionMismatch select_fgm_components(data; fold_ids=[1, 2])
    @test_throws ArgumentError select_fgm_components(data; fold_ids=ones(Int, length(grid)))
    @test_throws ArgumentError select_fgm_components(data; fold_ids=[isodd(t) ? 1 : 3 for t in eachindex(grid)])
end
