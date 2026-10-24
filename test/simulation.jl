@testset "Qiao generating models" begin
    for model in (:ar2, :alternating, :random)
        sim = simulate_qiao2019(MersenneTwister(120); n=12, p=25, ncomponents=3,
                                ntime=21, model, noise_sd=0)
        @test size(FDA.to_array(sim.data)) == (12, 25, 21)
        @test size(sim.scores) == (12, 25, 3)
        @test size(sim.precision) == (75, 75)
        @test issymmetric(sim.precision)
        @test isposdef(Symmetric(sim.precision))
        @test issymmetric(sim.adjacency)
        @test all(!sim.adjacency[j, j] for j in 1:25)
        @test sim.settings.paper_coefficient_dimension == false
        @test sim.settings.reproduction_status == :generating_design
        # Scores and precision use exactly the same node-major ordering.
        for j in 1:25, k in 1:25
            block = sim.precision[(3j-2):3j, (3k-2):3k]
            @test block ≈ block[1, 1] * Matrix{Float64}(I, 3, 3)
            if j != k
                @test sim.adjacency[j, k] == !iszero(block[1, 1])
            end
        end
        t = FDA.component(sim.data, 1).grid
        basis = hcat(ones(length(t)), sqrt(2) .* sin.(2pi .* t), sqrt(2) .* cos.(2pi .* t))
        @test FDA.to_array(sim.data)[:, 1, :] ≈ sim.scores[:, 1, :] * transpose(basis)
        repeated = simulate_qiao2019(MersenneTwister(120); n=12, p=25, ncomponents=3,
                                     ntime=21, model, noise_sd=0)
        @test FDA.to_array(repeated.data) == FDA.to_array(sim.data)
        @test repeated.precision == sim.precision
        if model == :ar2
            @test count(sim.adjacency) == 2 * (2 * 25 - 3)
            @test sim.precision[1, 4] == 0.4
            @test sim.precision[1, 7] == 0.2
        elseif model == :alternating
            @test count(sim.adjacency) == 2 * (17 + 7)
            @test all(iszero, sim.adjacency[11:20, :])
            @test !sim.adjacency[10, 11]
            @test sim.adjacency[21, 23]
        else
            degree = maximum(sum(sim.adjacency; dims=1))
            @test sim.precision[1, 1] ≈ 0.5degree + 0.1
        end
    end
    latent = simulate_qiao2019(MersenneTwister(4); n=2, p=2, ntime=10, noise_sd=0)
    noisy = simulate_qiao2019(MersenneTwister(4); n=2, p=2, ntime=10, noise_sd=0.5)
    @test latent.scores == noisy.scores
    @test FDA.to_array(latent.data) != FDA.to_array(noisy.data)
    @test latent.settings.paper_coefficient_dimension
    @test_throws ArgumentError simulate_qiao2019(MersenneTwister(1); model=:bad)
    @test_throws ArgumentError simulate_qiao2019(MersenneTwister(1); noise_sd=-1)
    @test_throws ArgumentError simulate_qiao2019(MersenneTwister(1); ntime=1)
    @test_throws ArgumentError simulate_qiao2019(MersenneTwister(1); random_diagonal_margin=0)
end

@testset "Undirected graph metrics" begin
    truth = falses(4, 4)
    estimate = falses(4, 4)
    truth[1, 2] = truth[2, 1] = true
    truth[2, 3] = truth[3, 2] = true
    estimate[1, 2] = estimate[2, 1] = true
    estimate[1, 4] = estimate[4, 1] = true
    m = graph_metrics(estimate, truth)
    @test (m.tp, m.fp, m.tn, m.fn) == (1, 1, 3, 1)
    @test m.tpr == m.recall == m.precision == m.f1 == 0.5
    @test m.fpr == 0.25
    @test m.specificity == 0.75
    @test m.accuracy ≈ 2 / 3
    empty = graph_metrics(falses(4, 4), falses(4, 4))
    @test isnan(empty.precision)
    @test isnan(empty.recall)
    @test isnan(empty.f1)
    @test empty.fpr == 0
    @test isnan(graph_metrics(falses(1, 1), falses(1, 1)).accuracy)
    @test_throws ArgumentError graph_metrics([0 1; 0 0], falses(2, 2))
    @test_throws ArgumentError graph_metrics([0 2; 2 0], falses(2, 2))
    @test_throws DimensionMismatch graph_metrics(falses(3, 3), falses(2, 2))
end
