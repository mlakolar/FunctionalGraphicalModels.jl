using Test
include("m4_compare.jl")

@testset "Neighborhood comparison sign alignment" begin
    n, p, m = 4, 3, 2
    bases = repeat(Matrix{Float64}(I, m, m), 1, 1, p)
    scores = reshape(collect(1.0:n*p*m*p), n, p, m, p)
    covariance = cat(([1.0 0.2; 0.2 2.0] .+ j for j in 1:p)...; dims=3)
    # Node-major predictor blocks, with nonzero cross-component coefficients.
    coefficients = reshape(collect(1.0:(p-1)*m*m*p), (p-1)*m, m, p)
    before = (; coordinates=:target_blocks, bases, weights=ones(m), scores,
        covariance=repeat(covariance, p, p, 1), precision=nothing, coefficients)
    after = deepcopy(before)
    signs = [-1.0 1.0; 1.0 -1.0; -1.0 -1.0]
    for target in 1:p
        D = Diagonal(signs[target, :])
        predictors = kron(Matrix{Float64}(I, p-1, p-1), D)
        all_nodes = kron(Matrix{Float64}(I, p, p), D)
        after.bases[:, :, target] .= before.bases[:, :, target] * D
        for node in 1:p
            after.scores[:, node, :, target] .= before.scores[:, node, :, target] * D
        end
        after.covariance[:, :, target] .= all_nodes * before.covariance[:, :, target] * all_nodes
        after.coefficients[:, :, target] .= predictors * before.coefficients[:, :, target] * D
    end
    unaligned = deepcopy(after)
    aligned_scores, aligned_covariance, precision, aligned_coefficients = aligned_arrays(before, after)
    @test aligned_scores == before.scores
    @test aligned_covariance == before.covariance
    @test precision === nothing
    @test aligned_coefficients == before.coefficients
    @test after == unaligned
    # A coefficient change must be measured even when bases and scores agree.
    after.coefficients[1, 2, 1] += 0.25
    _, _, _, changed_coefficients = aligned_arrays(before, after)
    @test relative_max_change(before.coefficients, changed_coefficients) ≈ 0.25 / maximum(abs, coefficients)
end
