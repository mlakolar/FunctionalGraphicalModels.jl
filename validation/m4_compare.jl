#!/usr/bin/env julia

# Reproduce the migration comparison on one Julia version. The original checkout
# is read only: package resolution and execution happen in a temporary copy.
# Usage: julia --project=. validation/m4_compare.jl NEW_OUTPUT_DIRECTORY [BASE_CHECKOUT]
# The caller bootstraps this checkout with scripts/instantiate.jl first.
using LinearAlgebra
using Serialization
using SHA

const M4_BASE_COMMIT = "6f0eac3d5cbe0fad93e88dcad28f389d89b61495"

function source_digest(root)
    files = [joinpath(root, "Project.toml")]
    for (directory, _, names) in walkdir(joinpath(root, "src"))
        append!(files, (joinpath(directory, name) for name in names if endswith(name, ".jl")))
    end
    io = IOBuffer()
    for file in sort(files)
        write(io, relpath(file, root), '\0', read(file), '\0')
    end
    bytes2hex(sha256(take!(io)))
end

csvcell(x) = string('"', replace(string(x), '"' => "\"\"", '\n' => ' ', '\r' => ' '), '"')
csvrow(io, xs...) = println(io, join(csvcell.(xs), ','))

function relative_max_change(before, after)
    before === nothing && return "NA"
    size(before) == size(after) || throw(DimensionMismatch("snapshot array dimensions changed"))
    maximum(abs, after - before) / max(maximum(abs, before), eps(Float64))
end

function aligned_arrays(before, after)
    before.coordinates == after.coordinates || error("comparison coordinates changed")
    after.scores === nothing && return (nothing, after.covariance, after.precision, after.coefficients)
    size(before.bases) == size(after.bases) || throw(DimensionMismatch("FPCA bases changed dimensions"))
    _, m, p = size(before.bases)
    signs = ones(p, m)
    for j in 1:p, l in 1:m
        overlap = dot(before.weights .* before.bases[:, l, j], after.bases[:, l, j])
        signs[j, l] = overlap < 0 ? -1.0 : 1.0
    end
    scores = copy(after.scores)
    covariance = copy(after.covariance)
    precision = after.precision === nothing ? nothing : copy(after.precision)
    coefficients = after.coefficients === nothing ? nothing : copy(after.coefficients)
    if before.coordinates == :target_blocks
        for target in 1:p, l in 1:m
            scores[:, :, l, target] .*= signs[target, l]
        end
        for target in 1:p
            orientation = repeat(signs[target, :], p)
            covariance[:, :, target] .*= orientation * transpose(orientation)
            # Every predictor and the response use this target's basis.
            # Coefficients have predictor components in rows, responses in columns.
            predictor_orientation = repeat(signs[target, :], p - 1)
            coefficients[:, :, target] .*= predictor_orientation * transpose(signs[target, :])
        end
    else
        for j in 1:p, l in 1:m
            scores[:, j, l] .*= signs[j, l]
        end
        if before.coordinates == :node_blocks
            orientation = vec(transpose(signs))
            covariance .*= orientation * transpose(orientation)
            precision .*= orientation * transpose(orientation)
        elseif before.coordinates == :components
            for l in 1:m
                orientation = signs[:, l]
                covariance[:, :, l] .*= orientation * transpose(orientation)
                precision[:, :, l] .*= orientation * transpose(orientation)
            end
        elseif before.coordinates != :pointwise
            error("unknown comparison coordinates $(before.coordinates)")
        end
    end
    scores, covariance, precision, coefficients
end

function comparison_main(args=ARGS)
    1 <= length(args) <= 2 || error("usage: julia --project=. validation/m4_compare.jl NEW_OUTPUT_DIRECTORY [BASE_CHECKOUT]")
    root = normpath(joinpath(@__DIR__, ".."))
    reference = length(args) == 2 ? abspath(args[2]) :
        joinpath(homedir(), "projects", "juliaPkg", "fgm-ref-bsplinekit")
    base_revision = readchomp(`git -C $reference rev-parse HEAD`)
    base_revision == M4_BASE_COMMIT || error("base checkout must be $M4_BASE_COMMIT; got $base_revision")
    current_revision = readchomp(`git -C $root rev-parse HEAD`)
    current_dirty = !isempty(readchomp(`git -C $root status --porcelain --untracked-files=normal`))
    output = abspath(args[1])
    ispath(output) && error("output already exists: $output")
    mkpath(output)
    # Snapshot the worker so both revisions execute identical comparison code.
    worker = joinpath(output, "worker.jl")
    cp(joinpath(@__DIR__, "m4_compare_worker.jl"), worker)
    base_output = joinpath(output, "base.jls")
    current_output = joinpath(output, "current.jls")
    executable = Base.julia_cmd()
    mktempdir() do temporary
        copied_base = joinpath(temporary, "fgm-base")
        cp(reference, copied_base)
        # Never use the copied worktree's pointer to mutate the reference git state.
        rm(joinpath(copied_base, ".git"); recursive=true)
        open(joinpath(output, "base-instantiate.log"), "w") do io
            run(pipeline(`$executable --startup-file=no --project=$copied_base $(joinpath(copied_base, "scripts", "instantiate.jl"))`; stdout=io, stderr=io))
        end
        open(joinpath(output, "base-worker.log"), "w") do io
            run(pipeline(`$executable --startup-file=no --project=$copied_base $worker $base_output`; stdout=io, stderr=io))
        end
    end
    current_source_digest = source_digest(root)
    open(joinpath(output, "current-worker.log"), "w") do io
        run(pipeline(`$executable --startup-file=no --project=$root $worker $current_output`; stdout=io, stderr=io))
    end
    source_digest(root) == current_source_digest || error("current source changed while the worker was running; rerun the comparison")
    before, after = deserialize(base_output), deserialize(current_output)
    before.julia_version == after.julia_version || error("comparison must use one Julia version")
    Set(keys(before.rows)) == Set(keys(after.rows)) || error("comparison fixture names changed")
    table = joinpath(output, "comparison.csv")
    open(table, "w") do io
        println(io, "# generated_by=validation/m4_compare.jl; worker_sha256=", bytes2hex(sha256(read(worker))))
        println(io, "# Julia=", after.julia_version, "; base_commit=", base_revision,
            "; current_HEAD=", current_revision, current_dirty ? "+dirty" : "",
            "; current_source_sha256=", current_source_digest)
        println(io, "# base_dependencies=", repr(before.dependencies))
        println(io, "# current_dependencies=", repr(after.dependencies))
        println(io, "# relative_change=max(abs(after_aligned-before))/max(max(abs(before)),eps(Float64)); alignment uses weighted basis inner products")
        println(io, "# covariance uses likelihood denominator n; adjacency counts are upper-triangle edges, summed over time for pointwise fits")
        println(io, "# neighborhood scores/covariances use fitted means and target_bases; coefficients are the estimator's own predictor-row/response-column blocks, sign-aligned in each target's basis")
        for name in sort(collect(keys(after.rows)))
            println(io, "# ", name, " fixture=", after.rows[name].fixture)
            println(io, "# ", name, " estimator_settings=", after.rows[name].settings)
            println(io, "# ", name, " base_fpca_settings=", before.rows[name].fpca_settings)
            println(io, "# ", name, " current_fpca_settings=", after.rows[name].fpca_settings)
        end
        csvrow(io, "estimator", "base_edges", "current_edges", "adjacency_equal",
            "scores_relative_change", "covariance_relative_change", "precision_relative_change",
            "coefficients_relative_change",
            "eigenvalues_relative_change", "eigenvalues_after_dense_rescaling_relative_change",
            "base_noise_variance", "current_noise_variance")
        for name in sort(collect(keys(after.rows)))
            a, b = before.rows[name], after.rows[name]
            scores, covariance, precision, coefficients = aligned_arrays(a, b)
            dense_rescaling = name in ("glasso_dense", "neighborhood_dense", "doubly_dense",
                "separable_dense", "separable_pooled_rank")
            adjusted_eigenvalues = b.eigenvalues === nothing ? nothing :
                dense_rescaling ? b.eigenvalues .* ((b.nsubjects - 1) / b.nsubjects) : b.eigenvalues
            csvrow(io, name, count(a.adjacency) ÷ 2, count(b.adjacency) ÷ 2,
                a.adjacency == b.adjacency, relative_max_change(a.scores, scores),
                relative_max_change(a.covariance, covariance),
                relative_max_change(a.precision, precision),
                relative_max_change(a.coefficients, coefficients),
                relative_max_change(a.eigenvalues, b.eigenvalues),
                relative_max_change(a.eigenvalues, adjusted_eigenvalues),
                a.noise_variance === nothing ? "NA" : repr(a.noise_variance),
                b.noise_variance === nothing ? "NA" : repr(b.noise_variance))
        end
    end
    println("Comparison written to $table")
end

if abspath(PROGRAM_FILE) == @__FILE__
    comparison_main()
end
