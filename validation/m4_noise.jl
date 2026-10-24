#!/usr/bin/env julia

# Diagnose the estimated-noise migration without solving zero-noise CE systems.
# The base checkout is copied before Julia evaluates any of its source.
# Usage: julia --project=. validation/m4_noise.jl OUTPUT_CSV [BASE_CHECKOUT]
using FunctionalGraphicalModels
import FunctionalDataAnalysis as FDA
using LinearAlgebra
using Pkg
using SHA

module M4NoiseReference
    using LinearAlgebra
    # The diagnostic only calls the base sparse FPCA routine; this type permits
    # loading its unused dense graph method without loading the base package.
    struct DenseFunctionalData end
end

const M4_NOISE_BASE_COMMIT = "6f0eac3d5cbe0fad93e88dcad28f389d89b61495"


function noise_source_digest(root)
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

function noise_fixture(kind)
    n = kind === :constant ? 60 : 40
    times = Matrix{Vector{Float64}}(undef, n, 1)
    values = similar(times)
    if kind === :constant
        for i in 1:n
            times[i, 1] = collect(range(0.02mod(i, 3), 1 - 0.01mod(i, 4);
                length=5 + mod(i, 3)))
            values[i, 1] = fill(sqrt(2) * sin(2pi * i / n), length(times[i, 1]))
        end
        return (; times, values, grid=collect(range(0, 1; length=9)),
            bandwidth=0.3, noise_variance=0.04)
    end
    schedule = [0.02, 0.06, 0.14, 0.27, 0.44, 0.73, 0.98]
    for i in 1:n
        a, b = sqrt(2) * sin(2pi * i / n), 0.7sqrt(2) * cos(2pi * i / n)
        times[i, 1] = copy(schedule)
        values[i, 1] = a .+ b .* (2schedule .- 1)
    end
    (; times, values, grid=collect(range(0, 1; length=17)),
        bandwidth=0.2, noise_variance=0.05)
end

function integrated_excess(grid, diagonal, covariance)
    a, b = first(grid), last(grid)
    lower, upper = a + (b-a)/4, b - (b-a)/4
    interior = vcat(lower, grid[lower .< grid .< upper], upper)
    curve = FDA.GridFunctionalData(reshape(diagonal - diag(covariance), 1, :), grid)
    excess = vec(FDA.evaluate(curve, interior, FDA.PiecewiseLinear()))
    2dot(FDA.quadrature_weights(interior), excess)/(b-a)
end

function noise_diagnostic(kind)
    f = noise_fixture(kind)
    n, q = size(f.times, 1), length(f.grid)
    old = M4NoiseReference.fit_sparse_fpca(
        M4NoiseReference.SparseFunctionalData(f.times, f.values);
        grid=f.grid, bandwidth=f.bandwidth, noise_variance=f.noise_variance,
        ncomponents=1)
    data = FDA.MultivariateFunctionalData([
        FDA.IrregularFunctionalData(f.times[:, 1], f.values[:, 1])])
    current = fit_doubly(data; lambda=0.1, grid=f.grid, bandwidth=f.bandwidth,
        noise_variance=f.noise_variance, ncomponents=1)
    node = FDA.component(current.fpca, 1)
    covariance = node.settings.covariance_input
    # Evaluate the old source's exact diagonal smoother with its own fitted
    # means, matching the deleted fit_sparse_fpca noise-estimation formula.
    pooled_times, squares, pointweights = Float64[], Float64[], Float64[]
    for i in 1:n
        ti = f.times[i, 1]
        centered = f.values[i, 1] .- [M4NoiseReference._dfgm_interpolate(
            f.grid, @view(old.mean[1, :]), t) for t in ti]
        append!(pooled_times, ti)
        append!(squares, centered .^ 2)
        append!(pointweights, fill(1/(n*length(ti)), length(ti)))
    end
    old_diagonal = [M4NoiseReference._dfgm_local_linear(pooled_times, squares,
        pointweights, t, f.bandwidth) for t in f.grid]
    middle = max(1, ceil(Int, q/4)):min(q, floor(Int, 3q/4))
    old_raw = sum(old_diagonal[t] - old.covariances[t, t, 1] for t in middle)/length(middle)
    old_integrated_psd_raw = integrated_excess(f.grid, old_diagonal, old.covariances[:, :, 1])
    current_integrated_psd_raw = integrated_excess(f.grid,
        covariance.noise_details.diagonal, node.settings.full_projected_covariance)
    current_integrated_raw = covariance.noise_details.raw_variance
    # Keep unclipped raw estimates: truncation at zero could conceal a change.
    (; fixture=string(kind), base_noise=max(old_raw, 0),
        integrated_base_psd_noise=max(old_integrated_psd_raw, 0),
        integrated_current_psd_noise=max(current_integrated_psd_raw, 0),
        fda_noise=covariance.noise_variance,
        base_raw=old_raw, integrated_base_psd_raw=old_integrated_psd_raw,
        integrated_current_psd_raw=current_integrated_psd_raw,
        fda_raw=current_integrated_raw,
        integration_delta_raw=old_integrated_psd_raw-old_raw,
        raw_vs_psd_delta_raw=current_integrated_raw-current_integrated_psd_raw,
        smoother_diagonal_max_change=maximum(abs, old_diagonal-covariance.noise_details.diagonal),
        projected_covariance_max_change=maximum(abs,
            old.covariances[:, :, 1]-node.settings.full_projected_covariance),
        raw_covariance_min_eigenvalue=minimum(eigvals(Symmetric(
            sqrt.(node.weights).*covariance.surface.*transpose(sqrt.(node.weights))))))
end

function noise_main(args=ARGS)
    1 <= length(args) <= 2 || error("usage: julia --project=. validation/m4_noise.jl OUTPUT_CSV [BASE_CHECKOUT]")
    reference = length(args) == 2 ? abspath(args[2]) :
        joinpath(homedir(), "projects", "juliaPkg", "fgm-ref-bsplinekit")
    revision = readchomp(`git -C $reference rev-parse HEAD`)
    revision == M4_NOISE_BASE_COMMIT || error("base checkout must be $M4_NOISE_BASE_COMMIT; got $revision")
    root = normpath(joinpath(@__DIR__, ".."))
    output = abspath(args[1])
    ispath(output) && error("output already exists: $output")
    rows = mktempdir() do temporary
        copied = joinpath(temporary, "fgm-base")
        cp(reference, copied)
        rm(joinpath(copied, ".git"); recursive=true)
        Base.include(M4NoiseReference, joinpath(copied, "src", "doubly.jl"))
        # include introduces methods after this function was compiled.
        [Base.invokelatest(noise_diagnostic, kind) for kind in (:constant, :full_modes)]
    end
    mkpath(dirname(output))
    open(output, "w") do io
        println(io, "# generated_by=validation/m4_noise.jl; script_sha256=", bytes2hex(sha256(read(@__FILE__))))
        println(io, "# Julia=", VERSION, "; base_commit=", revision,
            "; current_HEAD=", readchomp(`git -C $root rev-parse HEAD`),
            "; current_source_sha256=", noise_source_digest(root))
        for info in values(Pkg.dependencies())
            info.name in ("FunctionalDataAnalysis", "KernelSmoothing") || continue
            println(io, "# ", info.name, " version=", info.version,
                "; git_revision=", info.git_revision, "; git_source=", info.git_source,
                "; tree_hash=", info.tree_hash)
        end
        println(io, "# Both fits use the fixture's explicit positive noise; estimated-noise formulas are evaluated without CE solves.")
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(string.(values(row)), ','))
        end
    end
    println("Noise diagnostic written to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    noise_main()
end
