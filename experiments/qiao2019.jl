#!/usr/bin/env julia

# Run on an authorized compute allocation after bootstrapping the package.
# This driver fixes the spline dimension and FPCA truncation. It does not
# implement the original paper's joint five-fold choice of L and M or its
# complete competitor comparison, even when --paper is selected.

using FunctionalGraphicalModels
using LinearAlgebra
using Random
using TOML

function usage(io=stdout)
    println(io, """
    julia --project=. experiments/qiao2019.jl [--smoke | --paper] [options]

    --smoke                 Small pilot (default): n=40, p=6, T=41, M=3,
                            one replicate and six penalties.
    --paper                 Paper-size pilot: n=100, p=50, T=100, M=5,
                            100 replicates and 20 penalties. Tuning remains fixed.
    --n INT                 Subjects.
    --p INT                 Nodes.
    --ntime INT             Common-grid observations per curve.
    --components INT        Retained FPCA components; generating dimension is 5.
    --nbasis INT            Cubic B-spline dimension (default 7).
    --replicates INT         Number of independent datasets.
    --model NAME            ar2 (default), alternating, random, or all.
    --seed INT              First replicate seed (default 20261024).
    --nlambda INT           Number of logarithmically spaced penalties.
    --lambda-min-ratio REAL Smallest/largest penalty (default 0.001).
    --maxiter INT           Solver iteration limit (default 10000).
    --threshold REAL        Graph extraction threshold (default 1e-6).
    --out PATH              New output directory (default results/qiao2019-MODE).
    --help                  Print this message.

    Output: metrics.csv, edges.csv, reference_endpoints.csv, settings.toml,
    and failures.log. Set FGM_SOURCE_REV to the tested source revision before
    running from an archive. Failed and unconverged fits remain in metrics.csv
    with NaN graph metrics and a nonzero process exit status. Reference endpoints
    are explicit empty/complete graphs, not estimates or extra fitted penalties.
    """)
end

function options(args)
    any(==("--help"), args) && return nothing
    mode = "smoke"
    modes = filter(x -> x in ("--smoke", "--paper"), args)
    length(modes) <= 1 || throw(ArgumentError("choose one of --smoke and --paper"))
    !isempty(modes) && (mode = modes[1][3:end])
    opts = Dict{String,Any}("mode" => mode, "n" => (mode == "paper" ? 100 : 40),
        "p" => (mode == "paper" ? 50 : 6), "ntime" => (mode == "paper" ? 100 : 41),
        "components" => (mode == "paper" ? 5 : 3), "nbasis" => 7,
        "replicates" => (mode == "paper" ? 100 : 1), "model" => "ar2",
        "seed" => 20261024, "nlambda" => (mode == "paper" ? 20 : 6),
        "lambda-min-ratio" => 0.001, "maxiter" => 10000, "threshold" => 1e-6,
        "out" => joinpath("results", "qiao2019-" * mode))
    integer_keys = ("n", "p", "ntime", "components", "nbasis", "replicates",
                    "seed", "nlambda", "maxiter")
    real_keys = ("lambda-min-ratio", "threshold")
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--smoke", "--paper")
            i += 1
            continue
        end
        startswith(arg, "--") || throw(ArgumentError("unexpected argument $arg"))
        pair = split(arg[3:end], '='; limit=2)
        key = pair[1]
        haskey(opts, key) && key != "mode" || throw(ArgumentError("unknown option --$key"))
        if length(pair) == 2
            value = pair[2]
        else
            i < length(args) || throw(ArgumentError("missing value for --$key"))
            i += 1
            value = args[i]
        end
        opts[key] = key in integer_keys ? parse(Int, value) :
                    key in real_keys ? parse(Float64, value) : value
        i += 1
    end
    opts["n"] >= 2 && opts["p"] >= 2 && opts["ntime"] >= 2 ||
        throw(ArgumentError("n, p and ntime must be at least two"))
    1 <= opts["components"] <= min(opts["n"]-1, opts["ntime"], opts["nbasis"]) ||
        throw(ArgumentError("components must fit within subject, grid and spline ranks"))
    4 <= opts["nbasis"] <= opts["ntime"] ||
        throw(ArgumentError("cubic nbasis must be between 4 and ntime"))
    opts["replicates"] > 0 && opts["nlambda"] >= 2 && opts["maxiter"] > 0 ||
        throw(ArgumentError("replicates/maxiter must be positive and nlambda at least two"))
    opts["seed"] >= 0 || throw(ArgumentError("seed must be nonnegative"))
    0 < opts["lambda-min-ratio"] < 1 || throw(ArgumentError("lambda-min-ratio must lie in (0,1)"))
    isfinite(opts["threshold"]) && opts["threshold"] >= 0 ||
        throw(ArgumentError("threshold must be finite and nonnegative"))
    opts["model"] in ("ar2", "alternating", "random", "all") ||
        throw(ArgumentError("unknown model"))
    opts["out"] = abspath(opts["out"])
    opts
end

csvcell(x) = string('"', replace(string(x), '"' => "\"\"", '\n' => ' ', '\r' => ' '), '"')
csvrow(io, xs...) = (println(io, join(csvcell.(xs), ',')); flush(io))

function penalty_bounds(data, fpca)
    n, p, nt = size(data.values)
    m = size(fpca.scores, 3)
    S = score_covariance(fpca.scores)
    glasso_bound = maximum(norm(view(S, (j-1)*m+1:j*m, (k-1)*m+1:k*m))
                           for j in 1:p for k in j+1:p)
    # EJS requires projecting every predictor onto the target's common basis.
    # Marginal-score covariance bounds would describe a different regression.
    target_scores = zeros(n, p, m)
    neighborhood_bound = 0.0
    for target in 1:p
        projection = fpca.weights .* fpca.eigenfunctions[:, :, target]
        for node in 1:p
            centered = data.values[:, node, :] .- transpose(fpca.mean[node, :])
            target_scores[:, node, :] .= centered * projection
        end
        neighborhood_bound = max(neighborhood_bound,
            neighborhood_lambda_max(target_scores, target))
    end
    # A strict margin covers roundoff at the theoretical zero-graph endpoint.
    # Positive fallback also permits a path for exactly uncorrelated data.
    (; fglasso=1.05 * max(glasso_bound, eps(Float64)),
       neighborhood=1.05 * max(neighborhood_bound, eps(Float64)))
end

invalid_metrics() = (; tp=NaN, fp=NaN, tn=NaN, fn=NaN, tpr=NaN, fpr=NaN,
                      precision=NaN, recall=NaN, f1=NaN, accuracy=NaN)

function metric_row(io, model, replicate, seed, method, index, lambda, bound,
                    converged, kkt, primal, dual, iterations, edges, metrics,
                    status, message)
    csvrow(io, model, replicate, seed, method, index, lambda, bound, converged,
           kkt, primal, dual, iterations, edges, metrics.tp, metrics.fp, metrics.tn,
           metrics.fn, metrics.tpr, metrics.fpr, metrics.precision, metrics.recall,
           metrics.f1, metrics.accuracy, status, message)
end

function edge_rows(io, adjacency, model, replicate, seed, method, index, status)
    for j in 2:size(adjacency, 1), k in 1:j-1
        adjacency[j, k] && csvrow(io, model, replicate, seed, method, index, k, j, status)
    end
end

function record_failure(io, context, err, trace)
    println(io, context)
    showerror(io, err)
    println(io)
    Base.show_backtrace(io, trace)
    println(io)
    flush(io)
end

function run_experiment(opts, arguments)
    out = opts["out"]
    files = ("metrics.csv", "edges.csv", "reference_endpoints.csv", "settings.toml", "failures.log")
    any(name -> ispath(joinpath(out, name)), files) &&
        throw(ArgumentError("output files already exist in $out; choose a new --out directory"))
    mkpath(out)
    settings = copy(opts)
    settings["arguments"] = collect(arguments)
    settings["julia_version"] = string(VERSION)
    settings["package_version"] = string(Base.pkgversion(FunctionalGraphicalModels))
    settings["source_revision"] = get(ENV, "FGM_SOURCE_REV", "not_recorded")
    settings["generating_components"] = 5
    settings["noise_sd"] = 0.5
    settings["spline_degree"] = 3
    settings["spline_penalty"] = 0.0
    settings["random_diagonal_margin"] = 0.1
    settings["random_diagonal_rule"] = "0.5 * maximum_degree + 0.1"
    settings["basis_convention"] = "1,sqrt(2)sin(2pi*t),sqrt(2)cos(2pi*t),sqrt(2)sin(4pi*t),sqrt(2)cos(4pi*t)"
    settings["covariance_normalization"] = "1/n"
    settings["penalty_grid"] = "log-spaced from 1.05*method_lambda_max to lambda-min-ratio times that value"
    settings["neighborhood_rule"] = "or"
    settings["methods"] = ["fglasso", "neighborhood"]
    settings["reproduction_status"] = "fixed-tuning pilot; original joint CV and competitors are not reproduced"
    settings["reference_endpoints"] = "empty and complete graphs; not fitted estimates; stored separately"
    settings["paper_source"] = "https://personal.lse.ac.uk/qiaox/qiao.links/research_papers/2017_FGM.pdf"
    open(joinpath(out, "settings.toml"), "w") do io
        TOML.print(io, settings; sorted=true)
    end
    metrics_io = open(joinpath(out, "metrics.csv"), "w")
    edges_io = open(joinpath(out, "edges.csv"), "w")
    endpoints_io = open(joinpath(out, "reference_endpoints.csv"), "w")
    failures_io = open(joinpath(out, "failures.log"), "w")
    failures = unconverged = successful = 0
    try
        csvrow(metrics_io, "model", "replicate", "seed", "method", "lambda_index",
            "lambda", "lambda_max_with_margin", "converged", "kkt_residual",
            "primal_residual", "dual_residual", "iterations", "estimated_edges",
            "tp", "fp", "tn", "fn", "tpr", "fpr", "precision", "recall", "f1",
            "accuracy", "status", "message")
        csvrow(edges_io, "model", "replicate", "seed", "method", "lambda_index", "node1", "node2", "status")
        csvrow(endpoints_io, "model", "replicate", "seed", "graph", "status", "tp", "fp", "tn", "fn", "tpr", "fpr")
        models = opts["model"] == "all" ? (:ar2, :alternating, :random) : (Symbol(opts["model"]),)
        for model in models, replicate in 1:opts["replicates"]
            seed = opts["seed"] + replicate - 1
            local sim, smoothed, fpca, bounds
            try
                sim = simulate_qiao2019(MersenneTwister(seed); n=opts["n"], p=opts["p"],
                    ncomponents=5, ntime=opts["ntime"], model, noise_sd=0.5)
                smoothed = smooth_bspline(sim.data; nbasis=opts["nbasis"])
                fpca = fit_fpca(smoothed; ncomponents=opts["components"])
                bounds = penalty_bounds(smoothed, fpca)
                empty_graph = falses(opts["p"], opts["p"])
                complete_graph = trues(opts["p"], opts["p"])
                for j in 1:opts["p"]
                    complete_graph[j, j] = false
                end
                for (name, graph) in (("empty", empty_graph), ("complete", complete_graph))
                    m = graph_metrics(graph, sim.adjacency)
                    csvrow(endpoints_io, model, replicate, seed, name, "reference_endpoint",
                           m.tp, m.fp, m.tn, m.fn, m.tpr, m.fpr)
                end
            catch err
                err isa InterruptException && rethrow()
                failures += 1
                message = sprint(showerror, err)
                record_failure(failures_io, "model=$model replicate=$replicate seed=$seed preprocessing", err, catch_backtrace())
                for method in (:fglasso, :neighborhood), index in 1:opts["nlambda"]
                    metric_row(metrics_io, model, replicate, seed, method, index,
                        NaN, NaN, false, NaN, NaN, NaN, 0, NaN, invalid_metrics(),
                        "preprocessing_error", message)
                end
                continue
            end
            for method in (:fglasso, :neighborhood)
                bound = getproperty(bounds, method)
                lambdas = exp.(range(log(bound), log(bound * opts["lambda-min-ratio"]); length=opts["nlambda"]))
                for (index, lambda) in enumerate(lambdas)
                    try
                        if method == :fglasso
                            fit = functional_glasso(fpca.scores, lambda;
                                maxiter=opts["maxiter"], threshold=opts["threshold"])
                            d = fit.diagnostics
                            converged, kkt, primal, dual, iterations = d.converged,
                                d.kkt_residual, d.primal_residual, d.dual_residual, d.iterations
                        else
                            fit = functional_neighborhood(smoothed, lambda;
                                ncomponents=opts["components"], maxiter=opts["maxiter"],
                                threshold=opts["threshold"], rule=:or)
                            ds = fit.diagnostics
                            converged = all(d -> d.converged, ds)
                            kkt = maximum(d -> d.kkt_residual, ds)
                            primal = maximum(d -> d.primal_residual, ds)
                            dual = maximum(d -> d.dual_residual, ds)
                            iterations = maximum(d -> d.iterations, ds)
                        end
                        status = converged ? "ok" : "nonconverged"
                        message = converged ? "" : "Graph metrics withheld; inspect diagnostics and saved edge support."
                        m = converged ? graph_metrics(fit.adjacency, sim.adjacency) : invalid_metrics()
                        metric_row(metrics_io, model, replicate, seed, method, index,
                            lambda, bound, converged, kkt, primal, dual, iterations,
                            count(fit.adjacency) ÷ 2, m, status, message)
                        edge_rows(edges_io, fit.adjacency, model, replicate, seed, method, index, status)
                        converged ? (successful += 1) : (unconverged += 1)
                    catch err
                        err isa InterruptException && rethrow()
                        failures += 1
                        message = sprint(showerror, err)
                        record_failure(failures_io,
                            "model=$model replicate=$replicate seed=$seed method=$method lambda=$lambda", err, catch_backtrace())
                        metric_row(metrics_io, model, replicate, seed, method, index,
                            lambda, bound, false, NaN, NaN, NaN, 0, NaN,
                            invalid_metrics(), "solver_error", message)
                    end
                end
            end
            println("model=$model replicate=$replicate seed=$seed completed; successful=$successful unconverged=$unconverged failures=$failures")
            flush(stdout)
        end
    finally
        close(metrics_io)
        close(edges_io)
        close(endpoints_io)
        close(failures_io)
    end
    println("Saved results to $out. Successful fits: $successful; unconverged: $unconverged; failures: $failures.")
    failures == 0 && unconverged == 0 ? 0 : 1
end

function main(args=ARGS)
    opts = options(args)
    if opts === nothing
        usage()
        return 0
    end
    run_experiment(opts, args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
