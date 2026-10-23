#!/usr/bin/env julia

# Original implementation of the five-method comparison in Qiao, Guo and
# James, Section 5.1. Run on an authorized compute allocation. The default is
# a bounded paper-size pilot, not the paper's 100-replicate result.

using FunctionalGraphicalModels
using CovarianceSelection
using LinearAlgebra
using Random
using TOML

const COMPARISON_METHODS = ("fglasso", "ALL", "NEVER", "HALF", "PCA")

function comparison_usage(io=stdout)
    println(io, """
    julia --project=. experiments/qiao2019_comparison.jl [options]

    Five-method Qiao et al. comparison with joint time-fold CV of spline
    dimension L and functional PCA dimension M. Defaults: n=100, p=50,
    T=100, five generating Fourier components, Gaussian noise SD=0.5.

    --replicates INT         Independent datasets (default 1; paper uses 100).
    --nlambda INT            Positive penalties per method (default 8).
    --model NAME             ar2 (default), alternating, random, or all.
    --p INT                  Nodes (default 50; paper uses 50,100,150).
    --seed INT               First dataset seed (default 20261024).
    --out PATH               New output directory (default results/qiao2019-comparison).
    --lambda-min-ratio REAL  Smallest fraction of lambda_max (default 0.01).
    --fgm-rho REAL           Relative FGM ADMM penalty (default 1; objective unchanged).
    --maxiter INT            ADMM iteration limit per fit (default 10000).
    --threshold REAL         Absolute precision support threshold (default 1e-6).
    --n INT                  Subjects (default 100).
    --ntime INT              Common-grid measurements (default 100).
    --folds INT              Time-fold CV count (default 5).
    --nbasis-grid GRID       L candidates, e.g. 5:10 or 5,7,9 (default 5:10).
    --components-grid GRID   M candidates (default 1:6; candidates M>L are excluded).
    --help                   Print this message.

    ALL is the intersection of timepoint supports; NEVER is their union;
    HALF requires strictly more than T/2 selections. These three use raw
    noisy measurements and one common scalar penalty across all timepoints.
    PCA uses its own L cross-validation with M fixed to one. Diagonals are
    unpenalized in every method. Scalar fits use CovarianceSelection.covsel!.

    Output: metrics.csv, solver_diagnostics.csv, edges.csv, cv_losses.csv,
    cv_folds.csv, auc.csv, reference_endpoints.csv, settings.toml, and failures.log.
    CV losses include invalid-candidate reasons; CV folds save each tuning
    method's exact timepoint-to-fold assignments.
    Unconverged/invalid fits have withheld graph metrics and cause exit status 1.
    A full fitted AUC is reported only if valid fitted empty and complete
    graphs are both observed. Otherwise only the observed-FPR-span area is
    reported. Empty/complete reference graphs are saved separately and are
    never inserted as fitted path points. This driver does not reproduce the
    paper's original random-number stream or establish agreement with its table.
    Coarse path areas are numerical diagnostics and do not establish method
    rankings; refine the penalty grid and resolve failed fits before comparisons.
    """)
end

function comparison_grid(value)
    text = String(value)
    if occursin(':',text)
        bounds = parse.(Int,split(text,':'))
        length(bounds) == 2 || throw(ArgumentError("grid ranges must use start:stop"))
        bounds[1] <= bounds[2] || throw(ArgumentError("grid ranges must be increasing"))
        return collect(bounds[1]:bounds[2])
    end
    unique(parse.(Int,split(text,',')))
end

function comparison_options(args)
    any(==("--help"),args) && return nothing
    opts = Dict{String,Any}("replicates"=>1,"nlambda"=>8,"model"=>"ar2",
        "p"=>50,"seed"=>20261024,"out"=>joinpath("results","qiao2019-comparison"),
        "lambda-min-ratio"=>0.01,"fgm-rho"=>1.0,"maxiter"=>10000,"threshold"=>1e-6,
        "n"=>100,"ntime"=>100,"folds"=>5,"nbasis-grid"=>collect(5:10),
        "components-grid"=>collect(1:6))
    intkeys = ("replicates","nlambda","p","seed","maxiter","n","ntime","folds")
    realkeys = ("lambda-min-ratio","fgm-rho","threshold")
    gridkeys = ("nbasis-grid","components-grid")
    i = 1
    while i <= length(args)
        startswith(args[i],"--") || throw(ArgumentError("unexpected argument $(args[i])"))
        pair = split(args[i][3:end],'='; limit=2)
        key = pair[1]
        haskey(opts,key) || throw(ArgumentError("unknown option --$key"))
        if length(pair) == 2
            value = pair[2]
        else
            i < length(args) || throw(ArgumentError("missing value for --$key"))
            i += 1
            value = args[i]
        end
        opts[key] = key in intkeys ? parse(Int,value) : key in realkeys ?
                    parse(Float64,value) : key in gridkeys ? comparison_grid(value) : value
        i += 1
    end
    opts["n"] >= 2 && opts["p"] >= 2 && opts["ntime"] >= 3 ||
        throw(ArgumentError("n,p must be at least two and ntime at least three"))
    opts["replicates"] >= 1 && opts["nlambda"] >= 2 && opts["maxiter"] >= 1 ||
        throw(ArgumentError("replicates/maxiter must be positive and nlambda at least two"))
    0 <= opts["seed"] <= typemax(Int)-opts["replicates"]-10_000_000 ||
        throw(ArgumentError("seed must be nonnegative and leave room for replicate/tuning seeds"))
    2 <= opts["folds"] <= opts["ntime"] || throw(ArgumentError("folds must lie between 2 and ntime"))
    0 < opts["lambda-min-ratio"] < 1 || throw(ArgumentError("lambda-min-ratio must lie in (0,1)"))
    isfinite(opts["fgm-rho"]) && opts["fgm-rho"] > 0 ||
        throw(ArgumentError("fgm-rho must be finite and positive"))
    isfinite(opts["threshold"]) && opts["threshold"] >= 0 ||
        throw(ArgumentError("threshold must be finite and nonnegative"))
    !isempty(opts["nbasis-grid"]) && all(>=(4),opts["nbasis-grid"]) ||
        throw(ArgumentError("cubic spline dimensions must be at least four"))
    !isempty(opts["components-grid"]) && all(>(0),opts["components-grid"]) ||
        throw(ArgumentError("component candidates must be positive"))
    maximum(opts["nbasis-grid"]) <= opts["ntime"]-cld(opts["ntime"],opts["folds"]) ||
        throw(ArgumentError("every CV training fold must have at least max(L) timepoints"))
    opts["model"] in ("ar2","alternating","random","all") || throw(ArgumentError("unknown model"))
    opts["out"] = abspath(opts["out"])
    opts
end

comparison_csvcell(x) = string('"',replace(string(x),'"'=>"\"\"",'\n'=>' ','\r'=>' '),'"')
comparison_csvrow(io,xs...) = (println(io,join(comparison_csvcell.(xs),',')); flush(io))
comparison_invalid_metrics() = (;tp=NaN,fp=NaN,tn=NaN,fn=NaN,tpr=NaN,fpr=NaN,
                                 precision=NaN,recall=NaN,f1=NaN,accuracy=NaN)

function comparison_record_failure(io,context,err,trace)
    println(io,context)
    showerror(io,err)
    println(io)
    Base.show_backtrace(io,trace)
    println(io)
    flush(io)
end

function comparison_covariance(A)
    X = Matrix{Float64}(A)
    X .-= sum(X; dims=1)/size(X,1)
    S = transpose(X)*X/size(X,1)
    all(isfinite,S) || error("nonfinite sample covariance")
    Matrix(Symmetric(S))
end

function comparison_scalar_bound(S)
    p = size(S,1)
    maximum(abs(S[j,k]) for j in 2:p for k in 1:j-1)
end

function comparison_block_bound(S,m)
    p = div(size(S,1),m)
    maximum(norm(view(S,((j-1)*m+1):(j*m),((k-1)*m+1):(k*m)))
            for j in 2:p for k in 1:j-1)
end

function comparison_scalar_graph(Z,threshold)
    p = size(Z,1)
    graph = falses(p,p)
    for j in 2:p,k in 1:j-1
        graph[j,k] = graph[k,j] = abs((Z[j,k]+Z[k,j])/2) > threshold
    end
    graph
end

function comparison_scalar_kkt(S,Z,lambda)
    isposdef(Symmetric(Z)) || return Inf
    gradient = S-inv(Symmetric(Z))
    worst = maximum(abs,diag(gradient))
    for j in 2:size(S,1),k in 1:j-1
        z = (Z[j,k]+Z[k,j])/2
        residual = iszero(z) ? max(abs(gradient[j,k])-lambda,0) :
                   abs(gradient[j,k]+lambda*sign(z))
        worst = max(worst,residual)
    end
    worst
end

function comparison_diagnostic(io,model,replicate,seed,method,time,index,fraction,
                                lambda,bound,converged,pd,kkt,primal,dual,iterations,status,message)
    comparison_csvrow(io,model,replicate,seed,method,time,index,fraction,lambda,bound,
        converged,pd,kkt,primal,dual,iterations,status,message)
end

function comparison_metrics(io,edges_io,adjacency,truth,model,replicate,seed,method,
                             index,fraction,lambda,bound,valid,status,message)
    metrics = valid ? graph_metrics(adjacency,truth) : comparison_invalid_metrics()
    estimated_edges = adjacency === nothing ? NaN : count(adjacency)÷2
    comparison_csvrow(io,model,replicate,seed,method,index,fraction,lambda,bound,valid,
        estimated_edges,metrics.tp,metrics.fp,metrics.tn,metrics.fn,metrics.tpr,
        metrics.fpr,metrics.precision,metrics.recall,metrics.f1,metrics.accuracy,status,message)
    if adjacency !== nothing
        for j in 2:size(adjacency,1),k in 1:j-1
            adjacency[j,k] && comparison_csvrow(edges_io,model,replicate,seed,method,
                                               index,k,j,status)
        end
    end
    (;valid,tpr=metrics.tpr,fpr=metrics.fpr,
      empty=valid && iszero(estimated_edges),
      complete=valid && estimated_edges == size(truth,1)*(size(truth,1)-1)÷2)
end

function comparison_auc(io,points,model,replicate,seed,method,nlambda)
    valid = filter(x -> x.valid && isfinite(x.fpr) && isfinite(x.tpr),points)
    unique_points = Dict{Float64,Float64}()
    for point in valid
        unique_points[point.fpr] = max(get(unique_points,point.fpr,-Inf),point.tpr)
    end
    xs = sort!(collect(keys(unique_points)))
    area = length(xs) >= 2 ? sum((xs[i+1]-xs[i])*
        (unique_points[xs[i+1]]+unique_points[xs[i]])/2 for i in 1:length(xs)-1) : NaN
    empty_observed = any(x -> x.empty,valid)
    complete_observed = any(x -> x.complete,valid)
    all_valid = length(valid) == nlambda
    full_auc = all_valid && empty_observed && complete_observed ? area : NaN
    status = !all_valid ? "incomplete_or_failed_path" :
             empty_observed && complete_observed ? "fitted_endpoints_observed" : "observed_span_only"
    comparison_csvrow(io,model,replicate,seed,method,length(valid),nlambda,
        isempty(xs) ? NaN : first(xs),isempty(xs) ? NaN : last(xs),
        empty_observed,complete_observed,area,full_auc,status)
end

function comparison_cv_rows(io,tuning,model,replicate,seed,method)
    for (a,L) in pairs(tuning.nbasis_grid),(b,M) in pairs(tuning.ncomponents_grid)
        loss = tuning.losses[a,b]
        comparison_csvrow(io,model,replicate,seed,method,L,M,loss,
            L == tuning.nbasis && M == tuning.ncomponents,
            isfinite(loss) ? "ok" : "ineligible_or_failed",tuning.reasons[a,b])
    end
end

function comparison_cv_fold_rows(io,tuning,model,replicate,seed,method)
    for (time_index,fold) in pairs(tuning.fold_ids)
        comparison_csvrow(io,model,replicate,seed,method,time_index,fold)
    end
end

function comparison_scalar_path!(diagnostics_io,failures_io,Ss,bound,fractions,
                                 model,replicate,seed,opts;method="timepoint_glasso")
    p = size(Ss[1],1)
    counts = zeros(Int,p,p,length(fractions))
    valid = trues(length(fractions))
    for (time,S) in pairs(Ss)
        X,Z,U = zeros(p,p),zeros(p,p),zeros(p,p)
        warm_rho = 1.0
        for (index,fraction) in pairs(fractions)
            lambda = fraction*bound
            try
                info = CovarianceSelection.ADMMInfo()
                options = CovarianceSelection.ADMMOptions(; ρ=warm_rho,
                    maxiter=opts["maxiter"],abstol=1e-7,reltol=1e-6,adaptive=true)
                CovarianceSelection.covsel!(X,Z,U,S,lambda;
                    options,penalize_diag=false,info)
                pd = all(isfinite,Z) && isposdef(Symmetric(Z))
                admissible = info.converged && pd
                kkt = pd ? comparison_scalar_kkt(S,Z,lambda) : Inf
                status = !info.converged ? "nonconverged" : !pd ? "nonpositive_precision" : "ok"
                message = admissible ? "" : "Graph metrics withheld for aggregate methods at this path point."
                comparison_diagnostic(diagnostics_io,model,replicate,seed,method,time,index,
                    fraction,lambda,bound,info.converged,pd,kkt,info.r_norm,info.s_norm,
                    info.iterations,status,message)
                if admissible
                    counts[:,:,index] .+= comparison_scalar_graph(Z,opts["threshold"])
                    warm_rho = isfinite(info.ρ) && info.ρ > 0 ? info.ρ : 1.0
                else
                    valid[index] = false
                    fill!(Z,0); fill!(U,0); warm_rho=1.0
                end
            catch err
                err isa InterruptException && rethrow()
                valid[index] = false
                comparison_record_failure(failures_io,
                    "model=$model replicate=$replicate seed=$seed method=$method time=$time index=$index",err,catch_backtrace())
                comparison_diagnostic(diagnostics_io,model,replicate,seed,method,time,index,
                    fraction,lambda,bound,false,false,NaN,NaN,NaN,0,"solver_error",sprint(showerror,err))
                fill!(Z,0); fill!(U,0); warm_rho=1.0
            end
        end
    end
    counts,valid
end

function comparison_run(opts,arguments)
    out = opts["out"]
    names = ("metrics.csv","solver_diagnostics.csv","edges.csv","cv_losses.csv","cv_folds.csv",
             "auc.csv","reference_endpoints.csv","settings.toml","failures.log")
    any(name -> ispath(joinpath(out,name)),names) &&
        throw(ArgumentError("output files already exist in $out; use a new --out"))
    mkpath(out)
    settings = copy(opts)
    settings["arguments"] = collect(arguments)
    settings["julia_version"] = string(VERSION)
    settings["package_version"] = string(Base.pkgversion(FunctionalGraphicalModels))
    settings["covariance_selection_version"] = string(Base.pkgversion(CovarianceSelection))
    settings["source_revision"] = get(ENV,"FGM_SOURCE_REV","not_recorded")
    settings["paper_source"] = "https://personal.lse.ac.uk/qiaox/qiao.links/research_papers/2017_FGM.pdf"
    settings["methods"] = collect(COMPARISON_METHODS)
    settings["generating_components"] = 5
    settings["noise_sd"] = 0.5
    settings["spline_degree"] = 3
    settings["spline_penalty"] = 0.0
    settings["cv"] = "random time-point folds; L,M jointly selected for fglasso; separate M=1 PCA selection"
    settings["tuning_seed_rule"] = "dataset_seed + 10000000; same initial RNG for fglasso and PCA partitions"
    settings["scalar_timepoint_input"] = "raw noisy observations; centered; sample covariance denominator n"
    settings["scalar_diagonals"] = "unpenalized"
    settings["scalar_timepoint_lambda"] = "one common lambda for all T timepoints; bound=1.05*max_t,max_j!=k |S_t[j,k]|"
    settings["penalty_grid"] = "common log fractions from 1 to lambda-min-ratio, times each method bound"
    settings["fgm_rho_units"] = "relative ADMM penalty; solver multiplies rho by squared average covariance variance; changes iterations, not the convex objective"
    settings["ALL"] = "edge selected in every timepoint network"
    settings["NEVER"] = "edge selected in at least one timepoint network"
    settings["HALF"] = "edge selected in strictly more than half of timepoint networks"
    settings["auc_rule"] = "trapezoids sorted by FPR; duplicate FPR retains maximum TPR; full_auc requires all fits valid and fitted empty+complete graphs"
    settings["reference_endpoints"] = "stored separately; never inserted into estimated path or AUC"
    settings["random_diagonal_rule"] = "0.5 * maximum_degree + 0.1; paper leaves diagonal choice unspecified"
    settings["reproduction_status"] = "paper-procedure comparison driver; bounded pilot defaults; no claimed agreement with published results"
    settings["comparison_limit"] = "coarse path areas do not establish method rankings; penalty-grid refinement and resolution of failed fits are required"
    open(joinpath(out,"settings.toml"),"w") do io
        TOML.print(io,settings; sorted=true)
    end
    ios = Dict(name=>open(joinpath(out,name),"w") for name in names if name != "settings.toml")
    metrics_io,diagnostics_io,edges_io = ios["metrics.csv"],ios["solver_diagnostics.csv"],ios["edges.csv"]
    cv_io,auc_io,endpoints_io,failures_io = ios["cv_losses.csv"],ios["auc.csv"],ios["reference_endpoints.csv"],ios["failures.log"]
    folds_io = ios["cv_folds.csv"]
    failed_points = 0
    try
        comparison_csvrow(metrics_io,"model","replicate","seed","method","lambda_index","lambda_fraction",
            "lambda","lambda_max_with_margin","valid","estimated_edges","tp","fp","tn","fn","tpr","fpr",
            "precision","recall","f1","accuracy","status","message")
        comparison_csvrow(diagnostics_io,"model","replicate","seed","method","time_index","lambda_index",
            "lambda_fraction","lambda","lambda_max_with_margin","converged","positive_definite",
            "kkt_residual","primal_residual","dual_residual","iterations","status","message")
        comparison_csvrow(edges_io,"model","replicate","seed","method","lambda_index","node1","node2","status")
        comparison_csvrow(cv_io,"model","replicate","seed","method","nbasis","ncomponents","mean_cv_loss","selected","status","reason")
        comparison_csvrow(folds_io,"model","replicate","seed","method","time_index","fold")
        comparison_csvrow(auc_io,"model","replicate","seed","method","valid_points","requested_points","fpr_min","fpr_max",
            "fitted_empty_observed","fitted_complete_observed","observed_span_area","full_fitted_auc","status")
        comparison_csvrow(endpoints_io,"model","replicate","seed","graph","status","tp","fp","tn","fn","tpr","fpr")
        fractions = exp.(range(0,log(opts["lambda-min-ratio"]); length=opts["nlambda"]))
        models = opts["model"] == "all" ? (:ar2,:alternating,:random) : (Symbol(opts["model"]),)
        for model in models,replicate in 1:opts["replicates"]
            seed = opts["seed"]+replicate-1
            sim = simulate_qiao2019(MersenneTwister(seed);n=opts["n"],p=opts["p"],
                ncomponents=5,ntime=opts["ntime"],model,noise_sd=0.5)
            p = opts["p"]
            complete = trues(p,p)
            for j in 1:p
                complete[j,j] = false
            end
            for (label,graph) in (("empty",falses(p,p)),("complete",complete))
                mt = graph_metrics(graph,sim.adjacency)
                comparison_csvrow(endpoints_io,model,replicate,seed,label,"reference_only",
                    mt.tp,mt.fp,mt.tn,mt.fn,mt.tpr,mt.fpr)
            end
            # Reset the RNG so both CV calls use the same initial time partition.
            tuning = select_fgm_components(sim.data;nbasis_grid=opts["nbasis-grid"],
                ncomponents_grid=opts["components-grid"],folds=opts["folds"],rng=MersenneTwister(seed+10_000_000))
            pca_tuning = select_fgm_components(sim.data;nbasis_grid=opts["nbasis-grid"],
                ncomponents_grid=[1],folds=opts["folds"],rng=MersenneTwister(seed+10_000_000))
            comparison_cv_rows(cv_io,tuning,model,replicate,seed,"fglasso")
            comparison_cv_rows(cv_io,pca_tuning,model,replicate,seed,"PCA")
            comparison_cv_fold_rows(folds_io,tuning,model,replicate,seed,"fglasso")
            comparison_cv_fold_rows(folds_io,pca_tuning,model,replicate,seed,"PCA")
            println("model=$model replicate=$replicate CV selected fglasso L=$(tuning.nbasis), M=$(tuning.ncomponents); PCA L=$(pca_tuning.nbasis)")
            flush(stdout)
            points = Dict(method=>NamedTuple[] for method in COMPARISON_METHODS)
            S = score_covariance(tuning.fpca.scores)
            bound = 1.05 * max(comparison_block_bound(S,tuning.ncomponents),eps(Float64))
            for (index,fraction) in pairs(fractions)
                lambda = fraction*bound
                try
                    fit = functional_glasso(S,tuning.ncomponents,lambda;
                        rho=opts["fgm-rho"],maxiter=opts["maxiter"],threshold=opts["threshold"])
                    d = fit.diagnostics
                    pd = isposdef(Symmetric(fit.precision))
                    valid = d.converged && pd
                    status = valid ? "ok" : "nonconverged"
                    comparison_diagnostic(diagnostics_io,model,replicate,seed,"fglasso",0,index,
                        fraction,lambda,bound,d.converged,pd,d.kkt_residual,d.primal_residual,
                        d.dual_residual,d.iterations,status,"")
                    push!(points["fglasso"],comparison_metrics(metrics_io,edges_io,fit.adjacency,
                        sim.adjacency,model,replicate,seed,"fglasso",index,fraction,lambda,bound,
                        valid,status,valid ? "" : "Graph metrics withheld."))
                    valid || (failed_points += 1)
                catch err
                    err isa InterruptException && rethrow()
                    failed_points += 1
                    comparison_record_failure(failures_io,"model=$model replicate=$replicate seed=$seed fglasso index=$index",err,catch_backtrace())
                    comparison_diagnostic(diagnostics_io,model,replicate,seed,"fglasso",0,index,
                        fraction,lambda,bound,false,false,NaN,NaN,NaN,0,"solver_error",sprint(showerror,err))
                    push!(points["fglasso"],comparison_metrics(metrics_io,edges_io,nothing,
                        sim.adjacency,model,replicate,seed,"fglasso",index,fraction,lambda,bound,
                        false,"solver_error",sprint(showerror,err)))
                end
            end
            raw_covariances = [comparison_covariance(sim.data.values[:,:,time]) for time in 1:opts["ntime"]]
            scalar_bound = 1.05 * max(maximum(comparison_scalar_bound,raw_covariances),eps(Float64))
            counts,valid = comparison_scalar_path!(diagnostics_io,failures_io,raw_covariances,
                scalar_bound,fractions,model,replicate,seed,opts)
            for (index,fraction) in pairs(fractions),method in ("ALL","NEVER","HALF")
                adjacency = method == "ALL" ? counts[:,:,index] .== opts["ntime"] :
                            method == "NEVER" ? counts[:,:,index] .> 0 :
                            counts[:,:,index] .> opts["ntime"]/2
                status = valid[index] ? "ok" : "invalid_timepoint_fit"
                push!(points[method],comparison_metrics(metrics_io,edges_io,adjacency,
                    sim.adjacency,model,replicate,seed,method,index,fraction,
                    fraction*scalar_bound,scalar_bound,valid[index],status,
                    valid[index] ? "" : "At least one timepoint fit failed; aggregate graph metrics withheld."))
                valid[index] || (failed_points += 1)
            end
            pca_covariance = score_covariance(pca_tuning.fpca.scores)
            pca_bound = 1.05 * max(comparison_scalar_bound(pca_covariance),eps(Float64))
            pca_counts,pca_valid = comparison_scalar_path!(diagnostics_io,failures_io,[pca_covariance],
                pca_bound,fractions,model,replicate,seed,opts;method="PCA")
            for (index,fraction) in pairs(fractions)
                adjacency = pca_counts[:,:,index] .> 0
                status = pca_valid[index] ? "ok" : "invalid_scalar_fit"
                push!(points["PCA"],comparison_metrics(metrics_io,edges_io,adjacency,
                    sim.adjacency,model,replicate,seed,"PCA",index,fraction,fraction*pca_bound,
                    pca_bound,pca_valid[index],status,pca_valid[index] ? "" : "Graph metrics withheld."))
                pca_valid[index] || (failed_points += 1)
            end
            for method in COMPARISON_METHODS
                comparison_auc(auc_io,points[method],model,replicate,seed,method,opts["nlambda"])
            end
            println("model=$model replicate=$replicate completed; cumulative invalid method/path points=$failed_points")
            flush(stdout)
        end
    catch err
        err isa InterruptException || comparison_record_failure(failures_io,
            "unhandled preprocessing or comparison error",err,catch_backtrace())
        rethrow()
    finally
        foreach(close,values(ios))
    end
    println("Comparison results saved to $out; invalid method/path points=$failed_points")
    iszero(failed_points) ? 0 : 1
end

function comparison_main(args=ARGS)
    opts = comparison_options(args)
    if opts === nothing
        comparison_usage()
        return 0
    end
    comparison_run(opts,args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(comparison_main())
end
