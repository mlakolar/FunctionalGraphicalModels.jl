# One process per setting/replicate. Run ONLY inside a Slurm allocation.
using FunctionalGraphicalModels, LinearAlgebra, Random, TOML, Serialization, SHA
import FunctionalDataAnalysis as FDA
import KernelSmoothing
import CovarianceSelection
include("dfgm_common.jl")
include("dfgm_factorized.jl")

function option(name,default)
    i = findfirst(==("--"*name),ARGS)
    i===nothing ? default : ARGS[i+1]
end
const METHODS = (:aclime,:l2_aclime,:scaled_lasso,:DFGM_C,:DFGM_G,
    :DGM_C,:DGM_G,:GGL,:Naive_C,:Naive_G)
const THRESHOLD = 1e-6
csv(x) = "\""*replace(string(x),'"'=>"\"\"")*"\""
function row(io,args...)
    println(io,join(csv.(args),",")); flush(io)
end

function sample_covariances(data,grid)
    n,p = FDA.nobs(data),FDA.ncomponents(data)
    X = zeros(n,p,length(grid))
    for j in 1:p,i in 1:n
        ti,yi = FDA.component(data,j)[i]
        ti==grid || error("sample covariance requires common observations")
        X[i,j,:] .= yi
    end
    S = zeros(p,p,length(grid))
    for t in eachindex(grid)
        Y = X[:,:,t] .- sum(X[:,:,t];dims=1)./n
        S[:,:,t] .= transpose(Y)*Y/n
    end
    S
end

function reconstructed_covariances(data,grid;mode=:pilot,bandwidth=.25,work_points=31)
    n,p = FDA.nobs(data),FDA.ncomponents(data)
    workgrid = sort!(unique(vcat(collect(range(0,1;length=work_points)),grid)))
    X = zeros(n,p,length(grid))
    selections = NamedTuple[]
    cv_cache = Dict()
    for j in 1:p
        fd = FDA.component(data,j)
        selector = mode===:paper ? FDA.CurveCV([.1,.2,.3,.4];folds=n) : bandwidth
        # This is the literal whole-curve CV path. It can be expensive; pilot
        # mode fixes h and uses AIC for rank. Neither mode knows latent truth.
        fp = if all(fd[i][1]==fd[1][1] for i in 1:n)
            meanfit = FDA.fit(FDA.SmoothedMean,fd;grid=workgrid,bandwidth=selector,
                kernel=KernelSmoothing.GaussianKernel(),weighting=:subj)
            covariance = _dfgm_common_covariance(fd,meanfit,workgrid;bandwidth=selector,cache=cv_cache)
            FDA.fit(FDA.FPCA,fd;method=:pace,grid=workgrid,ncomponents=FDA.AIC(),
                maxK=5,mean=meanfit,covariance)
        else
            meanfit = FDA.fit(FDA.SmoothedMean,fd;grid=workgrid,bandwidth=selector,
                kernel=KernelSmoothing.GaussianKernel(),weighting=:subj)
            covariance = _dfgm_factorized_covariance(fd,meanfit,workgrid;bandwidth=selector)
            FDA.fit(FDA.FPCA,fd;method=:pace,grid=workgrid,ncomponents=FDA.AIC(),
                maxK=5,mean=meanfit,covariance)
        end
        basis = FDA.basismatrix(fp.basis,grid)
        X[:,j,:] .= fp.scores*transpose(basis)
        push!(selections,(;node=j,rank=length(fp.eigenvalues),settings=fp.settings))
    end
    S = zeros(p,p,length(grid))
    for t in eachindex(grid)
        S[:,:,t] .= transpose(X[:,:,t])*X[:,:,t]/n
    end
    S,selections
end

function fit_scalar_path(S,kind,penalty)
    p,_,q = size(S)
    P = zeros(p,p,q)
    for t in 1:q
        if penalty==0
            isposdef(Symmetric(S[:,:,t])) || error("zero-penalty covariance is not positive definite at time $t")
            P[:,:,t] .= inv(Symmetric(S[:,:,t]))
        elseif kind===:clime
            P[:,:,t] .= Matrix(CovarianceSelection.clime(S[:,:,t],penalty;ntasks=1))
        else
            f = functional_glasso(S[:,:,t],1,penalty;maxiter=10000)
            f.diagnostics.converged || error("glasso did not converge at time $t: $(f.diagnostics)")
            P[:,:,t] .= f.sparse_precision
        end
    end
    P
end

function metrics(P,truth)
    q = size(P,3)
    operator = sum(opnorm(P[:,:,t]-truth[:,:,t]) for t in 1:q)/q
    frobenius = sum(norm(P[:,:,t]-truth[:,:,t]) for t in 1:q)/q
    normalized_squared_frobenius = sum(sum(abs2,P[:,:,t]-truth[:,:,t])/size(P,1) for t in 1:q)/q
    operator,frobenius,normalized_squared_frobenius
end

function rocpoint(P,truth,t)
    p = size(P,1)
    tp,fp,pos,neg = 0,0,0,0
    for j in 1:p,k in j+1:p
        target = abs(truth[j,k,t])>THRESHOLD
        pred = abs(P[j,k,t])>THRESHOLD
        pos+=target; neg+=!target; tp+=target&&pred; fp+=!target&&pred
    end
    (neg==0 ? NaN : fp/neg,pos==0 ? NaN : tp/pos,pos,neg)
end

# Full AUROC requires fitted endpoints. Never append reference graphs to a path.
function observed_auc(points)
    usable = [(x,y) for (x,y) in points if isfinite(x)&&isfinite(y)]
    isempty(usable) && return (NaN,NaN,NaN,NaN)
    xs = sort!(unique(first.(usable)))
    ys = [maximum(y for (xx,y) in usable if xx==x) for x in xs]
    area = sum(((xs[i+1]-xs[i])*(ys[i+1]+ys[i])/2 for i in 1:length(xs)-1);init=0.)
    full = any(x==0&&y==0 for (x,y) in usable) && any(x==1&&y==1 for (x,y) in usable)
    (area,full ? area : NaN,first(xs),last(xs))
end

function worker(setting_index,replicate,out;mode=:pilot,n_override=0,p_override=0,
                integration_points=101,nlambda=7,work_points=31,evaluation_points=0,
                bias_constant=1.,bias_constants=nothing,lp_solver=:highs)
    mkpath(out)
    done = joinpath(out,"complete.toml")
    (isfile(done)||isfile(joinpath(out,"settings.toml"))) &&
        error("output already has results: $out; choose another attempt path")
    setting = dfgm2020_settings()[setting_index]
    seed = 2026102500+10000setting_index+replicate
    rng = MersenneTwister(seed)
    started = time()
    effective_n = n_override>0 ? n_override : setting.n
    effective_p = p_override>0 ? p_override : setting.p
    grid_override = evaluation_points>0 ? collect(range(0,1;length=evaluation_points)) : nothing
    sim = simulate_dfgm2020(rng;p=effective_p,n=effective_n,design=setting.design,
        schedule=setting.schedule,integration_points,evaluation_grid=grid_override)
    println("generated setting=$setting_index p=$effective_p n=$effective_n seconds=$(time()-started)");flush(stdout)
    # Common methods use actual common observations, irrespective of a reduced
    # evaluation grid. For paper settings evaluation is the observation grid.
    obsgrid = sim.data.components[1][1][1]
    if setting.design===:common && evaluation_points>0
        error("reduced evaluation grids are supported only for irregular pilots")
    end
    methods = setting.design===:common ? collect(METHODS) : collect(METHODS[1:5])
    nlambda>=2 || error("nlambda must be at least two")
    mode in (:pilot,:paper) || error("mode must be :pilot or :paper")
    penalties = vcat(0.,exp.(range(log(.02),log(1.2);length=nlambda-1)))
    # K2 is tuned on the same fixed menu for both proposals. It is not a theorem
    # constant in these exploratory paths. Fixed K1 controls the shared pilot.
    multipliers = exp.(range(log(.03),log(2.);length=nlambda))
    bias_menu = bias_constants===nothing ? [Float64(bias_constant)] : Float64.(collect(bias_constants))
    !isempty(bias_menu) && all(b->isfinite(b)&&b>0,bias_menu) && allunique(bias_menu) ||
        error("bias constants must be distinct, positive and finite")
    direct_candidates = [(;K2,bias_constant=b) for b in bias_menu for K2 in multipliers]
    hs = [.15,.25,.4]
    C = direct_covariance(sim.data,sim.grid;bandwidths=hs,
        design=setting.design===:common ? :stratified_synchronous : :independent,
        boundary=:local_linear)
    # Retain the original default separately. This empirical menu addresses
    # the observed sensitivity to sparse-visit covariance error. Paper-mode
    # selection uses the same observable pointwise BIC as the other paths;
    # it is distinct from the held-out protocol in the controlled study.
    scaled_values = Any[:quantile,.2,.35,.5,.75,
        sqrt(2log(effective_p)/effective_n)*(1+1/(C.min_visits*hs[2]))]
    push!(methods,:scaled_lasso_tuned)
    expected_candidates(method) = method===:scaled_lasso ? 1 :
        method===:scaled_lasso_tuned ? length(scaled_values) :
        method in (:aclime,:l2_aclime) ? length(direct_candidates) : nlambda
    println("direct covariance ready seconds=$(time()-started)");flush(stdout)
    settings = Dict("setting_index"=>setting_index,"replicate"=>replicate,"seed"=>seed,
        "design"=>string(setting.design),"schedule"=>string(setting.schedule),
        "p"=>effective_p,"n"=>effective_n,"mode"=>string(mode),
        "integration_points"=>integration_points,"work_points"=>work_points,
        "evaluation_grid"=>sim.grid,"bandwidths"=>hs,"K1"=>2.,"K2_menu"=>multipliers,
        "scalar_penalties"=>penalties,"bias_constant"=>bias_constant,"lp_solver"=>string(lp_solver),"threshold"=>THRESHOLD,
        "bias_constant_menu"=>bias_menu,
        "direct_selection"=>"pointwise BIC across bias/K2 menu; empirical when envelope not externally certified",
        "expected_candidates"=>Dict(string(m)=>expected_candidates(m) for m in methods),
        "noise_variance"=>.5,"julia_version"=>string(VERSION),
        "scaled_lasso_penalty_menu"=>string.(scaled_values),
        "scaled_lasso_tuned_selection"=>"pointwise BIC; empirical menu, no iid-theory transfer",
        "baseline_bandwidth_selection"=>mode===:paper ? "whole_curve_CV" : "fixed_0.25",
        "baseline_rank_selection"=>"FDA_AIC_maxK5",
        "truth_selection"=>"no_truth_used_for_tuning",
        "source_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f))))
            for f in ("Project.toml","src/FunctionalGraphicalModels.jl","src/direct.jl",
                "src/dfgm_simulation.jl","experiments/dfgm2020.jl",
                "experiments/dfgm_common.jl","experiments/dfgm_factorized.jl")))
    manifest=joinpath(@__DIR__,"..","Manifest.toml")
    settings["dependency_manifest_sha256"]=isfile(manifest) ? bytes2hex(sha256(read(manifest))) : "not_saved"
    open(joinpath(out,"settings.toml"),"w") do io; TOML.print(io,settings); end
    # Keep paired reporting inputs; no raw simulated data are needed to recreate.
    serialize(joinpath(out,"targets.bin"),(;grid=sim.grid,actual=sim.precision,
        nominal=sim.nominal_precision,graph=sim.graph))
    mean_mismatch = metrics(sim.precision,sim.nominal_precision)
    open(joinpath(out,"target_discrepancy.toml"),"w") do io
        TOML.print(io,Dict("operator"=>mean_mismatch[1],"frobenius"=>mean_mismatch[2],
            "normalized_squared_frobenius"=>mean_mismatch[3]))
    end
    paths = Dict{Symbol,Vector{Tuple{Int,Array{Float64,3}}}}(m=>[] for m in methods)
    failures = 0
    Srecon = nothing
    reconstruction_error = ""
    try
        Srecon,selections = reconstructed_covariances(sim.data,sim.grid;mode,work_points)
        serialize(joinpath(out,"fpca_selections.bin"),selections)
    catch e
        reconstruction_error = sprint(showerror,e,catch_backtrace())
    end
    println("FPCA ready success=$(Srecon!==nothing) seconds=$(time()-started)");flush(stdout)
    Sraw = setting.design===:common ? sample_covariances(sim.data,obsgrid) : nothing
    Sdgm = nothing
    if Sraw!==nothing
        T = length(obsgrid)
        h = max((log(effective_p)/(effective_n*T))^(1/3),T^(-4/5))
        Sdgm = similar(Sraw)
        for (t,u) in enumerate(obsgrid)
            w = exp.(-0.5 .* ((obsgrid.-u)./h).^2); w ./= sum(w)
            Sdgm[:,:,t] .= sum(w[r].*Sraw[:,:,r] for r in eachindex(obsgrid))
        end
    end
    # Expensive reconstruction inputs are reusable for later precision-only
    # sensitivity analyses. Their source, subjects and seeds are in settings.
    serialize(joinpath(out,"covariances.bin"),(;grid=sim.grid,direct=C,
        reconstructed=Srecon,raw=Sraw,dynamic=Sdgm))
    open(joinpath(out,"fits.csv"),"w") do io
        row(io,"method","candidate","penalty","bias_constant","status","seconds","pd_slices","max_constraint_violation","psd_displacement","error")
        for method in methods
            println("starting $method seconds=$(time()-started)");flush(stdout)
            for candidate in 1:expected_candidates(method)
                tick = time()
                scaled = method in (:scaled_lasso,:scaled_lasso_tuned)
                direct_parameter = method in (:aclime,:l2_aclime) ? direct_candidates[candidate] : nothing
                bvalue = direct_parameter===nothing ? bias_constant : direct_parameter.bias_constant
                value = direct_parameter!==nothing ? direct_parameter.K2 :
                    scaled ? scaled_values[candidate] : penalties[candidate]
                try
                    violation,displacement = 0.,0.
                    P = if method in (:aclime,:l2_aclime,:scaled_lasso,:scaled_lasso_tuned)
                        f = fit_direct(C;method=scaled ? :scaled_lasso : method,
                            pilot_index=2,bias_constant=bvalue,K2=scaled ? 1. : value,
                            lambda0=scaled ? value : :quantile,noise_variance=.5,lp_solver)
                        violation = maximum(d.max_constraint_violation for d in f.diagnostics)
                        displacement = maximum(d.covariance_psd_displacement for d in f.diagnostics)
                        f.precisions
                    elseif method in (:DFGM_C,:DFGM_G)
                        Srecon===nothing && error("FPCA failed: "*reconstruction_error)
                        fit_scalar_path(Srecon,method===:DFGM_C ? :clime : :glasso,value)
                    elseif method===:GGL
                        if value==0
                            fit_scalar_path(Sraw,:glasso,0.)
                        else
                            f = joint_graphical_lasso(Sraw;lambda=value,alpha=0.,maxiter=10000)
                            f.converged || error("group glasso failed: residual $(f.primal_residual), $(f.dual_residual)")
                            f.precision
                        end
                    else
                        S = method in (:DGM_C,:DGM_G) ? Sdgm : Sraw
                        kind = method in (:DGM_C,:Naive_C) ? :clime : :glasso
                        fit_scalar_path(S,kind,value)
                    end
                    pd = count(t->isposdef(Symmetric(P[:,:,t])),axes(P,3))
                    row(io,method,candidate,value,bvalue,"ok",time()-tick,pd,violation,displacement,"")
                    push!(paths[method],(candidate,P))
                catch e
                    failures += 1
                    row(io,method,candidate,value,bvalue,"failed",time()-tick,0,NaN,NaN,sprint(showerror,e))
                end
            end
        end
    end
    # BIC is computed per time using each method's observed covariance. Direct
    # estimates can be indefinite: those BIC candidates receive Inf, never a
    # hidden SPD precision repair. Target losses are computed after selection.
    selected_bias = Dict{Tuple{Symbol,Int},Float64}()
    open(joinpath(out,"metrics.csv"),"w") do io
        row(io,"method","target","time_index","candidate","bic","operator","frobenius","normalized_squared_frobenius","fpr","tpr","edges","nonedges")
        for method in methods
            path = paths[method]
            isempty(path) && continue
            observed = method in (:DFGM_C,:DFGM_G) ? Srecon :
                method in (:DGM_C,:DGM_G) ? Sdgm :
                method in (:GGL,:Naive_C,:Naive_G) ? Sraw : C.matrices[:,:,:,2]
            for t in eachindex(sim.grid)
                bic = [isposdef(Symmetric(P[:,:,t])) ? effective_n*(dot(P[:,:,t],observed[:,:,t])-
                    logdet(Symmetric(P[:,:,t])))+log(effective_n)*count(>(THRESHOLD),abs.(P[:,:,t])) : Inf
                    for (_,P) in path]
                best = argmin(bic)
                if method in (:aclime,:l2_aclime) && isfinite(bic[best])
                    selected_bias[method,t] = direct_candidates[path[best][1]].bias_constant
                end
                for (k,(candidate,P)) in enumerate(path), (target,truth) in
                        ((:actual,sim.precision),(:nominal,sim.nominal_precision))
                    loss = metrics(P[:,:,t:t],truth[:,:,t:t])
                    fp,tp,pos,neg = rocpoint(P,truth,t)
                    row(io,method,target,t,candidate,bic[k],loss...,fp,tp,pos,neg)
                end
                # Save selected index even if every candidate has invalid BIC.
                open(joinpath(out,"selected.csv"),"a") do sel
                    row(sel,method,t,isfinite(bic[best]) ? path[best][1] : 0,
                        isfinite(bic[best]) ? "selected" : "no_positive_definite_candidate")
                end
            end
        end
    end
    # Keep precision paths for numerical re-evaluation and paired comparisons.
    serialize(joinpath(out,"paths.bin"),paths)
    open(joinpath(out,"auc.csv"),"w") do io
        row(io,"method","target","time_index","observed_area","full_auc","minimum_fpr","maximum_fpr","successful_candidates","expected_candidates","selected_bias_constant")
        for method in methods,(target,truth) in ((:actual,sim.precision),(:nominal,sim.nominal_precision)),t in eachindex(sim.grid)
            path = paths[method]
            expected = expected_candidates(method)
            chosen_bias = NaN
            if method in (:aclime,:l2_aclime) && length(bias_menu)>1
                # Fix the bias constant using observable BIC before computing
                # the K2 ROC. Pooling all bias paths would use target TPR to
                # form an optimistic envelope across an extra tuning dimension.
                chosen_bias = get(selected_bias,(method,t),NaN)
                path = [(k,P) for (k,P) in path if direct_candidates[k].bias_constant==chosen_bias]
                expected = nlambda
            end
            points = [rocpoint(P,truth,t)[1:2] for (_,P) in path]
            area,full,lo,hi = observed_auc(points)
            length(path)==expected || (full=NaN)
            row(io,method,target,t,area,full,lo,hi,length(path),expected,chosen_bias)
        end
    end
    open(done,"w") do io
        TOML.print(io,Dict("status"=>failures==0 ? "completed" : "completed_with_failures",
            "fit_failures"=>failures,"elapsed_seconds"=>time()-started,
            "paper_reproduction"=>false,"mode"=>string(mode)))
    end
    println("setting=$setting_index replicate=$replicate failures=$failures seconds=$(time()-started)")
    failures==0
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run scientific work inside Slurm")
    index = parse(Int,option("setting","1"))
    replicate = parse(Int,option("replicate","1"))
    ok = worker(index,replicate,option("out","results/dfgm2020");
        mode=Symbol(option("mode","pilot")),n_override=parse(Int,option("n","0")),
        p_override=parse(Int,option("p","0")),integration_points=parse(Int,option("integration-points","101")),
        nlambda=parse(Int,option("nlambda","7")),work_points=parse(Int,option("work-points","31")),
        evaluation_points=parse(Int,option("evaluation-points","0")),
        bias_constant=parse(Float64,option("bias-constant","1")),
        bias_constants=let value=option("bias-constants",""); isempty(value) ? nothing : parse.(Float64,split(value,',')); end,
        lp_solver=Symbol(option("lp-solver","highs")))
    exit(ok ? 0 : 1)
end
