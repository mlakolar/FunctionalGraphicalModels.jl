# Program(P) extension; existing baseline fits are read, never recomputed.
using FunctionalGraphicalModels, LinearAlgebra, Random, Serialization, TOML, SHA, Dates, Sockets
const FGM=FunctionalGraphicalModels
const ROOT=get(ENV,"DFGM_ROOT","/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct")
const SOURCE=get(ENV,"DFGM_SPECTRAL_SOURCE",dirname(@__DIR__))
const RUNTIME=joinpath(ROOT,"common-fast-v2-source")
Base.include(FGM,joinpath(SOURCE,"src/spectral_direct.jl"))
include(joinpath(RUNTIME,"experiments/direct_stress.jl"))
const SCORE_MENU=[.01,.03,.1,.3,1.,2.]

function observation_hash(data)
    io=IOBuffer()
    for fd in data.components,i in 1:FDA.nobs(data)
        ts,ys=fd[i];write(io,Int64(length(ts)));write(io,ts);write(io,ys)
    end
    bytes2hex(sha256(take!(io)))
end
function subset_covariance(C,ix)
    FGM.DirectCovariance(C.grid[ix],C.bandwidths,C.matrices[:,:,ix,:],C.n,C.min_visits,C.pair_subjects[:,ix,:],C.settings)
end
function source_hashes()
    Dict(f=>bytes2hex(sha256(read(joinpath(SOURCE,f)))) for f in
        ("src/spectral_direct.jl","experiments/spectral_comparison.jl","docs/SPECTRAL_L2_ACLIME.md"))
end
function spectral_dataset(kind,index,rep,out;point_index=0,ntasks=1,resume_folder="")
    get(ENV,"SLURM_JOB_ID","")!="" || error("scientific work requires Slurm")
    1<=ntasks<=Threads.nthreads() || error("invalid time-task count")
    !ispath(out) || error("existing output: $out; reconcile before any retry")
    mkpath(out);started=time();seed=2026102500+10000index+rep
    regeneration=time();cache=nothing;validation=nothing;identity=Dict{String,Any}()
    if kind=="paper"
        1<=rep<=20 || error("paper spectral scope is fixed first20")
        folder=realpath(joinpath(ROOT,"outputs/paper-v4/setting-$index-rep-$rep"))
        prior=TOML.parsefile(joinpath(folder,"settings.toml"));targets=deserialize(joinpath(folder,"targets.bin"))
        prior["seed"]==seed && prior["noise_variance"]==.5 || error("baseline identity differs")
        for f in ("src/direct.jl","src/dfgm_simulation.jl")
            prior["source_sha256"][f]==bytes2hex(sha256(read(joinpath(RUNTIME,f)))) || error("baseline source differs: $f")
        end
        setting=dfgm2020_settings()[index]
        sim=simulate_dfgm2020(MersenneTwister(seed);p=prior["p"],n=prior["n"],design=setting.design,
            schedule=setting.schedule,integration_points=prior["integration_points"],noise_variance=.5)
        cache=deserialize(joinpath(folder,"covariances.bin")).direct
        sim.grid==cache.grid==targets.grid || error("regenerated grid differs")
        target_difference=maximum(abs.(sim.precision.-targets.actual))
        target_difference<=1e-10*(1+maximum(abs,targets.actual)) || error("regenerated target differs")
        data=sim.data;grid=sim.grid;actual=targets.actual;nominal=targets.nominal
        design=setting.design===:common ? :stratified_synchronous : :independent
        boundary=:local_linear;noise_variance=.5;L=1.;kappa_bound=10.
        identity["baseline_folder"]=realpath(folder);identity["target_max_difference"]=target_difference
        identity["baseline_settings_sha256"]=bytes2hex(sha256(read(joinpath(folder,"settings.toml"))))
        identity["baseline_covariance_sha256"]=bytes2hex(sha256(read(joinpath(folder,"covariances.bin"))))
        identity["calibration_status"]="L and kappa are prespecified empirical sensitivity inputs, not certified class bounds"
    elseif kind=="control"
        1<=rep<=100 || error("control scope is recorded100")
        configurations=[(;degree,m,curvature,noise_variance) for degree in (2,10) for m in (8,30) for curvature in (0.,4.) for noise_variance in (.5,2.)]
        config=configurations[index];sim=stress_data(MersenneTwister(seed);config...)
        data=sim.data[1:160];validation=sim.data[161:200];grid=[.5]
        actual=reshape(sim.precision,50,50,1);nominal=actual;noise_variance=config.noise_variance
        c=config.curvature;L=max(.01,.32c*(1+c/4),.16c^2);kappa_bound=2.
        design=:independent;boundary=:error
        identity["calibration_status"]="known control curvature bound; fixed class condition bound2; nonunit-noise extension"
        identity["train_subjects"]=160;identity["validation_subjects"]=40
        identity["generating_function_sha256"]=bytes2hex(sha256(read(joinpath(RUNTIME,"experiments/direct_stress.jl"))))
        identity["baseline_shared_covariance_available"]=false
    else
        error("unknown dataset family")
    end
    regeneration_seconds=time()-regeneration
    # Compare original .25 occupancy statistics before computing extra widths.
    shared_tick=time();shared_difference=NaN
    if cache!==nothing
        check_index=cld(length(grid),2)
        check=direct_covariance(data,grid[check_index:check_index];bandwidths=[.25],design,boundary)
        k=findfirst(==(.25),cache.bandwidths)
        shared_difference=maximum(abs.(check.matrices[:,:,1,1].-cache.matrices[:,:,check_index,k]))
        shared_difference<=1e-10*(1+maximum(abs,cache.matrices[:,:,check_index,k])) || error("regenerated observations fail shared covariance check")
        check.pair_subjects[:,1,1]==cache.pair_subjects[:,check_index,k] || error("shared occupancy counts differ")
        identity["shared_covariance_checked_time_indices"]=[check_index]
    end
    identity["shared_covariance_max_difference"]=shared_difference
    identity["shared_covariance_check_seconds"]=time()-shared_tick
    identity["observation_sha256"]=observation_hash(data)
    if point_index>0
        kind=="paper" || error("point pilot only supports paper")
        ix=point_index:point_index;grid=grid[ix];actual=actual[:,:,ix];nominal=nominal[:,:,ix];cache=subset_covariance(cache,ix)
    end
    calibration_tick=time();n=FDA.nobs(data);p=FDA.ncomponents(data)
    m=minimum(length(fd[i][1]) for fd in data.components for i in 1:n)
    cal=FGM.spectral_calibration(n,p,m;L,kappa_bound,h_max=.25,C_K=1.,Aconf=1.)
    calibration_seconds=time()-calibration_tick
    covariance_tick=time();built=FGM.spectral_covariance(data,grid,cal;cache,cache_mean_matches=cache!==nothing,design,boundary)
    C=built.covariance;covariance_seconds=time()-covariance_tick
    validation_tick=time();Cv=validation===nothing ? nothing : direct_covariance(validation,[.5];bandwidths=[.2])
    validation_covariance_seconds=time()-validation_tick
    serialize(joinpath(out,"covariances.bin"),(;direct=C,validation=Cv))
    serialize(joinpath(out,"targets.bin"),(;actual,nominal,grid))
    settings=Dict("kind"=>kind,"setting_index"=>index,"replicate"=>rep,"seed"=>seed,"n"=>n,"p"=>p,
        "minimum_visits"=>m,"sampling_schedule"=>kind=="paper" ? string(dfgm2020_settings()[index].schedule) : "fixed_m",
        "evaluation_grid"=>grid,"bandwidths"=>cal.bandwidths,"diagonal_bandwidth"=>cal.h_d,
        "grid_mode"=>"exact_theory_grid","L"=>L,"C_K"=>1.,"operator_bias_coefficient"=>L,"kappa_bound"=>kappa_bound,
        "Aconf"=>1.,"ell"=>cal.ell,"score_menu"=>SCORE_MENU,"fixed_Kscore"=>2.,"noise_variance"=>noise_variance,
        "domain_length"=>1.,"centering"=>"known_zero","boundary"=>string(boundary),"sampling_design"=>string(design),
        "primary_target"=>"actual","secondary_target"=>"nominal","point_pilot"=>point_index>0,
        "selection"=>kind=="paper" ? "pointwise observed-data BIC, PD candidates only" : "heldout covariance precision residual, finite optimal candidates",
        "threads"=>Threads.nthreads(),"time_tasks"=>ntasks,"solver_threads"=>1,"hostname"=>gethostname(),"cpu_model"=>Sys.cpu_info()[1].model,
        "slurm_job_id"=>ENV["SLURM_JOB_ID"],"julia_version"=>string(VERSION),"source_sha256"=>source_hashes(),
        "runtime_project"=>RUNTIME,"runtime_manifest_sha256"=>bytes2hex(sha256(read(joinpath(RUNTIME,"Manifest.toml")))),
        "path_reuse"=>true,"warmup"=>"none; first-call compilation included in stage timings",
        "timing_schema"=>"private_point_wall_clocks_v2",
        "precision_pilot"=>false,"psd_repair"=>false,"threshold"=>THRESHOLD,
        "candidate_count"=>6,"fixed_candidate_index"=>6,"identity"=>identity,
        "reused_covariance_widths"=>built.reused_widths,"additional_covariance_widths"=>built.computed_widths,
        "process_id"=>getpid(),"case_sequence"=>get(ENV,"DFGM_CASE_SEQUENCE","1"),
        "regeneration_seconds"=>regeneration_seconds,"calibration_seconds"=>calibration_seconds,"covariance_seconds"=>covariance_seconds,
        "validation_covariance_seconds"=>validation_covariance_seconds,"resume_folder"=>resume_folder)
    if !isempty(resume_folder)
        isfile(joinpath(resume_folder,"complete.toml")) && error("completed dataset must not be refitted")
        prior=TOML.parsefile(joinpath(resume_folder,"settings.toml"))
        old_job=prior["slurm_job_id"]
        occursin(r"^[0-9]+$",old_job) || error("invalid original Slurm identity")
        isempty(readchomp(`squeue -h -j $old_job -o %i`)) || error("original allocation still active")
        state=readchomp(`sacct -X -n -j $old_job --format=State --parsable2`)
        any(startswith(state,s) for s in ("COMPLETED","FAILED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","CANCELLED","PREEMPTED")) || error("original terminal state not established")
        for key in ("kind","setting_index","replicate","seed","n","p","evaluation_grid","bandwidths","L","C_K","kappa_bound","Aconf","score_menu","noise_variance","source_sha256")
            prior[key]==settings[key] || error("resume identity mismatch: $key")
        end
        prior["identity"]["observation_sha256"]==identity["observation_sha256"] || error("resume observations differ")
        old=deserialize(joinpath(resume_folder,"covariances.bin")).direct
        old.matrices==C.matrices && old.pair_subjects==C.pair_subjects || error("resume covariance differs")
    end
    open(joinpath(out,"settings.toml"),"w") do io;TOML.print(io,settings);end
    println("spectral $kind/$index/$rep grid=$(length(cal.bandwidths)) covariance_seconds=$covariance_seconds");flush(stdout)
    q=length(grid);paths=[fill(NaN,p,p,q) for _ in SCORE_MENU];columns=[fill(NaN,p,p,q) for _ in SCORE_MENU]
    results=Vector{Any}(undef,q)
    function at_time(t)
        if !isempty(resume_folder) && isfile(joinpath(resume_folder,"checkpoint-time-$t.bin"))
            cp=deserialize(joinpath(resume_folder,"checkpoint-time-$t.bin"))
            cp.source_sha256==source_hashes() && cp.observation_sha256==identity["observation_sha256"] && cp.time_index==t || error("checkpoint identity differs")
            for k in eachindex(SCORE_MENU)
                paths[k][:,:,t].=cp.precisions[k];columns[k][:,:,t].=cp.columns[k]
            end
            results[t]=merge(cp.result,(;restored=true))
            Base.cp(joinpath(resume_folder,"checkpoint-time-$t.bin"),joinpath(out,"checkpoint-time-$t.bin"))
            println("restored completed spectral time=$t/$q");flush(stdout)
            return
        end
        diagonal_tick=time();d=max.(diag(C.matrices[:,:,t,cal.diagonal_index]),n^-2);diagonal_seconds=time()-diagonal_tick
        point_path_start=time();rs=FGM.spectral_l2_aclime_path(C.matrices[:,:,t,:],SCORE_MENU;d,
            u=cal.u,v=cal.v,bias=cal.bias,kappa_bound,noise_variance)
        path_seconds=time()-point_path_start;stages=NamedTuple[]
        for (k,r) in enumerate(rs)
            columns[k][:,:,t].=r.columns;symmetry_seconds=0.;cap_seconds=0.;cap=NaN;status=r.status;err=""
            if all(isfinite,r.columns)
                try
                    tick=time();S=interval_symmetrize(r.columns);symmetry_seconds=time()-tick
                    tick=time();norm1=opnorm(S,1);cap=norm1==0 ? 1. : min(1.,n/norm1);paths[k][:,:,t].=cap.*S;cap_seconds=time()-tick
                catch e;status=:numerical_failure;err=sprint(showerror,e);end
            end
            push!(stages,(;status,error=err,cap,symmetry_seconds,cap_seconds,column_records=r.column_records))
        end
        # Bounded independent replay witnesses; all other columns retain checked residuals.
        witness=t==1 ? [(;Kscore=r.Kscore,columns=r.columns,budgets=r.budgets,bias_vectors=r.bias_vectors,d,B=r.B) for r in rs] : nothing
        results[t]=(;diagonal_seconds,path_seconds,stages,witness,restored=false)
        checkpoint=(;source_sha256=source_hashes(),observation_sha256=identity["observation_sha256"],time_index=t,
            precisions=[copy(P[:,:,t]) for P in paths],columns=[copy(W[:,:,t]) for W in columns],result=results[t])
        temporary=joinpath(out,"checkpoint-time-$t.tmp")
        serialize(temporary,checkpoint);mv(temporary,joinpath(out,"checkpoint-time-$t.bin"))
        println("spectral time=$t/$q seconds=$path_seconds statuses=$(getfield.(stages,:status))");flush(stdout)
    end
    precision_wall_start=time()
    if ntasks==1
        for t in 1:q;at_time(t);end
    else
        @sync for task in 1:min(ntasks,q)
            Threads.@spawn for t in task:ntasks:q;at_time(t);end
        end
    end
    fit_wall_seconds=time()-precision_wall_start
    serialize(joinpath(out,"paths.bin"),(;paths,columns,score_menu=SCORE_MENU))
    serialize(joinpath(out,"audit-witnesses.bin"),results[1].witness)
    serialize(joinpath(out,"diagnostics.bin"),results)
    failures=0
    open(joinpath(out,"fits.csv"),"w") do io
        row(io,"candidate","Kscore","time_index","status","positive_definite","precision_seconds","model_setup_seconds","solve_seconds","residual_seconds","symmetry_seconds","cap_seconds","cap","error")
        for t in 1:q,k in eachindex(SCORE_MENU)
            s=results[t].stages[k];failures+=s.status!==:optimal
            P=paths[k][:,:,t];pd=all(isfinite,P)&&isposdef(Symmetric(P));records=s.column_records
            row(io,k,SCORE_MENU[k],t,s.status,pd,sum(x.seconds for x in records),sum(x.model_setup_seconds for x in records),
                sum(x.solve_seconds for x in records),sum(x.residual_seconds for x in records),s.symmetry_seconds,s.cap_seconds,s.cap,s.error)
        end
    end
    open(joinpath(out,"column-status.csv"),"w") do io
        names=propertynames(results[1].stages[1].column_records[1]);row(io,"candidate","time_index",names...)
        for t in 1:q,k in eachindex(SCORE_MENU),record in results[t].stages[k].column_records;row(io,k,t,values(record)...);end
    end
    selected=Dict{Int,Int}();score_tick=time()
    open(joinpath(out,"metrics.csv"),"w") do io
        row(io,"method","target","time_index","candidate","bic","operator","frobenius","normalized_squared_frobenius","fpr","tpr","positive_edges","negative_edges","status","pd")
        for t in 1:q
            scores=fill(Inf,6)
            for k in 1:6
                P=paths[k][:,:,t];s=results[t].stages[k];finite=all(isfinite,P);pd=finite&&isposdef(Symmetric(P))
                if s.status===:optimal && finite
                    if Cv!==nothing
                        scores[k]=norm(Cv.matrices[:,:,1,1]*P-I)
                    elseif pd
                        # Same fixed .25 observed covariance as existing direct
                        # methods; h_d is only for diagonal calibration.
                        observed=C.matrices[:,:,t,findfirst(==(.25),C.bandwidths)]
                        scores[k]=n*(dot(P,observed)-logdet(Symmetric(P)))+log(n)*sum(abs.(P).>THRESHOLD)
                    end
                end
                for (target,truth) in ((:actual,actual),(:nominal,nominal))
                    loss=finite ? metrics(reshape(P,p,p,1),truth[:,:,t:t]) : (NaN,NaN,NaN)
                    graph=finite ? rocpoint(reshape(P,p,p,1),truth[:,:,t:t],1) : (NaN,NaN,0,0)
                    row(io,"spectral_l2_aclime_path",target,t,k,scores[k],loss...,graph...,s.status,pd)
                end
            end
            k=argmin(scores);selected[t]=isfinite(scores[k]) ? k : 0
        end
    end
    open(joinpath(out,"selected.csv"),"w") do io
        for t in 1:q
            row(io,"spectral_l2_aclime",t,6,results[t].stages[6].status===:optimal ? "returned_fixed" : "fixed_failed")
            row(io,"spectral_l2_aclime_tuned",t,selected[t],selected[t]>0 ? "selected" : "no_valid_candidate")
        end
    end
    open(joinpath(out,"auc.csv"),"w") do io
        row(io,"method","target","time_index","observed_area","full_auc","minimum_fpr","maximum_fpr","successful_candidates","expected_candidates")
        for t in 1:q,(target,truth) in ((:actual,actual),(:nominal,nominal)),(method,indices) in
                ((:spectral_l2_aclime,[6]),(:spectral_l2_aclime_tuned,collect(1:6)))
            good=[k for k in indices if results[t].stages[k].status===:optimal && all(isfinite,paths[k][:,:,t])]
            points=[rocpoint(paths[k],truth,t)[1:2] for k in good]
            area,full,lo,hi=observed_auc(points);length(good)==length(indices) || (full=NaN)
            row(io,method,target,t,area,full,lo,hi,length(good),length(indices))
        end
    end
    open(joinpath(out,"complete.toml"),"w") do io
        TOML.print(io,Dict("status"=>failures==0 ? "completed" : "completed_with_failures","candidate_time_failures"=>failures,
            "expected_candidate_time_fits"=>6q,"elapsed_seconds"=>time()-started,"fit_wall_seconds"=>fit_wall_seconds,
            "fit_path_sum_seconds"=>sum(x.path_seconds for x in results),"scoring_seconds"=>time()-score_tick,
            "new_fit_path_sum_seconds"=>sum((x.path_seconds for x in results if !x.restored);init=0.),
            "diagonal_calibration_sum_seconds"=>sum(x.diagonal_seconds for x in results),
            "restored_times"=>count(x->x.restored,results),
            "baseline_fits_repeated"=>0,"point_pilot"=>point_index>0,"source_sha256"=>source_hashes()))
    end
    failures==0
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=option("out","")
    try
        ok=spectral_dataset(option("kind","paper"),parse(Int,option("setting","7")),parse(Int,option("replicate","1")),out;
            point_index=parse(Int,option("point-index","0")),ntasks=parse(Int,option("ntasks","1")),resume_folder=option("resume-folder",""))
        exit(ok ? 0 : 1)
    catch e
        if isdir(out)
            open(joinpath(out,"pipeline-failure.txt"),"w") do io;showerror(io,e,catch_backtrace());end
        end
        rethrow()
    end
end
