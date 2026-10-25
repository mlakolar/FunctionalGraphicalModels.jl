# Same controlled datasets as direct_stress.jl; compare the iid default with
# observable held-out tuning. Preserve all default and failed outcomes.
include("direct_stress.jl")

function scaled_calibration(index,out;replicates=3)
    replicates>=1 || error("replicates must be positive")
    configurations=[(;degree,m,curvature,noise_variance) for degree in (2,10)
        for m in (8,30) for curvature in (0.,4.) for noise_variance in (.5,2.)]
    config=configurations[index]
    mkpath(out)
    isfile(joinpath(out,"fits.csv")) && error("existing calibration results; choose another attempt root")
    candidates=(:quantile,.2,.35,.5,.75,:visit_adjusted)
    errors=0
    open(joinpath(out,"settings.toml"),"w") do io
        TOML.print(io,Dict("index"=>index,"n"=>200,"p"=>50,"replicates"=>replicates,
            "degree"=>config.degree,"m"=>config.m,"curvature"=>config.curvature,
            "noise_variance"=>config.noise_variance,
            "candidates"=>collect(string.(candidates)),"selection"=>"heldout covariance precision residual",
            "visit_adjustment"=>"sqrt(2*log(p)/ntrain)*(1+1/(m*0.2)); empirical heuristic, not a theorem"))
    end
    open(joinpath(out,"fits.csv"),"w") do io
        row(io,"replicate","candidate","lambda0","status","seconds","validation_loss","selected",
            "operator","normalized_squared_frobenius","psd_displacement","cap","error")
        for rep in 1:replicates
            sim=stress_data(MersenneTwister(2026102500+10000index+rep);config...)
            Ct=direct_covariance(sim.data[1:160],[.5];bandwidths=[.1,.2,.3])
            Cv=direct_covariance(sim.data[161:200],[.5];bandwidths=[.2])
            fits=[]
            for (candidate,value) in enumerate(candidates)
                lambda0=value===:visit_adjusted ? sqrt(2log(50)/160)*(1+1/(config.m*.2)) : value
                tick=time()
                try
                    fit=fit_direct(Ct;method=:scaled_lasso,pilot_index=2,bias_constant=1.,lambda0,
                        noise_variance=config.noise_variance)
                    loss=norm(Cv.matrices[:,:,1,1]*fit.precisions[:,:,1]-I)
                    D=fit.precisions[:,:,1]-sim.precision
                    push!(fits,(;candidate,lambda0,seconds=time()-tick,loss,operator=opnorm(D),
                        frobenius=sum(abs2,D)/50,displacement=fit.diagnostics[1].covariance_psd_displacement,
                        cap=fit.diagnostics[1].cap))
                catch e
                    errors+=1
                    row(io,rep,candidate,lambda0,"failed",time()-tick,NaN,false,NaN,NaN,NaN,NaN,sprint(showerror,e))
                end
            end
            isempty(fits) && continue
            best=argmin([f.loss for f in fits])
            for (k,f) in enumerate(fits)
                row(io,rep,f.candidate,f.lambda0,"ok",f.seconds,f.loss,k==best,f.operator,
                    f.frobenius,f.displacement,f.cap,"")
            end
        end
    end
    open(joinpath(out,"complete.toml"),"w") do io
        TOML.print(io,Dict("status"=>errors==0 ? "completed" : "completed_with_failures",
            "expected_fits"=>6replicates,"fit_failures"=>errors))
    end
    println("scaled calibration index=$index failures=$errors")
    errors==0
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run scientific work inside Slurm")
    exit(scaled_calibration(parse(Int,option("setting","1")),option("out","results/scaled-calibration");
        replicates=parse(Int,option("replicates","3"))) ? 0 : 1)
end
