include("dfgm2020.jl")

function benchmark_common(out)
    mkpath(out)
    open(joinpath(out,"common.csv"),"w") do io
        row(io,"mode","p","n","schedule","seconds","minimum_eigenvalue","status","error")
        for (mode,p,n,schedule,design) in ((:pilot,50,100,:T25,:common),
                (:paper,10,20,:T10,:common),(:paper,10,20,:T6_9,:irregular))
            sim = simulate_dfgm2020(MersenneTwister(811);p,n,design,schedule,integration_points=51)
            tick=time()
            try
                S,details = reconstructed_covariances(sim.data,sim.grid;mode,work_points=31)
                row(io,mode,p,n,schedule,time()-tick,minimum(eigmin(Symmetric(S[:,:,t])) for t in axes(S,3)),"ok","")
                serialize(joinpath(out,"selection-$mode-$design.bin"),details)
            catch e
                row(io,mode,p,n,schedule,time()-tick,NaN,"failed",sprint(showerror,e))
            end
        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run scientific work inside Slurm")
    benchmark_common(option("out","results/common-benchmark"))
end
