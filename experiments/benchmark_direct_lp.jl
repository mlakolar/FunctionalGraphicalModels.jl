# Bounded LP backend comparison on the same observed covariance, not a refit of
# the full simulation campaign. Run under Slurm; retain both solver outcomes.
include("dfgm2020.jl")

function benchmark_lp(out)
    mkpath(out)
    sim = simulate_dfgm2020(MersenneTwister(9301);p=50,n=100,design=:irregular,
        schedule=:T6_9,integration_points=51,evaluation_grid=[.5])
    C = direct_covariance(sim.data,[.5];bandwidths=[.15,.25,.4])
    open(joinpath(out,"backend.csv"),"w") do io
        row(io,"lp_solver","method","repeat","seconds","status","max_constraint_violation","operator_loss","error")
        for lp_solver in (:highs,:clarabel),method in (:aclime,:l2_aclime),rep in 1:2
            tick=time()
            try
                f=fit_direct(C;method,pilot_index=2,bias_constant=1.,K2=.3,
                    noise_variance=.5,lp_solver)
                row(io,lp_solver,method,rep,time()-tick,"ok",
                    maximum(d.max_constraint_violation for d in f.diagnostics),
                    opnorm(f.precisions[:,:,1]-sim.precision[:,:,1]),"")
            catch e
                row(io,lp_solver,method,rep,time()-tick,"failed",NaN,NaN,sprint(showerror,e))
            end
        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run scientific work inside Slurm")
    benchmark_lp(option("out","results/direct-lp"))
end
