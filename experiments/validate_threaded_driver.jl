# Check deterministic independent-component/time parallelism, then exercise
# the full cached driver. All work belongs to one Slurm allocation.
include("dfgm2020.jl")
using Test

function validate_threaded(root,out)
    mkpath(out)
    @test Threads.nthreads()>=4
    @testset "Threaded covariance reconstruction equals serial" begin
        for (design,schedule) in ((:common,:T10),(:irregular,:T6_9))
            sim=simulate_dfgm2020(MersenneTwister(7206);n=20,p=10,design,schedule,integration_points=31)
            serial,ss=reconstructed_covariances(sim.data,sim.grid;mode=:paper,work_points=15,ntasks=1)
            parallel,ps=reconstructed_covariances(sim.data,sim.grid;mode=:paper,work_points=15,ntasks=4)
            @test parallel≈serial atol=1e-9 rtol=1e-9
            @test getproperty.(ps,:rank)==getproperty.(ss,:rank)
        end
    end
    covariance=deserialize(joinpath(root,"outputs","production-canary","setting-6-rep-1","covariances.bin"))
    C=covariance.direct
    center=argmin(abs.(C.grid.-.5)); selected=sort!(unique([center,min(center+1,length(C.grid))]))
    localC=DirectCovariance(C.grid[selected],C.bandwidths,C.matrices[:,:,selected,:],
        C.n,C.min_visits,C.pair_subjects[:,selected,:],C.settings)
    open(joinpath(out,"threaded-benchmark.csv"),"w") do io
        row(io,"p","method","time_points","candidates","serial_seconds","threaded_seconds","max_precision_difference","max_violation")
        @testset "Independent time-point solves on full-size input" begin
            for method in (:aclime,:l2_aclime)
                tick=time();serial=fit_direct_path(localC,[.1,1.];method,pilot_index=2,bias_constant=1.,noise_variance=.5)
                serial_seconds=time()-tick
                tick=time();parallel=fit_direct_path(localC,[.1,1.];method,pilot_index=2,bias_constant=1.,noise_variance=.5,ntasks=4)
                threaded_seconds=time()-tick
                @test all(isempty,serial.errors) && all(isempty,parallel.errors)
                delta=maximum(maximum(abs.(serial.fits[k].precisions-parallel.fits[k].precisions)) for k in 1:2)
                violation=maximum(d.max_constraint_violation for f in parallel.fits for d in f.diagnostics)
                @test delta<1e-4
                row(io,size(C.matrices,1),method,length(selected),2,serial_seconds,threaded_seconds,delta,violation)
            end
        end
    end
    @testset "Four-thread cached paper-mode pipeline" begin
        for setting_index in (1,7)
            folder=joinpath(out,"setting-$setting_index-rep-1")
            success=worker(setting_index,1,folder;mode=:paper,n_override=30,p_override=10,
                integration_points=31,nlambda=3,work_points=15,bias_constants=[.01,.1,1.,10.,20.],
                path_reuse=true,ntasks=4,
                covariance_cache=joinpath(root,"outputs","bias-menu-canary","setting-$setting_index-rep-1"))
            receipt=TOML.parsefile(joinpath(folder,"complete.toml"))
            @test success == (receipt["fit_failures"]==0)
            @test receipt["status"] == (success ? "completed" : "completed_with_failures")
            # A deterministic non-OPTIMAL column may remain excluded. Validate
            # failure accounting rather than treating it as a usable estimate.
            @test receipt["fit_failures"] <= 1
        end
    end
    open(joinpath(out,"tasks.tsv"),"w") do io
        for s in (1,7);println(io,"$s\t1\tpaper\t30\t10\t31\t3\t15");end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run under Slurm")
    validate_threaded(option("root",""),option("out","results/threaded-driver"))
end
