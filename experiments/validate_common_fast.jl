# Correctness and timing gate for exact common-grid product geometry.
# Full-size inputs are regenerated with the same seed; raw cached covariance
# targets and the original source remain untouched. Run only under Slurm.
include("dfgm2020.jl")
using Test

function validate_common_fast(root,out)
    mkpath(out)
    open(joinpath(out,"mean-benchmark.csv"),"w") do meanio
    row(meanio,"scope","visits","subjects","old_seconds","fast_seconds","max_mean_difference","CV_difference","same_bandwidth")
    open(joinpath(out,"benchmark.csv"),"w") do io
        row(io,"scope","visits","subjects","old_seconds","fast_seconds",
            "max_surface_difference","noise_difference","CV_difference","same_bandwidth")
        @testset "Exact common product geometry versus retained reference" begin
            rng=MersenneTwister(82104)
            for T in (6,10,25,50)
                times=sort!([(b+rand(rng))/T for b in 0:T-1])
                values=[randn(rng).*sinpi.(times)+.3randn(rng,T) for _ in 1:8]
                fd=FDA.IrregularFunctionalData([copy(times) for _ in 1:8],values;domain=(0.,1.))
                grid=sort!(unique(vcat(collect(range(0,1;length=15)),times)))
                selector=FDA.CurveCV([.1,.2,.3,.4];folds=8)
                tick=time();meanfit=FDA.fit(FDA.SmoothedMean,fd;grid,bandwidth=selector,
                    kernel=KernelSmoothing.GaussianKernel(),weighting=:subj);oldmeansec=time()-tick
                tick=time();fastmean=_dfgm_common_fast_mean(fd,grid;bandwidth=selector);fastmeansec=time()-tick
                md=maximum(abs.(meanfit.values-fastmean.values));mc=maximum(abs.(meanfit.selection.errors-fastmean.selection.errors))
                @test md<1e-10 && mc<1e-10
                @test meanfit.bandwidth==fastmean.bandwidth
                row(meanio,"small",T,8,oldmeansec,fastmeansec,md,mc,meanfit.bandwidth==fastmean.bandwidth)
                oldcache=Dict();fastcache=Dict()
                for bandwidth in (.2,FDA.CurveCV([.1,.2,.3,.4];folds=8))
                    tick=time();old=_dfgm_common_covariance(fd,meanfit,grid;bandwidth,cache=oldcache)
                    oldsec=time()-tick
                    tick=time();fast=_dfgm_common_fast_covariance(fd,meanfit,grid;bandwidth,cache=fastcache)
                    fastsec=time()-tick
                    delta=maximum(abs.(old.surface-fast.surface));noise=abs(old.noise_variance-fast.noise_variance)
                    cv=bandwidth isa FDA.CurveCV ? maximum(abs.(old.selection.errors-fast.selection.errors)) : 0.
                    @test delta<1e-9
                    @test noise<1e-9
                    @test cv<1e-9
                    @test old.bandwidth==fast.bandwidth
                    row(io,"small",T,8,oldsec,fastsec,delta,noise,cv,old.bandwidth==fast.bandwidth)
                end
            end
        end
        # Real common p100/T50 covariance from the superseded cold canary.
        seed=2026102500+10000*6+1
        sim=simulate_dfgm2020(MersenneTwister(seed);n=200,p=100,design=:common,
            schedule=:T50,integration_points=101)
        saved=deserialize(joinpath(root,"outputs","production-canary","setting-6-rep-1","covariances.bin"))
        @test sim.grid==saved.direct.grid
        fd=FDA.component(sim.data,1)
        workgrid=sort!(unique(vcat(collect(range(0,1;length=31)),sim.grid)))
        selector=FDA.CurveCV([.1,.2,.3,.4];folds=200)
        tick=time();meanfit=FDA.fit(FDA.SmoothedMean,fd;grid=workgrid,bandwidth=selector,
            kernel=KernelSmoothing.GaussianKernel(),weighting=:subj);oldmeansec=time()-tick
        tick=time();fastmean=_dfgm_common_fast_mean(fd,workgrid;bandwidth=selector);fastmeansec=time()-tick
        md=maximum(abs.(meanfit.values-fastmean.values));mc=maximum(abs.(meanfit.selection.errors-fastmean.selection.errors))
        @test md<1e-10 && mc<1e-10
        @test meanfit.bandwidth==fastmean.bandwidth
        row(meanio,"full_size_component",50,200,oldmeansec,fastmeansec,md,mc,meanfit.bandwidth==fastmean.bandwidth)
        @testset "Full-size component and AIC scores" begin
            tick=time();old=_dfgm_common_covariance(fd,meanfit,workgrid;bandwidth=selector);oldsec=time()-tick
            tick=time();fast=_dfgm_common_fast_covariance(fd,meanfit,workgrid;bandwidth=selector);fastsec=time()-tick
            delta=maximum(abs.(old.surface-fast.surface));noise=abs(old.noise_variance-fast.noise_variance)
            cv=maximum(abs.(old.selection.errors-fast.selection.errors))
            @test delta<1e-9 && noise<1e-9 && cv<1e-9
            @test old.bandwidth==fast.bandwidth
            oldfp=FDA.fit(FDA.FPCA,fd;method=:pace,grid=workgrid,ncomponents=FDA.AIC(),maxK=5,mean=meanfit,covariance=old)
            fastfp=FDA.fit(FDA.FPCA,fd;method=:pace,grid=workgrid,ncomponents=FDA.AIC(),maxK=5,mean=meanfit,covariance=fast)
            @test length(oldfp.eigenvalues)==length(fastfp.eigenvalues)
            @test oldfp.scores*transpose(FDA.basismatrix(oldfp.basis,sim.grid))≈
                fastfp.scores*transpose(FDA.basismatrix(fastfp.basis,sim.grid)) atol=1e-8 rtol=1e-8
            row(io,"full_size_component",50,200,oldsec,fastsec,delta,noise,cv,old.bandwidth==fast.bandwidth)
        end
        @testset "Full p100 reconstructed covariance and all ranks" begin
            tick=time();S,selections=reconstructed_covariances(sim.data,sim.grid;
                mode=:paper,work_points=31,ntasks=4,common_fast=true);seconds=time()-tick
            oldselections=deserialize(joinpath(root,"outputs","production-canary","setting-6-rep-1","fpca_selections.bin"))
            delta=maximum(abs.(S-saved.reconstructed))
            @test delta<1e-8
            @test getproperty.(selections,:rank)==getproperty.(oldselections,:rank)
            row(io,"full_p100_reconstruction",50,200,NaN,seconds,delta,NaN,NaN,true)
        end
    end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run under Slurm")
    validate_common_fast(option("root",""),option("out","results/common-fast-check"))
end
