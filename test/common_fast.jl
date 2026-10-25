include(joinpath(@__DIR__,"..","experiments","dfgm_common_fast.jl"))
@testset "Common-grid sufficient statistics preserve FDA mean CV and covariance" begin
    rng=MersenneTwister(82104)
    for T in (6,10)
        times=sort!([(b+rand(rng))/T for b in 0:T-1])
        values=[randn(rng).*sinpi.(times)+.3randn(rng,T) for _ in 1:8]
        fd=FDA.IrregularFunctionalData([copy(times) for _ in 1:8],values;domain=(0.,1.))
        grid=sort!(unique(vcat(collect(range(0,1;length=15)),times)))
        selector=FDA.CurveCV([.1,.2,.3,.4];folds=8)
        oldmean=FDA.fit(FDA.SmoothedMean,fd;grid,bandwidth=selector,
            kernel=KernelSmoothing.GaussianKernel(),weighting=:subj)
        fastmean=_dfgm_common_fast_mean(fd,grid;bandwidth=selector)
        @test oldmean.values==fastmean.values
        @test oldmean.selection.errors≈fastmean.selection.errors atol=1e-10 rtol=0
        @test oldmean.bandwidth==fastmean.bandwidth
        @test fastmean.folds==8 && fastmean.weighting==:subj
        oldcov=_dfgm_common_covariance(fd,oldmean,grid;bandwidth=selector)
        fastcov=_dfgm_common_fast_covariance(fd,fastmean,grid;bandwidth=selector)
        @test oldcov.surface≈fastcov.surface atol=1e-10 rtol=0
        @test oldcov.noise_variance≈fastcov.noise_variance atol=1e-10 rtol=0
        @test oldcov.selection.errors≈fastcov.selection.errors atol=1e-10 rtol=0
        @test oldcov.bandwidth==fastcov.bandwidth
    end
end
