@testset "Model-reusing ACLIME paths" begin
    S=[1. .3 0.;.3 1.5 .2;0. .2 .8]
    matrices=cat(S,S+.01I;dims=3)
    ks=[.1,.3,1.]
    for kind in (:linf,:l2),solver in (:highs,:clarabel)
        path=bias_aclime_path(matrices,ks;u=[.04,.05],v=[.01,.02],
            bias=[.02,.03],bias_norm=kind,lp_solver=solver)
        @test all(isempty,path.errors)
        for (index,K2) in enumerate(ks)
            cold=bias_aclime(matrices;u=[.04,.05],v=[.01,.02],bias=[.02,.03],
                bias_norm=kind,K2,lp_solver=solver)
            @test path.pilot≈cold.pilot atol=1e-6
            @test path.budgets[index]≈cold.budgets atol=3e-5 rtol=1e-5
            @test maximum(path.violations[index])<1e-5
            @test all(vec(sum(abs.(path.columns[index]);dims=1)).<=path.budgets[index].+1e-5)
        end
    end
    inverse=bias_aclime_path(reshape(S,3,3,1),ks;u=[0.],v=[0.],bias=[1e-8],bias_norm=:l2)
    @test all(isempty,inverse.errors)
    @test all(W->isapprox(W,inv(S);atol=2e-5),inverse.columns)
    zero=bias_aclime_path(zeros(3,3,2),ks;u=[.1,.2],v=[.1,.2],bias=[.02,.005],bias_norm=:l2)
    @test all(isempty,zero.errors)
    @test all(v->maximum(v)<1e-5,zero.violations)
    @test_throws ArgumentError bias_aclime_path(matrices,Float64[];u=[1.,1.],v=[1.,1.],bias=[.1,.1])
    # Vary both interval directions and retain diagonals across repeated solves.
    O=Matrix{Float64}(I,8,8); O[1,2:end].=.1; O[2:end,1].=.1
    columns=[copy(O),copy(O),copy(O)]
    columns[1][1,2:end].-=.04;columns[2][2:end,1].+=.06
    columns[3][1,2:end].+=.02
    symmetric=FunctionalGraphicalModels._interval_symmetrize_path(columns)
    @test all(isempty,symmetric.errors)
    for (W,Z) in zip(columns,symmetric.solutions)
        cold=interval_symmetrize(W)
        @test issymmetric(Z) && diag(Z)≈diag(W)
        @test opnorm(Z-W,1)≈opnorm(cold-W,1) atol=1e-6
        @test all(min.(W,transpose(W)).-1e-7 .<= Z .<= max.(W,transpose(W)).+1e-7)
        @test opnorm(Z-O,1)<=2opnorm(W-O,1)+1e-6
        @test sum(abs2,Z-O)<=2sum(abs2,W-O)+1e-6
    end
    diagonal=cat(Matrix(Diagonal([1.,2.,3.])),Matrix(Diagonal([1.01,2.01,3.01]));dims=3)
    C=DirectCovariance([.5],[.1,.2],reshape(diagonal,3,3,1,2),
        200,20,fill(200,3,1,2),(;design=:independent))
    for method in (:aclime,:l2_aclime)
        path=fit_direct_path(C,ks;method,pilot_index=2,bias_constant=.5)
        @test all(isempty,path.errors)
        for (index,K2) in enumerate(ks)
            cold=fit_direct(C;method,pilot_index=2,bias_constant=.5,K2)
            @test path.fits[index]!==nothing
            @test path.fits[index].precisions≈cold.precisions atol=2e-4
            @test path.fits[index].settings.K2==K2
            @test maximum(d.max_constraint_violation for d in path.fits[index].diagnostics)<1e-5
        end
    end
    multiple=DirectCovariance([.4,.5],[.1,.2],repeat(C.matrices;outer=(1,1,2,1)),
        200,20,fill(200,3,2,2),(;design=:independent))
    for method in (:aclime,:l2_aclime)
        serial=fit_direct_path(multiple,ks;method,pilot_index=2,bias_constant=.5)
        threaded=fit_direct_path(multiple,ks;method,pilot_index=2,bias_constant=.5,
            ntasks=min(4,Threads.nthreads()))
        @test all(isempty,threaded.errors)
        for index in eachindex(ks)
            @test threaded.fits[index].precisions≈serial.fits[index].precisions atol=1e-5
            @test threaded.fits[index].columns≈serial.fits[index].columns atol=1e-5
            @test length(threaded.fits[index].diagnostics)==2
        end
    end
end
