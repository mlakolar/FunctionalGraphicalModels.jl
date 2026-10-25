@testset "Spectral-calibrated multiscale program P" begin
    F=FunctionalGraphicalModels
    c=F.spectral_calibration(200,3,8;L=1.,kappa_bound=2.,h_max=.25)
    @test c.I==ceil(Int,log2(1600))+1
    @test length(c.bandwidths)<=c.Jstar
    @test c.h_d in c.bandwidths
    @test c.ell≈(3log(3)+log(600)+log(c.Jstar))/200
    @test c.q≈1 ./ (8 .*c.bandwidths)
    nonunit=F.spectral_calibration(200,3,8;L=1.,kappa_bound=2.,h_max=.5,domain_length=2.)
    @test nonunit.q≈2 ./ (8 .*nonunit.bandwidths)
    @test nonunit.h_d≈min(.5,(2sqrt(nonunit.ell)/8)^(1/3))
    @test_throws ArgumentError F.spectral_calibration(200,3,8;L=1.,kappa_bound=2.,h_max=.25,bandwidths=[.1])
    # Constant Gaussian covariance kernel, independent schedules, no bias/error:
    # true columns satisfy every constraint with b=0, at all supplied widths.
    S=[1. .2; .2 1.5];truth=inv(S);A=cat(S,S;dims=3)
    d=diag(S);B=4 ./d;u=[.001,.002];v=[.001,.002];bias=[.0001,.001]
    for j in 1:2,h in 1:2
        @test norm(truth[:,j])<=B[j]
        @test maximum(abs.(A[:,:,h]*truth[:,j]-Float64.(1:2 .== j)))<1e-12
    end
    path=F.spectral_l2_aclime_path(A,[.1,.2];d,u,v,bias,kappa_bound=2.,noise_variance=.5)
    @test all(r->r.status===:optimal,path)
    for r in path,j in 1:2
        beta=r.columns[:,j];t=r.budgets[j];bs=r.bias_vectors[:,:,j]
        @test norm(beta,1)<=t+1e-5
        @test norm(beta)<=B[j]+1e-5
        for h in 1:2
            lambda=r.Kscore.*sqrt.(d.*B[j]).*(u[h]+.5B[j]*v[h])
            @test all(abs.(A[:,:,h]*beta-Float64.(1:2 .== j)-bs[:,h]).<=lambda.+1e-5)
            @test norm(bs[:,h])<=bias[h]*t+1e-5
            @test norm(bs[:,h])<=bias[h]*B[j]+1e-5
        end
        @test r.column_records[j].max_relative_violation<1e-5
    end
    # Variance transport at fixed widths: A,d,bbar multiply by c², B divides
    # by c². sigma²*B and sqrt(d*B), hence the score allowance, stay unchanged.
    scale=3.;scaled=F.spectral_l2_aclime(scale.*A;Kscore=.1,d=scale.*d,u,v,
        bias=scale.*bias,kappa_bound=2.,noise_variance=scale*.5)
    @test scaled.status===:optimal
    @test scaled.columns≈path[1].columns./scale atol=2e-5
    indefinite=cat([-1. .1;.1 1.],[-1. .1;.1 1.];dims=3)
    r=F.spectral_l2_aclime(indefinite;d=[1.,1.],u,v,bias,kappa_bound=2.,Kscore=.1)
    @test r.status===:optimal
    @test r.columns[1,1]<0 # raw indefinite A was not clipped
    infeasible=F.spectral_l2_aclime(zeros(2,2,1);d=[1.,1.],u=[.001],v=[0.],bias=[0.],kappa_bound=2.,Kscore=.1)
    @test infeasible.status===:infeasible_zero_fallback
    @test all(x->x.termination=="INFEASIBLE" && x.dual_status=="INFEASIBILITY_CERTIFICATE",infeasible.column_records)
    @test all(iszero,infeasible.columns)
    analytic=F.spectral_l2_aclime(A;d,u=[2.,2.],v,bias,kappa_bound=2.,Kscore=2.)
    @test all(x->x.status===:optimal_analytic_zero,analytic.column_records)
    # Integration and diagonal calibration; legacy methods before/after new fit.
    cc=F.spectral_calibration(10000,2,100;L=1.,kappa_bound=2.,h_max=.25)
    C=F.DirectCovariance([.5],cc.bandwidths,reshape(repeat(A[:,:,1:1],1,1,length(cc.bandwidths)),2,2,1,length(cc.bandwidths)),
        10000,100,zeros(Int,2,1,length(cc.bandwidths)),(;domain=(0.,1.),design=:independent,boundary=:error,centering=:supplied_mean))
    before=[F.fit_direct(C;method,bias_constant=1.,K2=.1).precisions for method in (:aclime,:l2_aclime,:scaled_lasso)]
    f=F.fit_spectral_direct(C;L=1.,kappa_bound=2.,h_max=.25,Kscore=.1)
    @test f.settings.method===:spectral_l2_aclime
    @test f.diagnostics[1].diagonal==max.(diag(S),10000^-2)
    @test diag(f.precisions[:,:,1])≈f.diagnostics[1].cap.*diag(f.columns[:,:,1])
    @test opnorm(f.precisions[:,:,1],1)<=10000+1e-6
    @test f.diagnostics[1].covariance_psd_displacement==0
    after=[F.fit_direct(C;method,bias_constant=1.,K2=.1).precisions for method in (:aclime,:l2_aclime,:scaled_lasso)]
    @test before==after
    nodes=[FDA.IrregularFunctionalData([[.2,.4],[.3,.6]],[[1.,2.],[.5,1.]];domain=(0.,1.)) for _ in 1:2]
    data=FDA.MultivariateFunctionalData(nodes)
    small=F.spectral_calibration(2,2,2;L=1.,kappa_bound=2.,h_max=.25)
    cached=F.direct_covariance(data,[.5];bandwidths=[.25])
    @test_throws ArgumentError F.spectral_covariance(data,[.5],small;cache=cached)
    built=F.spectral_covariance(data,[.5],small;cache=cached,cache_mean_matches=true)
    @test built.reused_widths==[.25]
    @test !(.25 in built.computed_widths)
    @test built.covariance.matrices[:,:,:,1]==cached.matrices[:,:,:,1]
    wrongdomain=F.spectral_calibration(2,2,2;L=1.,kappa_bound=2.,h_max=.25,domain_length=2.)
    @test_throws ArgumentError F.spectral_covariance(data,[.5],wrongdomain)
end
