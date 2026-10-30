@testset "Native multiscale L2 workspace sequences" begin
    F=FunctionalGraphicalModels
    S=[1. .2 0.;.2 1.5 .1;0. .1 .8]
    ks=[.3,.1,.6] # arbitrary caller order must be preserved
    first_problem=(matrices=cat(S,S+.01I;dims=3),u=[.02,.04],v=[.01,.02],bias=[.01,.03],pilot_index=2)
    second_problem=(matrices=cat(S+.05I,S-.01I;dims=3),u=[.03,.02],v=[.02,.01],bias=[.02,.01],pilot_index=1)
    # Same shape but different sparsity of the covariance values.
    third_problem=(matrices=cat(Matrix(Diagonal(diag(S))),S;dims=3),u=[.02,.04],v=[.01,.02],bias=[.01,.03],pilot_index=2)
    problems=[first_problem,second_problem,third_problem]
    path=F.multiscale_l2_aclime_path(problems,ks;noise_variance=.5)
    independent=F.multiscale_l2_aclime_path(problems,ks;noise_variance=.5,reuse=false)
    @test path.diagnostics.workspace_builds==3
    @test path.diagnostics.matrix_updates==6
    @test path.diagnostics.rhs_updates==18
    @test !path.diagnostics.iterate_warm_start
    @test independent.diagnostics.workspace_builds==27
    for (a,spec) in enumerate(problems),(b,penalty) in enumerate(ks)
        r=path.results[a][b]
        cold=F.bias_aclime(spec.matrices;u=spec.u,v=spec.v,bias=spec.bias,
            pilot_index=spec.pilot_index,K2=penalty,noise_variance=.5,bias_norm=:l2)
        @test r.status===:optimal
        @test r.penalty==penalty
        @test r.pilot≈cold.pilot atol=1e-6
        @test r.budgets≈cold.budgets atol=5e-5 rtol=1e-5
        @test r.budgets≈independent.results[a][b].budgets atol=5e-5 rtol=1e-5
        @test all(x->x.max_relative_violation<=1e-5,r.column_records)
        @test issymmetric(F.interval_symmetrize(r.columns))
    end
    # Changing H requires a fresh structure; going back must also be safe.
    single=(matrices=first_problem.matrices[:,:,1:1],u=[.02],v=[.01],bias=[.01])
    shapes=F.multiscale_l2_aclime_path([first_problem,single,second_problem],ks)
    @test shapes.diagnostics.workspace_builds==9
    @test all(r->r.status===:optimal,Iterators.flatten(shapes.results))
    spectral=[(;matrices=s.matrices,u=s.u,v=s.v,bias=s.bias,d=diag(s.matrices[:,:,1]),kappa_bound=2.) for s in problems]
    path=F.multiscale_l2_aclime_path(spectral,ks;calibration=:spectral,noise_variance=.5)
    @test path.diagnostics.workspace_builds>=3 # failed columns rebuild before continuing
    for (a,spec) in enumerate(spectral),(b,penalty) in enumerate(ks)
        cold=F.spectral_l2_aclime(spec.matrices;d=spec.d,u=spec.u,v=spec.v,bias=spec.bias,
            kappa_bound=2.,noise_variance=.5,Kscore=penalty)
        r=path.results[a][b]
        @test (r.status===:optimal)==(cold.status===:optimal)
        for j in 1:3
            if cold.column_records[j].status===:optimal
                @test r.budgets[j]≈cold.budgets[j] atol=5e-5 rtol=1e-5
            else
                @test cold.column_records[j].fallback
                @test all(isnan,r.columns[:,j])
                @test r.column_records[j].termination=="PRIMAL_INFEASIBLE"
            end
        end
    end
    # Analytic zero is a feasible optimum. Infeasibility is NOT a zero estimate.
    impossible=(matrices=zeros(2,2,1),d=ones(2),u=[0.],v=[0.],bias=[0.],kappa_bound=2.)
    recovered=(;impossible...,u=[2.])
    failures=F.multiscale_l2_aclime_path([impossible,recovered],[.5,1.];calibration=:spectral)
    @test all(r->r.status===:failure,failures.results[1])
    @test all(r->all(isnan,r.columns),failures.results[1])
    @test all(r->r.status===:optimal,failures.results[2])
    @test all(r->maximum(abs,r.columns)<1e-5,failures.results[2])
    indefinite=(;first_problem...,matrices=cat([-1. .1;.1 1.],[-1. .1;.1 1.];dims=3))
    raw=F.multiscale_l2_aclime_path([indefinite],[.1])
    @test raw.results[1][1].status===:optimal
    @test_throws ArgumentError F.multiscale_l2_aclime_path(problems,[0.])
    @test_throws ArgumentError F.multiscale_l2_aclime_path(problems,[.1,.1])
    @test_throws ArgumentError F.multiscale_l2_aclime_path([(;first_problem...,bias=[0.,.1])],ks)
end

@testset "Explicit multiscale bandwidth and penalty grid" begin
    F=FunctionalGraphicalModels
    S=[1. .1;.1 1.2];hs=[.1,.2,.3]
    C=F.DirectCovariance([.5],hs,reshape(cat(S,S+.01I,S+.02I;dims=3),2,2,1,3),
        200,30,fill(200,2,1,3),(;domain=(0.,1.),design=:independent))
    sets=[[1,2],[2,3]];penalties=[.1,.3];constants=[.2,.5]
    grid=F.fit_multiscale_l2_path(C,penalties;bandwidth_sets=sets,bias_constants=constants)
    @test size(grid.fits)==(4,2)
    @test all(isempty,grid.errors)
    @test grid.solver_diagnostics[1].workspace_builds==2
    for a in 1:4,b in 1:2
        cfg=grid.configurations[a];indices=cfg.bandwidth_indices
        sub=F.DirectCovariance(C.grid,C.bandwidths[indices],C.matrices[:,:,:,indices],C.n,C.min_visits,C.pair_subjects[:,:,indices],C.settings)
        cold=F.fit_direct(sub;method=:l2_aclime,pilot_index=1,bias_constant=cfg.bias_constant,K2=penalties[b])
        @test grid.fits[a,b].precisions≈cold.precisions atol=1e-4 rtol=1e-4
        @test grid.fits[a,b].settings.bandwidth_indices==indices
    end
    @test_throws ArgumentError F.fit_multiscale_l2_path(C,penalties;bandwidth_sets=[[1]],pilot_indices=[2],bias_constants=constants)
end

@testset "Analytic zero and lazy covariance updates" begin
    F=FunctionalGraphicalModels
    # With theta=2*kappa/d=4 and u=1, zero is feasible exactly at penalty .5.
    a=(matrices=reshape([1. .1;.1 1.],2,2,1),d=ones(2),u=[1.],v=[0.],bias=[.01],kappa_bound=2.)
    b=(;a...,matrices=reshape([2. .3;.3 1.5],2,2,1))
    ks=[.5,.1,1.,.2] # skip before, between and after actual solves
    path=F.multiscale_l2_aclime_path([a,b],ks;calibration=:spectral)
    cold=F.multiscale_l2_aclime_path([a,b],ks;calibration=:spectral,reuse=false)
    @test path.diagnostics.analytic_zeros==8
    @test path.diagnostics.workspace_builds==2
    @test path.diagnostics.matrix_updates==2
    for i in 1:2,k in eachindex(ks)
        r=path.results[i][k]
        @test r.status===:optimal
        @test r.budgets≈cold.results[i][k].budgets atol=1e-5
        if k in (1,3)
            @test iszero(r.columns) && iszero(r.budgets)
            @test all(x->x.termination=="ANALYTIC_ZERO" && x.iterations==0,r.column_records)
        else
            @test all(x->x.termination=="SOLVED",r.column_records)
            @test all(>(0.),r.budgets)
        end
    end
    # Strict certificate: a sub-boundary allowance is solved, never rounded up.
    below=F.multiscale_l2_aclime_path([a],[prevfloat(.5)];calibration=:spectral)
    @test below.diagnostics.analytic_zeros==0
    # The same certificate applies with the data-based pilot calibration.
    pilot=(matrices=ones(1,1,1),u=[.1],v=[0.],bias=[.01])
    p=F.multiscale_l2_aclime_path([pilot],[100.,.1];lp_solver=:clarabel)
    @test p.results[1][1].column_records[1].termination=="ANALYTIC_ZERO"
    @test iszero(p.results[1][1].columns)
    @test p.results[1][2].column_records[1].termination=="SOLVED"
end

@testset "Unscaled primal guard rejects incompatible candidates" begin
    F=FunctionalGraphicalModels
    A=ones(1,1,1);base=fill(.01,1,1);bias=[.01]
    beta=[1.];budget=1.;g=ones(1,1)
    valid=F._l2_primal_residual(A,base,bias,1,1.,nothing,beta,budget,g)
    @test valid.relative<=1e-5
    # Independently violate auxiliary cone, score and column constraints.
    cone=F._l2_primal_residual(A,base,bias,1,1.,nothing,beta,budget,fill(1.0001,1,1))
    @test cone.relative>1e-5
    score=F._l2_primal_residual(A,zeros(1,1),[1.],1,1.,nothing,beta,budget,g)
    @test score.relative>1e-5
    column=F._l2_primal_residual(A,base,bias,1,1.,nothing,[2.],budget,g)
    @test column.relative>1e-5
end

@testset "Saved covariance numerical regression with strict retries" begin
    F=FunctionalGraphicalModels
    bytes=read(joinpath(@__DIR__,"fixtures","multiscale-id1-full8.f64"))
    @test length(bytes)==160000
    A=reshape(reinterpret(Float64,bytes),50,50,8)
    widths=[.025,.05,.075,.1,.15,.2,.3,.4]
    ell=(3log(50)+log(8))/200;q=1 ./ (8 .*widths)
    u=sqrt(ell).*(1 .+q).+ell.*(1 .+q.^2)
    v=sqrt(ell).*q.+ell.*(1 .+q.^2)
    problem=(;matrices=A,u,v,bias=0.1 .* widths.^2,pilot_index=5)
    native=F.multiscale_l2_aclime_path([problem],[.05];noise_variance=.5,
        diagonal_floor=200^-2,lp_solver=:clarabel)
    r=native.results[1][1]
    @test r.status===:optimal
    @test all(x->x.max_relative_violation<=1e-5,r.column_records)
    @test native.diagnostics.guard_retries==sum(x.retry_count for x in r.column_records)
    for record in r.column_records
        attempts=record.attempt_records
        @test length(attempts)==record.retry_count+1
        @test first(attempts).tolerance==1e-7
        @test last(attempts).termination=="SOLVED"
        @test last(attempts).max_relative_violation<=1e-5
        if record.retry_count>0
            @test first(attempts).termination=="SOLVED"
            @test first(attempts).max_relative_violation>1e-5
            @test all(i->attempts[i].tolerance<attempts[i-1].tolerance,2:length(attempts))
        end
    end
    cold=F.bias_aclime(A;u,v,bias=problem.bias,pilot_index=5,K2=.05,
        noise_variance=.5,bias_norm=:l2,diagonal_floor=200^-2,lp_solver=:clarabel)
    @test maximum(abs.(r.budgets-cold.budgets)./(1 .+abs.(cold.budgets)))<=1e-4
end
