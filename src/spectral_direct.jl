"""
    spectral_calibration(N, p, m; L, kappa_bound, h_max, C_K=1, Aconf=1,
                         domain_length=1, grid_mode=:theory, bandwidths=nothing)

Calibration for program (P), not Sun–Zhang scaled lasso. `L` bounds operator
curvature in physical time coordinates; `C_K*L*h^2` is the supplied operator
bias bound. Neither bound is estimated or certified by this function.
`grid_mode=:practical` uses an explicitly supplied reduced grid, augmented by
the same diagonal-calibration width. It is a different empirical protocol.
"""
function spectral_calibration(N::Integer,p::Integer,m::Integer;
        L::Real,kappa_bound::Real,h_max::Real,C_K::Real=1,Aconf::Real=1,
        domain_length::Real=1,grid_mode::Symbol=:theory,bandwidths=nothing)
    N>=2 && p>=1 && m>=2 || throw(ArgumentError("invalid subject/dimension/visit counts"))
    all(x->isfinite(x)&&x>0,(L,C_K,h_max,domain_length,Aconf)) &&
        isfinite(kappa_bound) && kappa_bound>=1 || throw(ArgumentError("invalid spectral calibration"))
    h_max<=domain_length || throw(ArgumentError("h_max exceeds domain length"))
    grid_mode in (:theory,:practical) || throw(ArgumentError("unknown grid mode"))
    grid_mode===:theory && bandwidths!==nothing &&
        throw(ArgumentError("theory grid cannot be replaced by supplied widths"))
    I=ceil(Int,log2(Float64(N)*m))+1; Jstar=I+2
    x=Aconf*log(Float64(p)*N)+log(Jstar); ell=(3log(p)+x)/N
    # In normalized coordinates L_unit=L*domain_length^2 and h_unit=h/T.
    # T*(sqrt(ell)/(m*L_unit))^(1/3) = (T*sqrt(ell)/(m*L))^(1/3).
    h_d=min(Float64(h_max),(domain_length*sqrt(ell)/(m*L))^(1/3))
    hs=grid_mode===:theory ? [h_max*2.0^(-i) for i in 0:I] :
        (bandwidths===nothing ? throw(ArgumentError("practical grid needs widths")) : Float64.(collect(bandwidths)))
    !isempty(hs) && all(h->isfinite(h)&&0<h<=h_max,hs) || throw(ArgumentError("invalid widths"))
    hs=sort!(unique(vcat(Float64.(hs),h_d));rev=true)
    q=domain_length./(m.*hs)
    u=sqrt(ell).*(1 .+q).+ell.*(1 .+q.^2)
    v=sqrt(ell).*q.+ell.*(1 .+q.^2)
    (;N,p,m,L=Float64(L),C_K=Float64(C_K),kappa_bound=Float64(kappa_bound),
      Aconf=Float64(Aconf),domain_length=Float64(domain_length),h_max=Float64(h_max),
      I,Jstar,x,ell,h_d,bandwidths=hs,diagonal_index=findfirst(==(h_d),hs),
      q,u,v,bias=C_K*L.*hs.^2,grid_mode)
end

"""Construct only missing widths of the occupancy covariance and reuse exact
 matching cached widths. The caller must establish observation/seed identity
and explicitly affirm the same supplied mean. Dimensions, schedules, domains,
centering tags and boundary rules are checked here.
"""
function spectral_covariance(data,grid,cal;cache=nothing,cache_mean_matches=false,design=:independent,
        mean=(j,t)->0.,boundary=:error,strata=5)
    hs=cal.bandwidths;g=Float64.(collect(grid));p=cal.p
    FDA.nobs(data)==cal.N && FDA.ncomponents(data)==p || throw(DimensionMismatch("calibration/data differ"))
    dom=FDA.domain(first(data.components))
    dom[2]-dom[1]==cal.domain_length || throw(ArgumentError("calibration domain length differs"))
    if cache!==nothing
        cache_mean_matches || throw(ArgumentError("cached covariance reuse requires affirmation of identical supplied mean"))
        cache.grid==g && cache.n==cal.N && cache.min_visits==cal.m && size(cache.matrices,1)==p ||
            throw(ArgumentError("cached covariance dimensions/grid differ"))
        cache.settings.design==design && cache.settings.boundary==boundary &&
            cache.settings.domain==FDA.domain(first(data.components)) &&
            cache.settings.centering===:supplied_mean || throw(ArgumentError("cache conventions differ"))
    end
    missing=[h for h in hs if cache===nothing || !(h in cache.bandwidths)]
    fresh=isempty(missing) ? nothing : direct_covariance(data,g;bandwidths=missing,design,mean,boundary,strata)
    fresh===nothing || fresh.min_visits==cal.m || throw(ArgumentError("visit resource differs"))
    A=zeros(p,p,length(g),length(hs));counts=zeros(Int,p,length(g),length(hs))
    for (i,h) in enumerate(hs)
        source=cache!==nothing && h in cache.bandwidths ? cache : fresh
        k=findfirst(==(h),source.bandwidths)
        A[:,:,:,i].=source.matrices[:,:,:,k];counts[:,:,i].=source.pair_subjects[:,:,k]
    end
    settings=cache===nothing ? fresh.settings : cache.settings
    C=DirectCovariance(g,copy(hs),A,cal.N,cal.m,counts,settings)
    (;covariance=C,reused_widths=setdiff(hs,missing),computed_widths=missing)
end

"""
    spectral_l2_aclime_path(A, Kscore_values; d, u, v, bias, kappa_bound,
                           noise_variance=1, lp_solver=:highs)

Program (P) on raw, possibly indefinite covariance matrices (p×p×H).
All widths and both bias SOCs occur in each column program. The fixed bounds
are `B=2*kappa_bound ./ d`. No precision pilot or PSD projection is used.
Models are reused only across the explicit score-constant path. Numerical
failures return NaN columns and errors, never a successful zero replacement.
Only an INFEASIBLE status with a dual infeasibility certificate invokes the
separately labeled zero-column fallback. All column statuses are retained.
"""
function spectral_l2_aclime_path(matrices::AbstractArray{<:Real,3},Kscore_values;
        d,u,v,bias,kappa_bound::Real,noise_variance::Real=1,lp_solver::Symbol=:highs)
    p,p2,H=size(matrices);ks=Float64.(collect(Kscore_values))
    p==p2 && p>0 && H>0 && all(isfinite,matrices) &&
        all(issymmetric(matrices[:,:,h]) for h in 1:H) || throw(ArgumentError("invalid symmetric covariances"))
    length(d)==p && all(x->isfinite(x)&&x>0,d) || throw(ArgumentError("invalid calibrated diagonals"))
    length(u)==length(v)==length(bias)==H || throw(DimensionMismatch("one resource per width"))
    all(x->isfinite(x)&&x>=0,u) && all(x->isfinite(x)&&x>=0,v) &&
        all(x->isfinite(x)&&x>=0,bias) || throw(ArgumentError("invalid resources/bias bounds"))
    !isempty(ks) && all(x->isfinite(x)&&x>0,ks) && allunique(ks) || throw(ArgumentError("invalid score path"))
    isfinite(kappa_bound) && kappa_bound>=1 && isfinite(noise_variance) && noise_variance>=0 ||
        throw(ArgumentError("invalid condition/noise bound"))
    A=Float64.(matrices);B=2kappa_bound./Float64.(d)
    all(isfinite,B) || throw(ArgumentError("column bounds overflow"))
    W=[fill(NaN,p,p) for _ in ks];budgets=[fill(NaN,p) for _ in ks]
    records=[NamedTuple[] for _ in ks];bias_vectors=[fill(NaN,p,H,p) for _ in ks]
    for j in 1:p
        base=[sqrt(d[k]*B[j])*(u[h]+noise_variance*B[j]*v[h]) for k in 1:p,h in 1:H]
        all(isfinite,base) || throw(ArgumentError("score allowances overflow"))
        model=nothing;beta=nothing;z=nothing;t=nothing;g=nothing
        upper=nothing;lower=nothing;scale=nothing
        for (index,Kscore) in enumerate(ks)
            started=time(); allowance=Kscore.*base
            # At t=0, beta=b=0 is feasible iff the jth residual is covered.
            # This is an exact global optimum, not a solver-failure fallback.
            if all(allowance[j,h]>=1 for h in 1:H)
                W[index][:,j].=0;budgets[index][j]=0;bias_vectors[index][:,:,j].=0
                push!(records[index],(;column=j,status=:optimal_analytic_zero,
                    termination="ANALYTIC_ZERO",primal_status="FEASIBLE_POINT",dual_status="NOT_RUN",
                    error="",seconds=time()-started,model_setup_seconds=0.,solve_seconds=0.,residual_seconds=0.,budget=0.,l1_violation=0.,ball_violation=0.,
                    score_violation=0.,bias_budget_violation=0.,bias_cap_violation=0.,
                    max_relative_violation=0.,fallback=false))
                continue
            end
            model_setup_seconds=0.;solve_seconds=0.;check_started=0.
            try
                if model===nothing
                    model=_precision_model(;cone=true,lp_solver)
                    JuMP.@variable(model,beta[1:p]);JuMP.@variable(model,z[1:p]>=0)
                    JuMP.@variable(model,0<=t<=sqrt(p)*B[j])
                    JuMP.@variable(model,g[1:p,1:H])
                    JuMP.@constraint(model,[k=1:p],beta[k]<=z[k])
                    JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
                    JuMP.@constraint(model,sum(z)<=t)
                    JuMP.@constraint(model,[B[j];beta] in JuMP.SecondOrderCone())
                    upper=Matrix{Any}(undef,p,H);lower=similar(upper)
                    # Rescale each row using a path-independent constant. The
                    # transformed bias g=b/bbar makes both SOC radii well-scaled.
                    scale=max.(1.,maximum(ks).*base)
                    for h in 1:H
                        JuMP.@constraint(model,[t;g[:,h]] in JuMP.SecondOrderCone())
                        JuMP.@constraint(model,[B[j];g[:,h]] in JuMP.SecondOrderCone())
                        for k in 1:p
                            r=JuMP.@expression(model,(sum(A[k,l,h]*beta[l] for l in 1:p)-bias[h]*g[k,h])/scale[k,h])
                            upper[k,h]=JuMP.@constraint(model,r<=(allowance[k,h]+(k==j))/scale[k,h])
                            lower[k,h]=JuMP.@constraint(model,-r<=(allowance[k,h]-(k==j))/scale[k,h])
                        end
                    end
                    JuMP.@objective(model,Min,t)
                else
                    for h in 1:H,k in 1:p
                        JuMP.set_normalized_rhs(upper[k,h],(allowance[k,h]+(k==j))/scale[k,h])
                        JuMP.set_normalized_rhs(lower[k,h],(allowance[k,h]-(k==j))/scale[k,h])
                    end
                end
                model_setup_seconds=time()-started;solve_started=time()
                JuMP.optimize!(model)
                solve_seconds=time()-solve_started;check_started=time()
                term=JuMP.termination_status(model);primal=string(JuMP.primal_status(model));dual=string(JuMP.dual_status(model))
                if term==JuMP.MOI.INFEASIBLE && JuMP.dual_status(model)==JuMP.MOI.INFEASIBILITY_CERTIFICATE
                    W[index][:,j].=0;budgets[index][j]=0;bias_vectors[index][:,:,j].=0
                    push!(records[index],(;column=j,status=:infeasible_zero_fallback,
                        termination=string(term),primal_status=primal,dual_status=dual,
                        error="solver-certified infeasible column; zero fallback is not a feasible optimum",
                        seconds=time()-started,model_setup_seconds,solve_seconds,residual_seconds=time()-check_started,budget=0.,l1_violation=NaN,ball_violation=NaN,
                        score_violation=NaN,bias_budget_violation=NaN,bias_cap_violation=NaN,
                        max_relative_violation=NaN,fallback=true))
                    continue
                end
                term==JuMP.MOI.OPTIMAL && JuMP.has_values(model) || error("spectral column status $term; primal=$primal; dual=$dual")
                column=JuMP.value.(beta);budget=JuMP.value(t);bs=JuMP.value.(g).*transpose(bias)
                l1=max(0.,norm(column,1)-budget);ball=max(0.,norm(column)-B[j])
                score=0.;bb=0.;bc=0.;relative=max(l1/(1+budget),ball/(1+B[j]))
                for h in 1:H
                    residual=A[:,:,h]*column-Float64.(1:p .== j)-bs[:,h]
                    sv=max.(0.,abs.(residual).-allowance[:,h]);bv=max(0.,norm(bs[:,h])-bias[h]*budget)
                    cv=max(0.,norm(bs[:,h])-bias[h]*B[j])
                    score=max(score,maximum(sv));bb=max(bb,bv);bc=max(bc,cv)
                    relative=max(relative,maximum(sv./(1 .+allowance[:,h])),bv/(1+bias[h]*budget),cv/(1+bias[h]*B[j]))
                end
                all(isfinite,column) && isfinite(relative) && relative<=1e-5 || error("spectral primal residual check failed: $relative")
                W[index][:,j].=column;budgets[index][j]=budget;bias_vectors[index][:,:,j].=bs
                push!(records[index],(;column=j,status=:optimal,termination=string(term),primal_status=primal,dual_status=dual,
                    error="",seconds=time()-started,model_setup_seconds,solve_seconds,residual_seconds=time()-check_started,budget,l1_violation=l1,ball_violation=ball,
                    score_violation=score,bias_budget_violation=bb,bias_cap_violation=bc,
                    max_relative_violation=relative,fallback=false))
            catch exception
                push!(records[index],(;column=j,status=:numerical_failure,termination=model===nothing ? "NOT_BUILT" : string(JuMP.termination_status(model)),
                    primal_status=model===nothing ? "NO_SOLUTION" : string(JuMP.primal_status(model)),
                    dual_status=model===nothing ? "NO_SOLUTION" : string(JuMP.dual_status(model)),
                    error=sprint(showerror,exception),seconds=time()-started,model_setup_seconds,solve_seconds,
                    residual_seconds=check_started==0 ? 0. : time()-check_started,budget=NaN,
                    l1_violation=NaN,ball_violation=NaN,score_violation=NaN,
                    bias_budget_violation=NaN,bias_cap_violation=NaN,max_relative_violation=NaN,fallback=false))
            end
        end
        model===nothing || JuMP.empty!(model)
    end
    [(;Kscore=ks[i],columns=W[i],budgets=budgets[i],bias_vectors=bias_vectors[i],
        column_records=records[i],B,diagonal=Float64.(d),
        status=any(x->x.status===:numerical_failure,records[i]) ? :numerical_failure :
            any(x->x.fallback,records[i]) ? :infeasible_zero_fallback : :optimal) for i in eachindex(ks)]
end

spectral_l2_aclime(A;Kscore::Real=2,kwargs...)=only(spectral_l2_aclime_path(A,[Kscore];kwargs...))

"""Fit program (P) using a `DirectCovariance` on the prescribed grid.
Failure diagnostics and original indefinite covariances are retained. Numerical
failures make the corresponding precision slice NaN. Infeasible-column zero
fallback estimates are explicitly labeled and are not optimal candidates.
"""
function fit_spectral_direct(C::DirectCovariance;L::Real,kappa_bound::Real,h_max::Real,
        C_K::Real=1,Aconf::Real=1,Kscore::Real=2,noise_variance::Real=1,
        grid_mode::Symbol=:theory,bandwidths=nothing,threshold::Real=1e-6,
        lp_solver::Symbol=:highs)
    isfinite(threshold) && threshold>=0 || throw(ArgumentError("invalid adjacency threshold"))
    p,_,q,H=size(C.matrices);domain=get(C.settings,:domain,(0.,1.))
    cal=spectral_calibration(C.n,p,C.min_visits;L,kappa_bound,h_max,C_K,Aconf,
        domain_length=domain[2]-domain[1],grid_mode,bandwidths)
    indices=[findfirst(==(h),C.bandwidths) for h in cal.bandwidths]
    all(i->i!==nothing,indices) || throw(ArgumentError("covariance cache lacks prescribed spectral widths"))
    indices=Int.(indices);P=fill(NaN,p,p,q);W=similar(P);diagnostics=NamedTuple[]
    for t in 1:q
        d=max.(diag(C.matrices[:,:,t,indices[cal.diagonal_index]]),C.n^-2)
        started=time();r=spectral_l2_aclime(C.matrices[:,:,t,indices];Kscore,d,
            u=cal.u,v=cal.v,bias=cal.bias,kappa_bound,noise_variance,lp_solver)
        fitting_seconds=time()-started;W[:,:,t].=r.columns;symmetry_seconds=0.;cap_seconds=0.;cap=NaN
        status=r.status;error=""
        if all(isfinite,r.columns)
            try
                started=time();S=interval_symmetrize(r.columns;lp_solver);symmetry_seconds=time()-started
                started=time();size1=opnorm(S,1);cap=size1==0 ? 1. : min(1.,C.n/size1)
                P[:,:,t].=cap.*S;cap_seconds=time()-started
            catch exception
                status=:numerical_failure;error="symmetrization: "*sprint(showerror,exception)
            end
        end
        push!(diagnostics,(;status,error,column_records=r.column_records,diagonal=d,B=r.B,
            fitting_seconds,symmetry_seconds,cap_seconds,cap,covariance_psd_displacement=0.,
            max_constraint_violation=maximum(x.max_relative_violation for x in r.column_records)))
    end
    adjacency=abs.(P).>threshold
    for t in 1:q,j in 1:p;adjacency[j,j,t]=false;end
    positive=BitVector([all(isfinite,P[:,:,t]) && isposdef(Symmetric(P[:,:,t])) for t in 1:q])
    settings=(;method=:spectral_l2_aclime,calibration=cal,Kscore,noise_variance,
        threshold,symmetrization=:interval_l1,lp_solver,socp_solver=:clarabel,
        socp_tolerance=1e-7,psd_repair=false,precision_pilot=false)
    DirectPrecisionFit(copy(C.grid),P,W,adjacency,positive,diagnostics,C,settings)
end
