function _cold_aclime_column(A,j,base,bias,K2,bias_norm,lp_solver)
    p,_,H=size(A)
    model=_precision_model(;cone=bias_norm===:l2,lp_solver)
    JuMP.@variable(model,beta[1:p])
    JuMP.@variable(model,z[1:p]>=0)
    JuMP.@expression(model,t,sum(z))
    JuMP.@constraint(model,[k=1:p],beta[k]<=z[k])
    JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
    b=bias_norm===:l2 ? JuMP.@variable(model,[1:p,1:H]) : nothing
    for h in 1:H
        if bias_norm===:l2
            JuMP.@constraint(model,[bias[h]*t;b[:,h]] in JuMP.SecondOrderCone())
        end
        for k in 1:p
            r=sum(A[k,l,h]*beta[l] for l in 1:p)-(k==j)
            if bias_norm===:l2
                JuMP.@constraint(model,r-b[k,h]<=K2*base[k,h])
                JuMP.@constraint(model,-r+b[k,h]<=K2*base[k,h])
            else
                JuMP.@constraint(model,r<=K2*base[k,h]+bias[h]*t)
                JuMP.@constraint(model,-r<=K2*base[k,h]+bias[h]*t)
            end
        end
    end
    JuMP.@objective(model,Min,t)
    _solve_checked!(model)
    column=JuMP.value.(beta);budget=JuMP.value(t)
    violation=max(0,norm(column,1)-budget)
    for h in 1:H
        r=A[:,:,h]*column-Float64.(1:p .== j)
        if bias_norm===:l2
            bh=JuMP.value.(b[:,h])
            violation=max(violation,norm(bh)-bias[h]*budget,maximum(abs.(r-bh)-K2.*base[:,h]))
        else
            violation=max(violation,maximum(abs.(r)-K2.*base[:,h].-bias[h]*budget))
        end
    end
    JuMP.empty!(model)
    (;column,budget,violation)
end

"""
    bias_aclime_path(matrices, K2_values; u, v, bias, kwargs...)

Solve the same programs as `bias_aclime` along a finite K2 path. The pilot
is independent of K2 and is solved once. Each second-stage column model is
built once, then only its normalized right-hand sides change. Failed path
members have an error string and must not be used as estimates. No failed
member is replaced by a reference graph. Solver tolerances are unchanged.
Nonunique optimal solutions can differ from independent cold solves.
"""
function bias_aclime_path(matrices::AbstractArray{<:Real,3}, K2_values;
        u,v,bias,pilot_index::Integer=1,K1::Real=2,noise_variance::Real=1,
        bias_norm::Symbol=:linf,diagonal_floor::Real=1e-8,lp_solver::Symbol=:highs)
    ks=Float64.(collect(K2_values)); p,p2,H=size(matrices)
    p==p2 && p>0 && H>0 && all(isfinite,matrices) || throw(ArgumentError("invalid covariance matrices"))
    all(issymmetric(matrices[:,:,h]) for h in 1:H) || throw(ArgumentError("matrices must be symmetric"))
    1<=pilot_index<=H || throw(ArgumentError("invalid pilot_index"))
    length(u)==length(v)==length(bias)==H || throw(DimensionMismatch("one resource/bias per bandwidth"))
    !isempty(ks) && all(x->isfinite(x)&&x>0,ks) || throw(ArgumentError("invalid K2 path"))
    all(x->isfinite(x)&&x>=0,u) && all(x->isfinite(x)&&x>=0,v) &&
        all(x->isfinite(x)&&x>0,bias) || throw(ArgumentError("invalid resources or bias envelopes"))
    all(x->isfinite(x)&&x>0,(K1,diagonal_floor)) && isfinite(noise_variance) && noise_variance>=0 ||
        throw(ArgumentError("invalid calibration constants"))
    bias_norm in (:linf,:l2) || throw(ArgumentError("invalid bias_norm"))
    A=Float64.(matrices); d=max.(diag(A[:,:,pilot_index]),diagonal_floor)
    theta=zeros(p); pilot_budgets=zeros(p)
    W=[zeros(p,p) for _ in ks]
    budgets=[zeros(p) for _ in ks]; violations=[zeros(p) for _ in ks]
    errors=fill("",length(ks))
    restart_reasons=fill("",length(ks))
    restart_counts=zeros(Int,length(ks))
    for j in 1:p
        model=_precision_model(;lp_solver)
        JuMP.@variable(model,beta[1:p])
        JuMP.@variable(model,z[1:p]>=0)
        JuMP.@expression(model,t,sum(z))
        JuMP.@constraint(model,[k=1:p],beta[k]<=z[k])
        JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
        JuMP.@constraint(model,beta[j]>=0)
        for k in 1:p
            c=K1*(max(d[k],d[j])*u[pilot_index]+noise_variance*v[pilot_index])
            r=sum(A[k,l,pilot_index]*beta[l] for l in 1:p)-(k==j)
            JuMP.@constraint(model,r<=c*beta[j]+bias[pilot_index]*t)
            JuMP.@constraint(model,-r<=c*beta[j]+bias[pilot_index]*t)
        end
        JuMP.@objective(model,Min,t)
        _solve_checked!(model)
        theta[j]=max(0,JuMP.value(beta[j])); pilot_budgets[j]=JuMP.value(t)
        JuMP.empty!(model)
    end
    for j in 1:p
        model=_precision_model(;cone=bias_norm===:l2,lp_solver)
        JuMP.@variable(model,beta[1:p])
        JuMP.@variable(model,z[1:p]>=0)
        JuMP.@expression(model,t,sum(z))
        JuMP.@constraint(model,[k=1:p],beta[k]<=z[k])
        JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
        b=bias_norm===:l2 ? JuMP.@variable(model,[1:p,1:H]) : nothing
        base=[sqrt(d[k]*theta[j])*(u[h]+noise_variance*theta[j]*v[h]) for k in 1:p,h in 1:H]
        upper=Matrix{Any}(undef,p,H); lower=similar(upper)
        for h in 1:H
            if bias_norm===:l2
                JuMP.@constraint(model,[bias[h]*t;b[:,h]] in JuMP.SecondOrderCone())
            end
            for k in 1:p
                r=sum(A[k,l,h]*beta[l] for l in 1:p)-(k==j)
                if bias_norm===:l2
                    upper[k,h]=JuMP.@constraint(model,r-b[k,h]<=ks[1]*base[k,h])
                    lower[k,h]=JuMP.@constraint(model,-r+b[k,h]<=ks[1]*base[k,h])
                else
                    upper[k,h]=JuMP.@constraint(model,r<=ks[1]*base[k,h]+bias[h]*t)
                    lower[k,h]=JuMP.@constraint(model,-r<=ks[1]*base[k,h]+bias[h]*t)
                end
            end
        end
        JuMP.@objective(model,Min,t)
        for (index,K2) in enumerate(ks)
            isempty(errors[index]) || continue
            try
                for h in 1:H,k in 1:p
                    # r includes -e_j. JuMP has moved this constant to RHS.
                    JuMP.set_normalized_rhs(upper[k,h],K2*base[k,h]+(k==j))
                    JuMP.set_normalized_rhs(lower[k,h],K2*base[k,h]-(k==j))
                end
                _solve_checked!(model)
                W[index][:,j] .= JuMP.value.(beta)
                budgets[index][j]=JuMP.value(t)
                violation=max(0,norm(W[index][:,j],1)-budgets[index][j])
                for h in 1:H
                    r=A[:,:,h]*W[index][:,j]-Float64.(1:p .== j)
                    if bias_norm===:l2
                        bh=JuMP.value.(b[:,h])
                        violation=max(violation,norm(bh)-bias[h]*budgets[index][j],
                            maximum(abs.(r-bh)-K2.*base[:,h]))
                    else
                        violation=max(violation,maximum(abs.(r)-K2.*base[:,h].-bias[h]*budgets[index][j]))
                    end
                end
                violations[index][j]=violation
            catch exception
                initial="column $j: "*sprint(showerror,exception)
                if occursin("ALMOST_OPTIMAL",initial)
                    # One fresh solve of precisely this column, retaining the
                    # same diagonal pilot, allowances and solver tolerances.
                    restart_reasons[index]*=isempty(restart_reasons[index]) ? initial : "; "*initial
                    restart_counts[index]+=1
                    try
                        fresh=_cold_aclime_column(A,j,base,bias,K2,bias_norm,lp_solver)
                        W[index][:,j] .= fresh.column
                        budgets[index][j]=fresh.budget
                        violations[index][j]=fresh.violation
                    catch cold_exception
                        errors[index]=initial*"; independent cold restart: "*sprint(showerror,cold_exception)
                    end
                else
                    errors[index]=initial
                end
            end
        end
        JuMP.empty!(model)
    end
    for index in eachindex(ks)
        if isempty(errors[index]) && maximum(violations[index])>1e-5*(1+maximum(budgets[index]))
            errors[index]="excessive precision constraint violation"
        end
    end
    (;columns=W,pilot=theta,pilot_budgets,budgets,violations,errors,restart_reasons,restart_counts)
end

"""Reuse the interval minimax-column LP across a path of raw column matrices."""
function _interval_symmetrize_path(columns;lp_solver=:highs)
    isempty(columns) && return (solutions=Matrix{Float64}[],errors=String[])
    p=size(first(columns),1)
    all(W->size(W)==(p,p)&&all(isfinite,W),columns) || throw(ArgumentError("invalid column path"))
    model=_precision_model(;lp_solver)
    JuMP.@variable(model,S[1:p,1:p],Symmetric)
    JuMP.@variable(model,d[1:p,1:p]>=0)
    JuMP.@variable(model,r>=0)
    diagonal=[JuMP.@constraint(model,S[j,j]==0) for j in 1:p]
    positive=Matrix{Any}(undef,p,p); negative=similar(positive)
    lo=Matrix{Any}(undef,p,p); hi=similar(lo)
    for j in 1:p
        JuMP.@constraint(model,sum(d[:,j])<=r)
        for k in 1:p
            positive[k,j]=JuMP.@constraint(model,d[k,j]>=S[k,j])
            negative[k,j]=JuMP.@constraint(model,d[k,j]>=-S[k,j])
        end
        for k in j+1:p
            lo[k,j]=JuMP.@constraint(model,S[k,j]>=0)
            hi[k,j]=JuMP.@constraint(model,S[k,j]<=0)
        end
    end
    JuMP.@objective(model,Min,r)
    solutions=[zeros(p,p) for _ in columns]; errors=fill("",length(columns))
    for (index,W) in enumerate(columns)
        try
            for j in 1:p
                JuMP.set_normalized_rhs(diagonal[j],W[j,j])
                for k in 1:p
                    JuMP.set_normalized_rhs(positive[k,j],-W[k,j])
                    JuMP.set_normalized_rhs(negative[k,j],W[k,j])
                end
                for k in j+1:p
                    JuMP.set_normalized_rhs(lo[k,j],min(W[k,j],W[j,k]))
                    JuMP.set_normalized_rhs(hi[k,j],max(W[k,j],W[j,k]))
                end
            end
            _solve_checked!(model)
            solutions[index] .= JuMP.value.(S)
        catch exception
            errors[index]=sprint(showerror,exception)
        end
    end
    JuMP.empty!(model)
    (;solutions,errors)
end

"""
    fit_direct_path(C, K2_values; method=:aclime, bias_constant, kwargs...)

Model-reusing K2 path for either ACLIME proposal at one fixed bias constant.
Returns `(fits, errors)`. A failed candidate has `fits[k] === nothing` and
an error string. Programs, feasibility checks, interval restrictions and
induced-norm caps match `fit_direct`; nonunique minimizers can differ.
"""
function fit_direct_path(C::DirectCovariance,K2_values;
        method::Symbol=:aclime,pilot_index::Integer=1,bias_constant::Real,
        K1::Real=2,confidence::Real=0,noise_variance::Real=1,threshold::Real=1e-6,
        lp_solver::Symbol=:highs,ntasks::Integer=1)
    method in (:aclime,:l2_aclime) || throw(ArgumentError("path supports the two ACLIME methods"))
    ks=Float64.(collect(K2_values))
    !isempty(ks) && all(x->isfinite(x)&&x>0,ks) || throw(ArgumentError("invalid K2 path"))
    isfinite(bias_constant) && bias_constant>0 || throw(ArgumentError("invalid bias constant"))
    isfinite(confidence) && confidence>=0 && isfinite(threshold) && threshold>=0 ||
        throw(ArgumentError("invalid confidence or threshold"))
    isfinite(K1) && K1>0 && isfinite(noise_variance) && noise_variance>=0 ||
        throw(ArgumentError("invalid stochastic calibration"))
    lp_solver in (:highs,:clarabel) || throw(ArgumentError("invalid lp_solver"))
    1<=ntasks<=Threads.nthreads() || throw(ArgumentError("ntasks exceeds Julia threads"))
    1<=pilot_index<=length(C.bandwidths) || throw(ArgumentError("invalid pilot_index"))
    p,_,q,H=size(C.matrices)
    ell=(3log(max(p,2))+confidence+log(H))/C.n
    domain=get(C.settings,:domain,(0.,1.))
    rates=(domain[2]-domain[1])./(C.min_visits.*C.bandwidths)
    u=sqrt(ell).*(1 .+rates).+ell.*(1 .+rates.^2)
    v=sqrt(ell).*rates.+ell.*(1 .+rates.^2)
    P=[zeros(p,p,q) for _ in ks]; W=[zeros(p,p,q) for _ in ks]
    diagnostics=[NamedTuple[] for _ in ks]; errors=fill("",length(ks))
    selected=method===:aclime ? [pilot_index] : collect(1:H)
    slices=Vector{Any}(undef,q)
    function solve_slice(time_index)
        try
            point_result=bias_aclime_path(C.matrices[:,:,time_index,selected],ks;u=u[selected],v=v[selected],
                    bias=bias_constant.*C.bandwidths[selected].^2,
                    pilot_index=findfirst(==(pilot_index),selected),K1,noise_variance,
                    bias_norm=method===:aclime ? :linf : :l2,diagonal_floor=C.n^-2,lp_solver)
            point_valid=findall(isempty,point_result.errors)
            point_symmetric=_interval_symmetrize_path(point_result.columns[point_valid];lp_solver)
            (;result=point_result,valid=point_valid,symmetric=point_symmetric,error="")
        catch exception
            (;result=nothing,valid=Int[],symmetric=nothing,error=sprint(showerror,exception))
        end
    end
    if ntasks==1
        for time_index in 1:q
            slices[time_index]=solve_slice(time_index)
        end
    else
        Threads.@threads for task in 1:ntasks
            for time_index in task:ntasks:q
                slices[time_index]=solve_slice(time_index)
            end
        end
    end
    # Assemble in original time order; threads share no JuMP models or writes
    # to diagnostic/error vectors during optimization.
    for time_index in 1:q
        slice=slices[time_index]
        if !isempty(slice.error)
            for index in eachindex(ks)
                isempty(errors[index]) && (errors[index]="time $time_index: "*slice.error)
            end
            continue
        end
        result,valid,symmetric=slice.result,slice.valid,slice.symmetric
        for index in eachindex(ks)
            if !isempty(result.errors[index]) && isempty(errors[index])
                errors[index]="time $time_index: "*result.errors[index]
            end
        end
        for (position,index) in enumerate(valid)
            if !isempty(symmetric.errors[position])
                errors[index]="time $time_index symmetry: "*symmetric.errors[position]
                continue
            end
            W[index][:,:,time_index] .= result.columns[index]
            S=symmetric.solutions[position]
            cap=min(1.0,C.n/max(opnorm(S,1),eps()))
            P[index][:,:,time_index] .= cap.*S
            push!(diagnostics[index],(;status=:optimal,max_constraint_violation=maximum(result.violations[index]),
                covariance_psd_displacement=0.,cap,pilot_minimum=minimum(result.pilot),
                cold_restart_reason=result.restart_reasons[index],cold_restart_count=result.restart_counts[index]))
        end
    end
    fits=Vector{Union{Nothing,DirectPrecisionFit}}(nothing,length(ks))
    for index in eachindex(ks)
        isempty(errors[index]) || continue
        adjacency=abs.(P[index]).>threshold
        for time_index in 1:q,j in 1:p
            adjacency[j,j,time_index]=false
        end
        positive=BitVector([isposdef(Symmetric(P[index][:,:,t])) for t in 1:q])
        settings=(;method,pilot_index,bias_constant,K1,K2=ks[index],confidence,noise_variance,
            lambda0=:quantile,psd_floor=1e-8,threshold,symmetrization=:interval_l1,
            lp_solver,socp_solver=:clarabel,lp_tolerance=lp_solver===:highs ? 1e-8 : 1e-7,
            socp_tolerance=1e-7,ell,path_reuse=:rhs_models,ntasks)
        fits[index]=DirectPrecisionFit(copy(C.grid),P[index],W[index],adjacency,positive,
            diagnostics[index],C,settings)
    end
    (;fits,errors)
end
