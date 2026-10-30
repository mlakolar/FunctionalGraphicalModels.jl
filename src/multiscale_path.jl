# Native Clarabel workspaces retain the sparse structure and symbolic analysis.
# solve! initializes new interior-point iterates; this is not iterate warm start.
mutable struct _L2Workspace
    solver::Clarabel.Solver{Float64}
    matrix::SparseArrays.SparseMatrixCSC{Float64,Int}
    rhs::Vector{Float64}
    score_rows::Matrix{Int}
end

function _l2_problem_data(A, base, bias, j, penalty, bound)
    p,_,H=size(A); n=2p+1+p*H; t=2p+1
    rows=Int[]; cols=Int[]; vals=Float64[]; rhs=Float64[]
    cones=Clarabel.SupportedCone[]
    function row!(pairs, b)
        push!(rhs,b); row=length(rhs)
        for (col,value) in pairs
            push!(rows,row);push!(cols,col);push!(vals,value)
        end
        row
    end
    for k in 1:p
        row!([(p+k,-1.)],0.)
        row!([(k,1.),(p+k,-1.)],0.)
        row!([(k,-1.),(p+k,-1.)],0.)
    end
    row!(vcat([(p+k,1.) for k in 1:p],[(t,-1.)]),0.)
    bound===nothing || row!([(t,1.)],sqrt(p)*bound)
    score_rows=zeros(Int,p,H)
    for h in 1:H,k in 1:p
        # Explicit zeros keep the CSC pattern identical between bandwidth sets.
        pairs=vcat([(l,A[k,l,h]) for l in 1:p],[(t+k+(h-1)*p,-bias[h])])
        score_rows[k,h]=row!(pairs,penalty*base[k,h]+(k==j))
        row!([(l,-v) for (l,v) in pairs],penalty*base[k,h]-(k==j))
    end
    push!(cones,Clarabel.NonnegativeConeT(length(rhs)))
    for h in 1:H
        row!([(t,-1.)],0.)
        for k in 1:p;row!([(t+k+(h-1)*p,-1.)],0.);end
        push!(cones,Clarabel.SecondOrderConeT(p+1))
        if bound!==nothing
            row!(Tuple{Int,Float64}[],bound)
            for k in 1:p;row!([(t+k+(h-1)*p,-1.)],0.);end
            push!(cones,Clarabel.SecondOrderConeT(p+1))
        end
    end
    if bound!==nothing
        row!(Tuple{Int,Float64}[],bound)
        for k in 1:p;row!([(k,-1.)],0.);end
        push!(cones,Clarabel.SecondOrderConeT(p+1))
    end
    matrix=SparseArrays.sparse(rows,cols,vals,length(rhs),n)
    objective=zeros(n);objective[t]=1.
    (;matrix,rhs,cones,objective,score_rows)
end

# Approximate infeasibility is a numerical outcome, not a certificate. Route it
# through recovery; only precise termination plus the original checks may pass.
function _l2_recoverable(status,guard_failed=false)
    status in ("PRIMAL_INFEASIBLE","DUAL_INFEASIBLE") && return false
    guard_failed || status in ("ALMOST_SOLVED","ALMOST_PRIMAL_INFEASIBLE",
        "ALMOST_DUAL_INFEASIBLE","ALMOST_INFEASIBLE","INFEASIBLE_OR_UNBOUNDED",
        "NUMERICAL_ERROR","INSUFFICIENT_PROGRESS","MAX_ITERATIONS",
        "ALMOST_OPTIMAL","SLOW_PROGRESS","ITERATION_LIMIT","OTHER_ERROR")
end

function _l2_workspace(data; tolerance=1e-8)
    n=length(data.objective)
    settings=Clarabel.Settings(verbose=false,tol_gap_abs=tolerance,tol_gap_rel=tolerance,
        tol_feas=tolerance,max_threads=1,presolve_enable=false,
        input_sparse_dropzeros=false,chordal_decomposition_enable=false)
    solver=Clarabel.Solver(SparseArrays.spzeros(n,n),data.objective,data.matrix,
        data.rhs,data.cones,settings)
    _L2Workspace(solver,data.matrix,data.rhs,data.score_rows)
end

function _l2_update!(workspace,data)
    old=workspace.matrix;new=data.matrix
    size(old)==size(new) && old.colptr==new.colptr && old.rowval==new.rowval ||
        throw(ArgumentError("workspace sparsity differs"))
    Clarabel.update_A!(workspace.solver,new)
    Clarabel.update_b!(workspace.solver,data.rhs)
    workspace.matrix=new;workspace.rhs=data.rhs;workspace.score_rows=data.score_rows
end

function _l2_penalty!(workspace,base,j,penalty)
    for h in axes(base,2),k in axes(base,1)
        row=workspace.score_rows[k,h]
        workspace.rhs[row]=penalty*base[k,h]+(k==j)
        workspace.rhs[row+1]=penalty*base[k,h]-(k==j)
    end
    Clarabel.update_b!(workspace.solver,workspace.rhs)
end

# Evaluate the original coordinate, auxiliary cone and column constraints.
# Auxiliary vectors are certificates; changing them does not change beta or t.
function _l2_primal_residual(A,base,bias,j,penalty,bound,beta,budget,g)
    p,_,H=size(A)
    violation=max(0.,norm(beta,1)-budget)
    relative=violation/(1+abs(budget))
    if bound!==nothing
        violation=max(violation,norm(beta)-bound,budget-sqrt(p)*bound)
        relative=max(relative,max(0.,norm(beta)-bound)/(1+bound),max(0.,budget-sqrt(p)*bound)/(1+sqrt(p)*bound))
    end
    for h in 1:H
        normalized_ball=max(0.,norm(g[:,h])-budget)
        relative=max(relative,normalized_ball/(1+abs(budget)))
        if bound!==nothing
            relative=max(relative,max(0.,norm(g[:,h])-bound)/(1+bound))
        end
        b=bias[h].*g[:,h]
        residual=A[:,:,h]*beta-Float64.(1:p .== j)-b
        score=max.(0.,abs.(residual).-penalty.*base[:,h])
        ball=max(0.,norm(b)-bias[h]*budget)
        violation=max(violation,maximum(score),ball)
        relative=max(relative,maximum(score./(1 .+penalty.*base[:,h])),ball/(1+abs(bias[h]*budget)))
        if bound!==nothing
            cap=max(0.,norm(b)-bias[h]*bound)
            violation=max(violation,cap);relative=max(relative,cap/(1+bias[h]*bound))
        end
    end
    (;violation,relative)
end

# Equivalent cold formulation in the original bias units. Keeping b explicit
# avoids the large normalized auxiliary variables that can stall the native path.
function _l2_unscaled_candidate(A,base,bias,j,penalty,bound;tolerance=1e-8)
    p,_,H=size(A);started=time_ns()
    model=_precision_model(;cone=true,lp_solver=:clarabel)
    for name in ("tol_gap_abs","tol_gap_rel","tol_feas")
        JuMP.set_optimizer_attribute(model,name,tolerance)
    end
    JuMP.@variable(model,beta[1:p]);JuMP.@variable(model,z[1:p]>=0)
    JuMP.@variable(model,t>=0)
    JuMP.@variable(model,b[1:p,1:H])
    JuMP.@constraint(model,[k=1:p],beta[k]<=z[k])
    JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
    JuMP.@constraint(model,sum(z)<=t)
    if bound!==nothing
        JuMP.@constraint(model,t<=sqrt(p)*bound)
        JuMP.@constraint(model,[bound;beta] in JuMP.SecondOrderCone())
    end
    for h in 1:H
        JuMP.@constraint(model,[bias[h]*t;b[:,h]] in JuMP.SecondOrderCone())
        if bound!==nothing
            JuMP.@constraint(model,[bias[h]*bound;b[:,h]] in JuMP.SecondOrderCone())
        end
        for k in 1:p
            r=JuMP.@expression(model,sum(A[k,l,h]*beta[l] for l in 1:p)-(k==j)-b[k,h] )
            JuMP.@constraint(model,r<=penalty*base[k,h])
            JuMP.@constraint(model,-r<=penalty*base[k,h])
        end
    end
    JuMP.@objective(model,Min,t)
    setup_seconds=(time_ns()-started)/1e9;tick=time_ns()
    JuMP.optimize!(model)
    solve_seconds=(time_ns()-tick)/1e9
    raw_status=string(JuMP.termination_status(model))
    termination=raw_status=="INFEASIBLE" ? "PRIMAL_INFEASIBLE" : raw_status
    iterations_available=true
    iterations=try
        Int(JuMP.barrier_iterations(model))
    catch
        iterations_available=false;0
    end
    values=JuMP.has_values(model)
    coefficient=values ? JuMP.value.(beta) : fill(NaN,p)
    budget=values ? JuMP.value(t) : NaN
    witness=values ? JuMP.value.(b) : fill(NaN,p,H)
    g=zeros(p,H)
    for h in 1:H
        # At zero bias the exact witness is zero. The unchanged score guard
        # must still pass after discarding any numerical residue in b.
        bias[h]>0 && (g[:,h].=witness[:,h]./bias[h])
    end
    JuMP.empty!(model)
    (;beta=coefficient,budget,g,termination,raw_status,iterations,iterations_available,
        setup_seconds,solve_seconds)
end

# A precise infeasibility certificate concerns only A*x+s=b, s in K, so
# removing the objective cannot change its validity for the original problem.
# A feasible probe is deliberately NOT returned as an optimized column.
function _l2_infeasibility_probe(A,base,bias,j,penalty,bound;tolerance=1e-8)
    p,_,H=size(A);started=time_ns()
    original=_l2_problem_data(A,base,bias,j,penalty,bound)
    data=(;original...,objective=zeros(length(original.objective)))
    workspace=_l2_workspace(data;tolerance)
    setup_seconds=(time_ns()-started)/1e9;tick=time_ns()
    solution=Clarabel.solve!(workspace.solver)
    solve_seconds=(time_ns()-tick)/1e9;raw_status=string(solution.status)
    termination=raw_status=="PRIMAL_INFEASIBLE" ? raw_status : "INFEASIBILITY_PROBE_"*raw_status
    (;beta=copy(solution.x[1:p]),budget=solution.x[2p+1],
        g=copy(reshape(solution.x[2p+2:end],p,H)),termination,raw_status,
        iterations=Int(solution.iterations),iterations_available=true,setup_seconds,solve_seconds,
        precision_estimate=false)
end

function _l2_pilot(A,d,u,v,bias,pilot_index,K1,noise_variance,lp_solver)
    p=size(A,1);theta=zeros(p)
    for j in 1:p
        model=_precision_model(;lp_solver)
        JuMP.@variable(model,beta[1:p]);JuMP.@variable(model,z[1:p]>=0)
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
        JuMP.@objective(model,Min,t);_solve_checked!(model)
        theta[j]=max(0.,JuMP.value(beta[j]));JuMP.empty!(model)
    end
    theta
end

"""
    multiscale_l2_aclime_path(problems, penalties; calibration=:pilot, ...)

Solve an explicit sequence of multiscale covariance problems, crossing each
with the supplied penalty path in caller order. Each problem is a named tuple
with `matrices` (p×p×H), `u`, `v`, `bias`. For `:pilot`, optional `pilot_index`
(default 1) selects the precision pilot. For `:spectral`, supply positive `d`
and `kappa_bound`; the column bounds are `2*kappa_bound ./ d`.

A native Clarabel workspace per column is reused across penalties and across
problems with the same p, H and calibration. Changing H builds new workspaces.
No iterate warm start is claimed. Inputs may be indefinite. No tuning,
clipping, thresholding, symmetrization or precision norm cap occurs here.
`results[a][b]` corresponds to problem a and penalty b. Numerical failures
remain NaN and have explicit column records; certified infeasibility is
reported as failure rather than substituted by zero. If the score allowances
contain the unit vector at every width, the exact optimum is returned as zero
without building or solving a conic program (`termination="ANALYTIC_ZERO"`).
Numerical rejections use cold native and equivalent unscaled JuMP/Clarabel
recovery at the same or stricter tolerances, with all attempts retained.
Use `interval_symmetrize` on successful columns when constructing the final precision estimator.
"""
function multiscale_l2_aclime_path(problems,penalties;calibration::Symbol=:pilot,
        K1::Real=2,noise_variance::Real=1,diagonal_floor::Real=1e-8,
        lp_solver::Symbol=:highs,reuse::Bool=true)
    calibration in (:pilot,:spectral) || throw(ArgumentError("invalid calibration"))
    ks=Float64.(collect(penalties));specs=collect(problems)
    !isempty(ks) && all(x->isfinite(x)&&x>0,ks) && allunique(ks) || throw(ArgumentError("invalid penalties"))
    !isempty(specs) || throw(ArgumentError("empty problem sequence"))
    isfinite(K1)&&K1>0 && isfinite(diagonal_floor)&&diagonal_floor>0 &&
        isfinite(noise_variance)&&noise_variance>=0 || throw(ArgumentError("invalid calibration constants"))
    lp_solver in (:highs,:clarabel) || throw(ArgumentError("invalid LP solver"))
    # Validate the entire sequence before fitting any member.
    for spec in specs
        A=spec.matrices;p,p2,H=size(A)
        p==p2 && p>0 && H>0 && all(isfinite,A) &&
            all(issymmetric(A[:,:,h]) for h in 1:H) || throw(ArgumentError("invalid covariance sequence"))
        length(spec.u)==length(spec.v)==length(spec.bias)==H || throw(DimensionMismatch("resource lengths"))
        all(x->isfinite(x)&&x>=0,spec.u) && all(x->isfinite(x)&&x>=0,spec.v) &&
            all(x->isfinite(x)&&(calibration===:pilot ? x>0 : x>=0),spec.bias) || throw(ArgumentError("invalid resource/bias values"))
        if calibration===:pilot
            1<=get(spec,:pilot_index,1)<=H || throw(ArgumentError("invalid pilot index"))
        else
            length(spec.d)==p && all(x->isfinite(x)&&x>0,spec.d) &&
                isfinite(spec.kappa_bound)&&spec.kappa_bound>=1 || throw(ArgumentError("invalid spectral bounds"))
        end
    end
    workspaces=Any[];previous_shape=nothing;results=Vector{Any}()
    builds=0;matrix_updates=0;rhs_updates=0;analytic_zeros=0;recovery_attempts=0;cold_restarts=0;unscaled_solves=0;infeasibility_probes=0
    for spec in specs
        A=Float64.(spec.matrices);p,_,H=size(A);pilot_index=get(spec,:pilot_index,1)
        d=calibration===:pilot ? max.(diag(A[:,:,pilot_index]),diagonal_floor) : Float64.(spec.d)
        theta=try
            calibration===:pilot ? _l2_pilot(A,d,spec.u,spec.v,spec.bias,pilot_index,K1,noise_variance,lp_solver) :
                2spec.kappa_bound./d
        catch exception
            message="pilot: "*sprint(showerror,exception)
            push!(results,[(;penalty=k,columns=fill(NaN,p,p),budgets=fill(NaN,p),pilot=fill(NaN,p),
                column_records=[(;column=j,status=:failure,termination="PILOT_FAILURE",iterations=0,
                    setup_seconds=0.,solve_seconds=0.,max_violation=NaN,max_relative_violation=NaN,retry_count=0,backend=:pilot,iterations_available=false,attempt_records=NamedTuple[],error=message) for j in 1:p],
                status=:failure) for k in ks])
            workspaces=Any[];previous_shape=nothing
            continue
        end
        bases=[Float64[sqrt(d[k]*theta[j])*(spec.u[h]+noise_variance*theta[j]*spec.v[h]) for k in 1:p,h in 1:H] for j in 1:p]
        all(b->all(isfinite,b),bases) && all(isfinite,theta) || throw(ArgumentError("calibration overflow"))
        shape=(p,H)
        if previous_shape!=shape || !reuse
            workspaces=Any[nothing for _ in 1:p];previous_shape=shape
        end
        columns=[fill(NaN,p,p) for _ in ks];budgets=[fill(NaN,p) for _ in ks]
        records=[NamedTuple[] for _ in ks]
        for j in 1:p
            base=bases[j];bound=calibration===:spectral ? theta[j] : nothing
            data_current=false # skipped analytic points must not hide a new covariance problem
            for (index,penalty) in enumerate(ks)
                # beta=g=t=0 is feasible iff every diagonal score allowance is >=1.
                # Since t>=norm(beta,1)>=0, this is a global optimum, not a fallback.
                if all(penalty*base[j,h]>=1. for h in 1:H)
                    columns[index][:,j].=0.;budgets[index][j]=0.;analytic_zeros+=1
                    push!(records[index],(;column=j,status=:optimal,termination="ANALYTIC_ZERO",iterations=0,
                        setup_seconds=0.,solve_seconds=0.,max_violation=0.,max_relative_violation=0.,retry_count=0,backend=:analytic,iterations_available=true,attempt_records=NamedTuple[],error=""))
                    continue
                end
                setup_seconds=0.;solve_seconds=0.;attempt_records=NamedTuple[]
                initially_reused=reuse && workspaces[j]!==nothing
                steps=[(:native,1e-8)]
                initially_reused && push!(steps,(:native,1e-8))
                push!(steps,(:jump_unscaled,1e-8))
                bound!==nothing && push!(steps,(:native_certificate,1e-8))
                append!(steps,[(:jump_unscaled,1e-9),(:jump_unscaled,1e-10)])
                for (attempt,(backend,tolerance)) in enumerate(steps)
                    attempt>1 && (recovery_attempts+=1)
                    started=time_ns();status="NOT_RUN";raw_status=status;iterations=0;iterations_available=true
                    violation=relative=NaN;guard_failed=false;setup_this=solve_this=0.
                    candidate=nothing;mode=:fresh
                    try
                        if backend===:native
                            if attempt>1 || !initially_reused
                                data=_l2_problem_data(A,base,spec.bias,j,penalty,bound)
                                workspaces[j]=_l2_workspace(data;tolerance);builds+=1
                                attempt>1 && (cold_restarts+=1)
                            elseif !data_current
                                data=_l2_problem_data(A,base,spec.bias,j,penalty,bound)
                                _l2_update!(workspaces[j],data);matrix_updates+=1;mode=:matrix_update
                            else
                                _l2_penalty!(workspaces[j],base,j,penalty);rhs_updates+=1;mode=:rhs_update
                            end
                            data_current=true
                            setup_this=(time_ns()-started)/1e9;setup_seconds+=setup_this;tick=time_ns()
                            solution=Clarabel.solve!(workspaces[j].solver)
                            solve_this=(time_ns()-tick)/1e9;solve_seconds+=solve_this
                            status=raw_status=string(solution.status);iterations=Int(solution.iterations)
                            beta=copy(solution.x[1:p]);budget=solution.x[2p+1]
                            g=copy(reshape(solution.x[2p+2:end],p,H))
                        else
                            if backend===:native_certificate
                                infeasibility_probes+=1
                                cold=_l2_infeasibility_probe(A,base,spec.bias,j,penalty,bound;tolerance)
                            else
                                unscaled_solves+=1
                                cold=_l2_unscaled_candidate(A,base,spec.bias,j,penalty,bound;tolerance)
                            end
                            setup_this=cold.setup_seconds;solve_this=cold.solve_seconds
                            setup_seconds+=setup_this;solve_seconds+=solve_this
                            status=cold.termination;raw_status=cold.raw_status
                            iterations=cold.iterations;iterations_available=cold.iterations_available
                            beta=cold.beta;budget=cold.budget;g=cold.g
                        end
                        candidate=(;beta,budget,g)
                        status in ("SOLVED","OPTIMAL") || error("$backend status $status")
                        check=_l2_primal_residual(A,base,spec.bias,j,penalty,bound,beta,budget,g)
                        violation=check.violation;relative=check.relative
                        guard_failed=!(all(isfinite,beta)&&isfinite(budget)&&isfinite(relative)&&relative<=1e-5)
                        guard_failed && error("primal residual check failed: $relative")
                        push!(attempt_records,(;backend,mode,termination=status,raw_status,iterations,iterations_available,tolerance,
                            setup_seconds=setup_this,solve_seconds=solve_this,max_violation=violation,
                            max_relative_violation=relative,rejected_candidate=nothing,error=""))
                        columns[index][:,j].=beta;budgets[index][j]=budget
                        push!(records[index],(;column=j,status=:optimal,termination=status,backend,
                            iterations=sum(r.iterations for r in attempt_records),
                            iterations_available=all(r.iterations_available for r in attempt_records),setup_seconds,solve_seconds,
                            max_violation=violation,max_relative_violation=relative,retry_count=attempt-1,attempt_records,error=""))
                        break
                    catch exception
                        message=sprint(showerror,exception)
                        push!(attempt_records,(;backend,mode,termination=status,raw_status,iterations,iterations_available,tolerance,
                            setup_seconds=setup_this,solve_seconds=solve_this,max_violation=violation,
                            max_relative_violation=relative,rejected_candidate=status=="PRIMAL_INFEASIBLE" ? nothing : candidate,error=message))
                        workspaces[j]=nothing
                        recoverable=_l2_recoverable(status,guard_failed) ||
                            (backend===:native_certificate && status!="PRIMAL_INFEASIBLE")
                        if recoverable && attempt<length(steps)
                            continue
                        end
                        push!(records[index],(;column=j,status=:failure,termination=status,backend,
                            iterations=sum(r.iterations for r in attempt_records),
                            iterations_available=all(r.iterations_available for r in attempt_records),setup_seconds,solve_seconds,
                            max_violation=violation,max_relative_violation=relative,retry_count=attempt-1,attempt_records,error=message))
                        break
                    end
                end
            end
        end
        push!(results,[(;penalty=ks[i],columns=columns[i],budgets=budgets[i],pilot=copy(theta),
            column_records=records[i],status=all(r->r.status===:optimal,records[i]) ? :optimal : :failure) for i in eachindex(ks)])
    end
    (;results,diagnostics=(;workspace_builds=builds,matrix_updates,rhs_updates,
        reuse,iterate_warm_start=false,calibration,analytic_zeros,recovery_attempts,cold_restarts,unscaled_solves,infeasibility_probes))
end

"""
    fit_multiscale_l2_path(C, penalties; bandwidth_sets=[eachindex(C.bandwidths)],
                          bias_constants, pilot_indices=first.(bandwidth_sets), ...)

Fit pilot-based multiscale ℓ₂-ACLIME over a fixed Cartesian grid of bandwidth
index sets, bias constants, and score penalties. All widths within one set are
intersected in ONE program; singleton sets give separate one-width programs.
`pilot_indices` refers to indices in C, one per set, and each must belong to its
set. Confidence resources are recalculated for each set's size. A precision
pilot is solved once per time/set/bias constant, independent of the penalty.
Native conic workspaces are reused between equal-size sets. Symmetrization
uses the existing interval LP path, followed by the same induced-norm cap n as
`fit_direct`. No data-dependent selection is performed.

Return `fits[configuration, penalty]` (nothing on any failed time slice),
matching `errors`, explicit `configurations`, and per-time `solver_diagnostics`.
Spectral-calibrated programs are available through `multiscale_l2_aclime_path`
with explicit spectral calibration problems; this wrapper uses precision pilots.
"""
function fit_multiscale_l2_path(C::DirectCovariance,penalties;
        bandwidth_sets=[collect(eachindex(C.bandwidths))],bias_constants,
        pilot_indices=nothing,K1::Real=2,confidence::Real=0,noise_variance::Real=1,
        threshold::Real=1e-6,lp_solver::Symbol=:highs,reuse::Bool=true)
    sets=[Int.(collect(s)) for s in bandwidth_sets];ks=Float64.(collect(penalties))
    bs=Float64.(collect(bias_constants));p,_,q,H=size(C.matrices)
    !isempty(sets)&&all(s->!isempty(s)&&allunique(s)&&all(i->1<=i<=H,s),sets) || throw(ArgumentError("invalid bandwidth sets"))
    pilots=pilot_indices===nothing ? first.(sets) : Int.(collect(pilot_indices))
    length(pilots)==length(sets)&&all(pilots[i] in sets[i] for i in eachindex(sets)) || throw(ArgumentError("invalid pilot indices"))
    !isempty(bs)&&all(x->isfinite(x)&&x>0,bs)&&allunique(bs) || throw(ArgumentError("invalid bias constants"))
    !isempty(ks)&&all(x->isfinite(x)&&x>0,ks)&&allunique(ks) || throw(ArgumentError("invalid penalties"))
    isfinite(confidence)&&confidence>=0&&isfinite(threshold)&&threshold>=0 || throw(ArgumentError("invalid confidence/threshold"))
    configurations=[(;bandwidth_indices=copy(sets[a]),bandwidths=copy(C.bandwidths[sets[a]]),
        pilot_index=pilots[a],bias_constant=b) for b in bs for a in eachindex(sets)]
    count=length(configurations);P=[zeros(p,p,q) for _ in 1:count,_ in ks]
    W=[zeros(p,p,q) for _ in 1:count,_ in ks]
    records=[NamedTuple[] for _ in 1:count,_ in ks];errors=fill("",count,length(ks))
    solver_diagnostics=NamedTuple[]
    domain=get(C.settings,:domain,(0.,1.));domain_length=domain[2]-domain[1]
    for time_index in 1:q
        problems=map(configurations) do cfg
            indices=cfg.bandwidth_indices;hs=C.bandwidths[indices]
            ell=(3log(max(p,2))+confidence+log(length(indices)))/C.n
            rates=domain_length./(C.min_visits.*hs)
            u=sqrt(ell).*(1 .+rates).+ell.*(1 .+rates.^2)
            v=sqrt(ell).*rates.+ell.*(1 .+rates.^2)
            (;matrices=C.matrices[:,:,time_index,indices],u,v,bias=cfg.bias_constant.*hs.^2,
                pilot_index=findfirst(==(cfg.pilot_index),indices))
        end
        path=multiscale_l2_aclime_path(problems,ks;K1,noise_variance,diagonal_floor=C.n^-2,lp_solver,reuse)
        push!(solver_diagnostics,path.diagnostics)
        valid=[(a,b) for a in 1:count for b in eachindex(ks) if path.results[a][b].status===:optimal]
        symmetric=_interval_symmetrize_path([path.results[a][b].columns for (a,b) in valid];lp_solver)
        for a in 1:count,b in eachindex(ks)
            r=path.results[a][b]
            if r.status!==:optimal
                errors[a,b]="time $time_index: "*join([x.error for x in r.column_records if !isempty(x.error)],"; ")
            end
        end
        for (position,(a,b)) in enumerate(valid)
            if !isempty(symmetric.errors[position])
                errors[a,b]="time $time_index symmetrization: "*symmetric.errors[position];continue
            end
            r=path.results[a][b];S=symmetric.solutions[position]
            cap=min(1.,C.n/max(opnorm(S,1),eps()))
            P[a,b][:,:,time_index].=cap.*S;W[a,b][:,:,time_index].=r.columns
            push!(records[a,b],(;status=:optimal,cap,column_records=r.column_records,
                max_constraint_violation=maximum(x.max_violation for x in r.column_records),covariance_psd_displacement=0.))
        end
    end
    fits=Matrix{Union{Nothing,DirectPrecisionFit}}(nothing,count,length(ks))
    for a in 1:count,b in eachindex(ks)
        isempty(errors[a,b]) || continue
        adjacency=abs.(P[a,b]).>threshold
        for t in 1:q,j in 1:p;adjacency[j,j,t]=false;end
        positive=BitVector([isposdef(Symmetric(P[a,b][:,:,t])) for t in 1:q])
        settings=(;method=:l2_aclime,configurations[a]...,K1,K2=ks[b],confidence,noise_variance,
            threshold,lp_solver,socp_solver=:clarabel,path_reuse=:native_clarabel,
            iterate_warm_start=false,numerical_recovery=:cold_native_then_unscaled_jump,
            native_tolerance=1e-8,symmetrization=:interval_l1)
        fits[a,b]=DirectPrecisionFit(copy(C.grid),P[a,b],W[a,b],adjacency,positive,records[a,b],C,settings)
    end
    (;fits,errors,configurations,penalties=ks,solver_diagnostics)
end
