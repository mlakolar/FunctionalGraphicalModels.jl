"""Direct, noise-diagonal-free covariance estimates: node × node × time × bandwidth.
`settings` records the design, centering and boundary conventions. These matrices
can be indefinite. `pair_subjects` counts subjects contributing a diagonal pair.
"""
struct DirectCovariance
    grid::Vector{Float64}
    bandwidths::Vector{Float64}
    matrices::Array{Float64,4}
    n::Int
    min_visits::Int
    pair_subjects::Array{Int,3}
    settings::NamedTuple
end

"""Direct pointwise precision fits, raw columns and numerical diagnostics.
No precision PSD projection is applied. The covariance projection used by scaled
lasso is recorded separately in `diagnostics`.
"""
struct DirectPrecisionFit
    grid::Vector{Float64}
    precisions::Array{Float64,3}
    columns::Array{Float64,3}
    adjacency::BitArray{3}
    positive_definite::BitVector
    diagnostics::Vector{NamedTuple}
    covariance::DirectCovariance
    settings::NamedTuple
end

# Stable binomial inclusion probabilities; recurrence avoids cancellation for
# tiny windows and N>=2, unlike 1-(1-a)^m-m*a*(1-a)^(m-1).
function _inclusion(m, a)
    a == 1 && return (1.0, m >= 2 ? 1.0 : 0.0)
    p0, p1, p2 = 1.0, 0.0, 0.0
    for _ in 1:m
        p2 += a*p1
        p1 = (1-a)*p1 + a*p0
        p0 *= 1-a
    end
    (-expm1(m*log1p(-a)), p2)
end

function _equivalent_kernel(t, h, domain, boundary)
    a, b = domain
    lo, hi = max(-1.0, (a-t)/h), min(1.0, (b-t)/h)
    boundary === :error && (lo > -1 || hi < 1) && throw(ArgumentError(
        "bandwidth crosses the domain boundary; use interior targets or boundary=:local_linear"))
    k = KernelSmoothing.EpanechnikovKernel()
    M0(u) = .75*(u-u^3/3)
    M1(u) = .75*(u^2/2-u^4/4)
    M2(u) = .75*(u^3/3-u^5/5)
    s0, s1, s2 = M0(hi)-M0(lo), M1(hi)-M1(lo), M2(hi)-M2(lo)
    det = s0*s2-s1^2
    det > 0 || throw(ArgumentError("degenerate smoothing window"))
    # The factor hi-lo changes a conditional-uniform window average into
    # the equivalent-kernel integral. At interior points this is exactly 2K.
    weight(u) = (hi-lo)*k(u)*(s2-s1*u)/det
    weight, h*(hi-lo)/(b-a)
end

"""
    direct_covariance(data, grid; bandwidths, design=:independent,
                      mean=(j,t)->0.0, boundary=:error)

Count-normalized product-kernel covariance from raw FDA observations. Times must
be independent uniform draws on the shared declared domain. `:independent`
means independent schedules across subjects AND nodes. Off-diagonal entries use
products of inclusion-corrected within-window averages; diagonal entries exclude
repeated visits, removing independent measurement-noise variance without knowing
it. Every curve needs at least two visits. Unequal visit counts use their own
exact binomial inclusion probabilities (the theory's common-m bound uses the
minimum count). The supplied mean is evaluated at original times; zero means
known centering, not empirical centering. Estimate means on independent subjects
and pass their evaluation function when the mean is unknown.

`:synchronous` requires one complete shared random grid and uses distinct-visit
pairs for every entry. Its expectation has the same smoothed target, but the
independent-schedule concentration theorem does not apply. It also needs at least
two common grid points inside each window; empty windows remain zero and are
reported through `pair_subjects`. A fixed deterministic grid is not the uniform
sampling model. Stratified grids in the paper are an explicitly empirical use.

`boundary=:local_linear` uses the analytically normalized boundary equivalent
Epanechnikov kernel, which can be signed. Interior points recover the report's
estimator exactly. Boundary fits are an extension; the interior stochastic
theorem is not asserted for them. Inputs are copied and no PSD repair is made.
"""
function direct_covariance(data::FDA.MultivariateFunctionalData, grid;
        bandwidths, design::Symbol=:independent, mean=(j,t)->0.0,
        boundary::Symbol=:error)
    design in (:independent, :synchronous) || throw(ArgumentError("unknown design"))
    boundary in (:error, :local_linear) || throw(ArgumentError("unknown boundary rule"))
    g, hs = Float64.(collect(grid)), Float64.(collect(bandwidths))
    !isempty(g) && all(isfinite,g) && issorted(g) && all(diff(g).>0) ||
        throw(ArgumentError("grid must be finite, nonempty and strictly increasing"))
    !isempty(hs) && all(x->isfinite(x)&&x>0,hs) && allunique(hs) ||
        throw(ArgumentError("bandwidths must be positive, finite and unique"))
    nodes = data.components
    all(fd->fd isa Union{FDA.GridFunctionalData,FDA.IrregularFunctionalData},nodes) ||
        throw(ArgumentError("direct covariance requires raw observations"))
    n, p = FDA.nobs(data), FDA.ncomponents(data)
    n >= 2 || throw(ArgumentError("at least two subjects are required"))
    dom = FDA.domain(first(nodes))
    dom[1] < dom[2] && all(fd->FDA.domain(fd)==dom,nodes) ||
        throw(ArgumentError("nodes must declare the same nondegenerate sampling domain"))
    all(t->dom[1]<=t<=dom[2],g) || throw(ArgumentError("targets outside sampling domain"))
    observations = [[fd[i] for i in 1:n] for fd in nodes]
    mmin = minimum(length(observations[j][i][1]) for j in 1:p,i in 1:n)
    mmin >= 2 || throw(ArgumentError("every curve needs at least two visits"))
    if design === :synchronous
        times = observations[1][1][1]
        all(observations[j][i][1]==times for j in 1:p,i in 1:n) ||
            throw(ArgumentError("synchronous design requires exactly one shared grid"))
    end
    # Cache residuals once; caller-supplied means may be expensive.
    residuals = [[observations[j][i][2] .- mean.(j,observations[j][i][1])
                  for i in 1:n] for j in 1:p]
    all(all(isfinite,y) for node in residuals for y in node) ||
        throw(ArgumentError("centered observations must be finite"))
    A = zeros(p,p,length(g),length(hs))
    counts = zeros(Int,p,length(g),length(hs))
    for (l,h) in enumerate(hs), (v,t) in enumerate(g)
        weight, prob = _equivalent_kernel(t,h,dom,boundary)
        sums, pairs = zeros(p), zeros(p)
        for i in 1:n
            R = Vector{Vector{Float64}}(undef,p)
            N = zeros(Int,p)
            pi1, pi2 = zeros(p), zeros(p)
            for j in 1:p
                times = observations[j][i][1]
                idx = findall(u->abs(u-t)<=h,times)
                N[j] = length(idx)
                pi1[j],pi2[j] = _inclusion(length(times),prob)
                R[j] = [weight((times[r]-t)/h)*residuals[j][i][r] for r in idx]
                sums[j] = sum(R[j])
                pairs[j] = N[j]>=2 ? (sums[j]^2-sum(abs2,R[j])) /
                    (pi2[j]*N[j]*(N[j]-1)) : 0.0
                counts[j,v,l] += N[j]>=2
            end
            if design === :independent
                Q = [N[j]>0 ? sums[j]/(pi1[j]*N[j]) : 0.0 for j in 1:p]
                A[:,:,v,l] .+= Q*transpose(Q)
                for j in 1:p
                    A[j,j,v,l] += pairs[j]-Q[j]^2
                end
            elseif N[1]>=2
                for j in 1:p,k in 1:p
                    A[j,k,v,l] += (sums[j]*sums[k]-dot(R[j],R[k])) /
                        (pi2[1]*N[1]*(N[1]-1))
                end
            end
        end
    end
    A ./= n
    DirectCovariance(g,hs,A,n,mmin,counts,
        (;design,boundary,centering=:supplied_mean,kernel=:epanechnikov,domain=dom))
end

function _conic_model()
    model = JuMP.Model(Clarabel.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model,"tol_gap_abs",1e-7)
    JuMP.set_optimizer_attribute(model,"tol_gap_rel",1e-7)
    JuMP.set_optimizer_attribute(model,"tol_feas",1e-7)
    model
end

function _solve_checked!(model)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    status == JuMP.MOI.OPTIMAL || error("precision program failed: $status")
    JuMP.has_values(model) || error("precision program has no primal solution")
    string(status)
end

"""Symmetric minimax-column-ℓ₁ LP with pairwise interval constraints.
This is the corrected common symmetrization for induced and Frobenius losses.
The diagonal is retained; positive definiteness is not imposed.
"""
function interval_symmetrize(W::AbstractMatrix{<:Real})
    p = size(W,1)
    p>0 && size(W,2)==p && all(isfinite,W) || throw(ArgumentError("W must be finite and square"))
    model = _conic_model()
    JuMP.@variable(model,S[1:p,1:p],Symmetric)
    JuMP.@variable(model,d[1:p,1:p]>=0)
    JuMP.@variable(model,r>=0)
    for j in 1:p
        JuMP.@constraint(model,S[j,j]==W[j,j])
        JuMP.@constraint(model,sum(d[:,j])<=r)
        for k in 1:p
            JuMP.@constraint(model,d[k,j]>=S[k,j]-W[k,j])
            JuMP.@constraint(model,d[k,j]>=W[k,j]-S[k,j])
        end
        for k in j+1:p
            JuMP.@constraint(model,S[k,j]>=min(W[k,j],W[j,k]))
            JuMP.@constraint(model,S[k,j]<=max(W[k,j],W[j,k]))
        end
    end
    JuMP.@objective(model,Min,r)
    _solve_checked!(model)
    Matrix(JuMP.value.(S))
end

"""
    bias_aclime(matrices; u, v, bias, pilot_index=1, K1=2, K2=2,
                noise_variance=1, bias_norm=:linf, diagonal_floor=1e-8)

Two-stage bias-aware ACLIME. `matrices` is p×p×H. The pilot uses one bandwidth;
stage two intersects all H constraints. With `:linf`, the residual allowance is
`K2*sqrt(d[k]*theta[j])*(u[h]+noise_variance*theta[j]*v[h])+bias[h]*t`.
With `:l2`, a separate vector b[h] obeys `norm(b[h])<=bias[h]*t`, while
the stochastic residual obeys the coordinatewise allowance. Both minimize the
epigraph t of the column ℓ₁ norm. **This is a bias budget, not a ridge penalty.**
Positive bias envelopes ensure an all-outcome feasible witness without imposing
the previously incorrect O(n²) budget. Matrices need not be PSD. Settings are
explicit calibration inputs; finite multipliers alone do not certify coverage.
Returns raw columns, diagonal pilots, budgets and checked primal violations.
"""
function bias_aclime(matrices::AbstractArray{<:Real,3}; u, v, bias,
        pilot_index::Integer=1, K1::Real=2,K2::Real=2,noise_variance::Real=1,
        bias_norm::Symbol=:linf,diagonal_floor::Real=1e-8)
    p,p2,H = size(matrices)
    p==p2 && p>0 && H>0 && all(isfinite,matrices) || throw(ArgumentError("invalid covariance matrices"))
    all(issymmetric(matrices[:,:,h]) for h in 1:H) || throw(ArgumentError("matrices must be symmetric"))
    1<=pilot_index<=H || throw(ArgumentError("invalid pilot_index"))
    length(u)==length(v)==length(bias)==H || throw(DimensionMismatch("one resource/bias per bandwidth"))
    all(x->isfinite(x)&&x>=0,u) && all(x->isfinite(x)&&x>=0,v) &&
        all(x->isfinite(x)&&x>0,bias) || throw(ArgumentError("invalid resources or bias envelopes"))
    all(x->isfinite(x)&&x>0,(K1,K2,diagonal_floor)) && isfinite(noise_variance) && noise_variance>=0 ||
        throw(ArgumentError("invalid calibration constants"))
    bias_norm in (:linf,:l2) || throw(ArgumentError("bias_norm must be :linf or :l2"))
    A = Float64.(matrices)
    d = max.(diag(A[:,:,pilot_index]),diagonal_floor)
    theta, budgets, W, violations = zeros(p),zeros(p),zeros(p,p),zeros(p)
    for j in 1:p
        model = _conic_model()
        JuMP.@variable(model,beta[1:p])
        JuMP.@variable(model,z[1:p]>=0)
        JuMP.@expression(model,t,sum(z))
        JuMP.@constraint(model,[k=1:p], beta[k]<=z[k])
        JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
        JuMP.@constraint(model,beta[j]>=0)
        for k in 1:p
            c = K1*(max(d[k],d[j])*u[pilot_index]+noise_variance*v[pilot_index])
            r = sum(A[k,l,pilot_index]*beta[l] for l in 1:p)-(k==j)
            JuMP.@constraint(model,r<=c*beta[j]+bias[pilot_index]*t)
            JuMP.@constraint(model,-r<=c*beta[j]+bias[pilot_index]*t)
        end
        JuMP.@objective(model,Min,t)
        _solve_checked!(model)
        theta[j] = max(0,JuMP.value(beta[j]))
    end
    for j in 1:p
        model = _conic_model()
        JuMP.@variable(model,beta[1:p])
        JuMP.@variable(model,z[1:p]>=0)
        JuMP.@expression(model,t,sum(z))
        JuMP.@constraint(model,[k=1:p], beta[k]<=z[k])
        JuMP.@constraint(model,[k=1:p],-beta[k]<=z[k])
        b = bias_norm===:l2 ? JuMP.@variable(model,[1:p,1:H]) : nothing
        allowances = [K2*sqrt(d[k]*theta[j])*(u[h]+noise_variance*theta[j]*v[h]) for k in 1:p,h in 1:H]
        for h in 1:H
            if bias_norm===:l2
                JuMP.@constraint(model,[bias[h]*t; b[:,h]] in JuMP.SecondOrderCone())
            end
            for k in 1:p
                r = sum(A[k,l,h]*beta[l] for l in 1:p)-(k==j)
                if bias_norm===:l2
                    JuMP.@constraint(model,r-b[k,h]<=allowances[k,h])
                    JuMP.@constraint(model,-r+b[k,h]<=allowances[k,h])
                else
                    JuMP.@constraint(model,r<=allowances[k,h]+bias[h]*t)
                    JuMP.@constraint(model,-r<=allowances[k,h]+bias[h]*t)
                end
            end
        end
        JuMP.@objective(model,Min,t)
        _solve_checked!(model)
        W[:,j] .= JuMP.value.(beta)
        budgets[j] = JuMP.value(t)
        violations[j] = max(0,norm(W[:,j],1)-budgets[j])
        for h in 1:H
            r = A[:,:,h]*W[:,j]-Float64.(1:p .== j)
            if bias_norm===:l2
                bh = JuMP.value.(b[:,h])
                violations[j] = max(violations[j],norm(bh)-bias[h]*budgets[j],
                    maximum(abs.(r-bh)-allowances[:,h]))
            else
                violations[j] = max(violations[j],maximum(abs.(r)-allowances[:,h].-bias[h]*budgets[j]))
            end
        end
    end
    maximum(violations)<=1e-5*(1+maximum(budgets)) || error("excessive precision constraint violation")
    (;columns=W,pilot=theta,budgets,violations)
end

"""
    fit_direct(covariance; method=:aclime, pilot_index=1, bias_constant,
               K1=2, K2=2, confidence=0, noise_variance=1,
               lambda0=:quantile, psd_floor=1e-8, threshold=1e-6)
    fit_direct(data, grid; bandwidths, bias_constant, kwargs...)

`:aclime` uses the pilot bandwidth and coordinatewise smoothing-bias budget;
`:l2_aclime` intersects all supplied bandwidths with ℓ₂ bias budgets. Each uses
the report's u,v resources with minimum visit count. `bias_constant*h²` must
bound the relevant covariance bias (entrywise or operator, respectively).
Constants/smoothness are supplied, never estimated from truth by this routine.

`:scaled_lasso` applies CovarianceSelection.slassoinv to the pilot covariance,
after an explicit eigenvalue projection with positive floor. The projection's
operator displacement is returned. `n` counts subjects, never visits. This is
an empirical adaptation to a smoothed covariance, not the iid theorem of Sun
and Zhang. Raw columns use their noise/residual-variance normalization; the
corrected interval LP provides a common symmetrization for all three methods.
The final induced-ℓ₁ cap is n. Neither positive definiteness nor graph selection
consistency is guaranteed. All solver failures throw and must be recorded by
experiment drivers; no failed fit is replaced by a reference graph.
"""
function fit_direct(C::DirectCovariance;method::Symbol=:aclime,pilot_index::Integer=1,
        bias_constant::Real,K1::Real=2,K2::Real=2,confidence::Real=0,
        noise_variance::Real=1,lambda0=:quantile,psd_floor::Real=1e-8,threshold::Real=1e-6)
    method in (:aclime,:l2_aclime,:scaled_lasso) || throw(ArgumentError("unknown precision method"))
    isfinite(bias_constant) && bias_constant>0 || throw(ArgumentError("bias_constant must be positive"))
    isfinite(confidence) && confidence>=0 || throw(ArgumentError("confidence must be nonnegative"))
    isfinite(threshold) && threshold>=0 || throw(ArgumentError("threshold must be nonnegative"))
    isfinite(psd_floor) && psd_floor>0 || throw(ArgumentError("psd_floor must be positive"))
    1<=pilot_index<=length(C.bandwidths) || throw(ArgumentError("invalid pilot_index"))
    p,_,q,H = size(C.matrices)
    ell = (3log(max(p,2))+confidence+log(H))/C.n
    domain = get(C.settings,:domain,(0.,1.))
    rates = (domain[2]-domain[1]) ./ (C.min_visits .* C.bandwidths)
    u = sqrt(ell).*(1 .+rates) .+ell.*(1 .+rates.^2)
    v = sqrt(ell).*rates .+ell.*(1 .+rates.^2)
    P,W = zeros(p,p,q),zeros(p,p,q)
    diagnostics = NamedTuple[]
    for t in 1:q
        displacement, violation = 0.0,0.0
        if method===:scaled_lasso
            S = C.matrices[:,:,t,pilot_index]
            E = eigen(Symmetric(S))
            Splus = E.vectors*Diagonal(max.(E.values,psd_floor))*transpose(E.vectors)
            displacement = max(0,psd_floor-minimum(E.values))
            W[:,:,t] .= Matrix(CovarianceSelection.slassoinv(Splus,C.n;λ₀=lambda0,symmetrize=false,ntasks=1))
            theta = diag(W[:,:,t])
        else
            selected = method===:aclime ? [pilot_index] : collect(1:H)
            result = bias_aclime(C.matrices[:,:,t,selected];u=u[selected],v=v[selected],
                bias=bias_constant.*C.bandwidths[selected].^2,
                pilot_index=findfirst(==(pilot_index),selected),K1,K2,noise_variance,
                bias_norm=method===:aclime ? :linf : :l2,diagonal_floor=C.n^-2)
            W[:,:,t] .= result.columns
            theta = result.pilot
            violation = maximum(result.violations)
        end
        S = interval_symmetrize(W[:,:,t])
        cap = min(1.0,C.n/max(opnorm(S,1),eps()))
        P[:,:,t] .= cap.*S
        push!(diagnostics,(;status=:optimal,max_constraint_violation=violation,
            covariance_psd_displacement=displacement,cap,pilot_minimum=minimum(theta)))
    end
    adjacency = abs.(P).>threshold
    for t in 1:q,j in 1:p
        adjacency[j,j,t] = false
    end
    positive = BitVector([isposdef(Symmetric(P[:,:,t])) for t in 1:q])
    settings = (;method,pilot_index,bias_constant,K1,K2,confidence,noise_variance,
        lambda0,psd_floor,threshold,symmetrization=:interval_l1,solver_tolerance=1e-7,ell)
    DirectPrecisionFit(copy(C.grid),P,W,adjacency,positive,diagnostics,C,settings)
end

function fit_direct(data::FDA.MultivariateFunctionalData,grid;
        bandwidths,design=:independent,mean=(j,t)->0.0,boundary=:error,kwargs...)
    C = direct_covariance(data,grid;bandwidths,design,mean,boundary)
    fit_direct(C;kwargs...)
end
