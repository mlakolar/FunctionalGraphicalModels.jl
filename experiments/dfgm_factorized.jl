# Exact marginal product-kernel sufficient statistics for irregular visits.
# For each subject, sum_{r!=s} a_r b_s = sum(a)*sum(b)-sum(a.*b).
# Applied to all six normal moments and three response moments, this reduces
# the work from pair enumeration to kernel matrices and BLAS products. It is
# the same local-linear objective and fixed-mean whole-curve CV as FDA.
import KernelSmoothing
isdefined(@__MODULE__, :_dfgm_local_polynomial) || include(joinpath(@__DIR__,"dfgm_kernel_api.jl"))

function _dfgm_residual_pool(fd,meanfit)
    times,residual,subjects = Float64[],Float64[],Int[]
    weights = Float64[]
    for i in 1:FDA.nobs(fd)
        ti,yi=fd[i]
        length(ti)>=2 || error("covariance requires two visits per curve")
        append!(times,ti);append!(residual,yi-FDA.evaluate(meanfit,ti))
        append!(subjects,fill(i,length(ti)))
        push!(weights,1/(length(ti)*(length(ti)-1)))
    end
    (;times,residual,subjects,weights)
end

function _dfgm_factor_surface(pool,queries,h;exclude=0)
    q,n = length(queries),length(pool.weights)
    kernel=KernelSmoothing.GaussianKernel()
    U = (transpose(pool.times).-queries)./h
    K0 = kernel.(U); K1=K0.*U; K2=K1.*U
    S0,S1,S2,Y0,Y1 = [zeros(q,n) for _ in 1:5]
    for a in eachindex(pool.times)
        i=pool.subjects[a];v=pool.residual[a]
        @views S0[:,i] .+= K0[:,a]
        @views S1[:,i] .+= K1[:,a]
        @views S2[:,i] .+= K2[:,a]
        @views Y0[:,i] .+= v.*K0[:,a]
        @views Y1[:,i] .+= v.*K1[:,a]
    end
    sw=copy(pool.weights)
    exclude==0 || (sw[exclude]=0.)
    ow=sw[pool.subjects]
    yw=ow.*pool.residual.^2
    product(A,B,w)=(A.*transpose(w))*transpose(B)
    A00=product(S0,S0,sw)-product(K0,K0,ow)
    A10=product(S1,S0,sw)-product(K1,K0,ow)
    A20=product(S2,S0,sw)-product(K2,K0,ow)
    A11=product(S1,S1,sw)-product(K1,K1,ow)
    B00=product(Y0,Y0,sw)-product(K0,K0,yw)
    B10=product(Y1,Y0,sw)-product(K1,K0,yw)
    surface=zeros(q,q)
    for s in 1:q,t in 1:q
        A=[A00[s,t] A10[s,t] A10[t,s];
           A10[s,t] A20[s,t] A11[s,t];
           A10[t,s] A11[s,t] A20[t,s]]
        b=[B00[s,t],B10[s,t],B10[t,s]]
        C=cholesky(Symmetric(A);check=false)
        isposdef(C) && minimum(abs2,diag(C.U))>100eps(Float64)*maximum(abs,A) ||
            error("factorized covariance has a singular local fit at ($(queries[s]),$(queries[t]))")
        surface[s,t]=(C\b)[1]
    end
    all(isfinite,surface) || error("nonfinite factorized covariance")
    surface
end

function _dfgm_factorized_covariance(fd,meanfit,grid;bandwidth)
    pool=_dfgm_residual_pool(fd,meanfit)
    n=FDA.nobs(fd)
    selection=nothing
    h=if bandwidth isa FDA.CurveCV
        bandwidth.folds==n && bandwidth.pilot===nothing || error("explicit leave-one-curve CV required")
        errors=Float64[]
        for candidate in bandwidth.bandwidths
            score=0.
            for i in 1:n
                ti,yi=fd[i]
                ri=yi-FDA.evaluate(meanfit,ti)
                predicted=_dfgm_factor_surface(pool,ti,candidate;exclude=i)
                score+=sum((abs2(ri[r]*ri[s]-predicted[r,s])
                    for r in eachindex(ti) for s in eachindex(ti) if r!=s);init=0.)*
                    pool.weights[i]
            end
            push!(errors,score/n)
        end
        selection=(;method=:exact_factorized_curve_cv,bandwidths=copy(bandwidth.bandwidths),
            errors,assignments=collect(1:n))
        bandwidth.bandwidths[argmin(errors)]
    else
        Float64(bandwidth)
    end
    surface=_dfgm_factor_surface(pool,grid,h)
    surface=(surface+transpose(surface))/2
    m=FDA.nobs_per_curve(fd)
    diagonal_weights=[1/m[i] for i in pool.subjects]
    kernel=KernelSmoothing.GaussianKernel()
    diagonal=_dfgm_local_polynomial(pool.times,pool.residual.^2,grid;
        bandwidth=h,kernel,degree=1,weights=diagonal_weights)
    lower,upper=first(grid)+(last(grid)-first(grid))/4,last(grid)-(last(grid)-first(grid))/4
    integration_grid=vcat(lower,grid[lower .< grid .< upper],upper)
    differences=FDA.GridFunctionalData(reshape(diagonal-diag(surface),1,:),grid)
    integrand=vec(FDA.evaluate(differences,integration_grid,FDA.PiecewiseLinear()))
    raw=2dot(quadrature_weights(integration_grid),integrand)/(last(grid)-first(grid))
    details=(;bandwidth=h,degree=1,minobs=0,folds=0,selection=nothing,
        diagonal_source=:covariance,raw_variance=raw,truncated=raw<0,
        bounds=(lower,upper),integration_grid,diagonal,aggregation=:exact_product_factorization)
    FDA.SmoothedCovariance(collect(grid),surface,quadrature_weights(grid),h,h,kernel,
        :subj,nothing,1,FDA.sampling_design(fd),selection===nothing ? 0 : n,
        selection,max(0.,raw),details)
end
