# Exact shared-grid product-kernel algebra. The old dfgm_common.jl helper is
# retained as an independent reference. No pair involving the same visit is
# included, and the covariance/mean/CV/noise conventions are unchanged.
import KernelSmoothing
isdefined(@__MODULE__, :_dfgm_local_polynomial) || include(joinpath(@__DIR__,"dfgm_kernel_api.jl"))


function _dfgm_common_fast_mean(fd,grid;bandwidth,cache=Dict())
    n=FDA.nobs(fd);times=fd[1][1]
    all(fd[i][1]==times for i in 1:n) || error("common mean requires one shared grid")
    kernel=KernelSmoothing.GaussianKernel()
    selection=nothing
    h=if bandwidth isa FDA.CurveCV
        bandwidth.folds==n && bandwidth.pilot===nothing || error("explicit leave-one-curve CV required")
        n>=2 || error("curve CV requires at least two subjects")
        Y=reduce(hcat,[fd[i][2] for i in 1:n])
        average=vec(sum(Y;dims=2)./n)
        errors=Float64[]
        for candidate in bandwidth.bandwidths
            key=(:exact_common_mean_weights,Tuple(times),candidate)
            H=get!(cache,key) do
                _dfgm_smoother_matrix(times,times;bandwidth=candidate,kernel,degree=1)
            end
            predicted=(n.*(H*average).-H*Y)./(n-1)
            push!(errors,sum(abs2,Y-predicted)/length(Y))
        end
        all(isfinite,errors) || error("nonfinite common mean CV criterion")
        selection=(;method=:exact_common_mean_curve_cv,bandwidths=copy(bandwidth.bandwidths),
            errors,assignments=collect(1:n))
        bandwidth.bandwidths[argmin(errors)]
    else
        Float64(bandwidth)
    end
    # Preserve FDA's public final fit and all version-specific metadata. Only
    # the selection calculation changes; this also supports its newer fields.
    fixed=FDA.fit(FDA.SmoothedMean,fd;grid,bandwidth=h,kernel,weighting=:subj)
    selection===nothing && return fixed
    FDA.SmoothedMean((name===:selection ? selection : name===:folds ? n :
        getproperty(fixed,name) for name in fieldnames(typeof(fixed)))...)
end

function _dfgm_common_geometry(times,queries,h)
    isfinite(h) && h>0 || error("bandwidth must be positive")
    q,T=length(queries),length(times)
    T>=3 || error("common local-linear covariance needs at least three visits")
    U=(transpose(times).-queries)./h
    K0=KernelSmoothing.GaussianKernel().(U)
    K1=K0.*U;K2=K1.*U
    s0=vec(sum(K0;dims=2));s1=vec(sum(K1;dims=2));s2=vec(sum(K2;dims=2))
    a,b,c=[zeros(q,q) for _ in 1:3]
    correction=zeros(q*q,T)
    for t in 1:q,s in 1:q
        # Product sums minus repeated-visit terms give the six normal moments.
        A00=s0[s]*s0[t]-dot(K0[s,:],K0[t,:])
        A10=s1[s]*s0[t]-dot(K1[s,:],K0[t,:])
        A01=s0[s]*s1[t]-dot(K0[s,:],K1[t,:])
        A20=s2[s]*s0[t]-dot(K2[s,:],K0[t,:])
        A02=s0[s]*s2[t]-dot(K0[s,:],K2[t,:])
        A11=s1[s]*s1[t]-dot(K1[s,:],K1[t,:])
        A=[A00 A10 A01;A10 A20 A11;A01 A11 A02]
        F=cholesky(Symmetric(A);check=false)
        isposdef(F) && minimum(abs2,diag(F.U))>100eps(Float64)*maximum(abs,A) ||
            error("common product geometry has a singular local fit at ($(queries[s]),$(queries[t]))")
        coefficients=F\[1.,0.,0.]
        a[s,t],b[s,t],c[s,t]=coefficients
        index=s+(t-1)*q
        for r in 1:T
            correction[index,r]=K0[s,r]*K0[t,r]*
                (a[s,t]+b[s,t]*U[s,r]+c[s,t]*U[t,r])
        end
    end
    (;K0,K1,a,b,c,correction)
end

function _dfgm_common_product_fits(geometry,R)
    q=size(geometry.K0,1);n=size(R,1)
    Y0=geometry.K0*transpose(R);Y1=geometry.K1*transpose(R)
    fits=-(geometry.correction*transpose(R.^2))
    for i in 1:n,t in 1:q,s in 1:q
        fits[s+(t-1)*q,i]+=geometry.a[s,t]*Y0[s,i]*Y0[t,i]+
            geometry.b[s,t]*Y1[s,i]*Y0[t,i]+geometry.c[s,t]*Y0[s,i]*Y1[t,i]
    end
    fits
end

function _dfgm_common_product_surface(geometry,S)
    q=size(geometry.K0,1)
    geometry.a.*(geometry.K0*S*transpose(geometry.K0))+
        geometry.b.*(geometry.K1*S*transpose(geometry.K0))+
        geometry.c.*(geometry.K0*S*transpose(geometry.K1))-
        reshape(geometry.correction*diag(S),q,q)
end

function _dfgm_common_fast_covariance(fd,meanfit,grid;bandwidth,cache=Dict())
    n=FDA.nobs(fd);n>=2 || error("curve CV requires at least two subjects")
    times=fd[1][1];T=length(times)
    all(fd[i][1]==times for i in 1:n) || error("common aggregation requires one shared grid")
    mu=FDA.evaluate(meanfit,times)
    R=reduce(vcat,[transpose(fd[i][2]-mu) for i in 1:n])
    pairs=[(r,s) for r in 1:T for s in 1:T if r!=s]
    indices=[r+(s-1)*T for (r,s) in pairs]
    products=[R[i,r]*R[i,s] for (r,s) in pairs,i in 1:n]
    function geometry(points,h)
        key=(:exact_common_product_geometry,Tuple(times),Tuple(points),h)
        get!(() -> _dfgm_common_geometry(times,points,h),cache,key)
    end
    selection=nothing
    h=if bandwidth isa FDA.CurveCV
        bandwidth.folds==n && bandwidth.pilot===nothing || error("explicit leave-one-curve CV required")
        errors=Float64[]
        for candidate in bandwidth.bandwidths
            fitted=_dfgm_common_product_fits(geometry(times,candidate),R)[indices,:]
            average=vec(sum(fitted;dims=2)./n)
            predicted=(n.*average.-fitted)./(n-1)
            push!(errors,sum(abs2,products-predicted)/length(products))
        end
        selection=(;method=:exact_common_product_curve_cv,bandwidths=copy(bandwidth.bandwidths),
            errors,assignments=collect(1:n))
        bandwidth.bandwidths[argmin(errors)]
    else
        Float64(bandwidth)
    end
    surface=_dfgm_common_product_surface(geometry(grid,h),transpose(R)*R/n)
    surface=(surface+transpose(surface))/2
    squares=vec(sum(abs2,R;dims=1)./n)
    kernel=KernelSmoothing.GaussianKernel()
    diagonal=_dfgm_local_polynomial(times,squares,grid;bandwidth=h,kernel,degree=1)
    lower,upper=first(grid)+(last(grid)-first(grid))/4,last(grid)-(last(grid)-first(grid))/4
    integration_grid=vcat(lower,grid[lower .< grid .< upper],upper)
    differences=FDA.GridFunctionalData(reshape(diagonal-diag(surface),1,:),grid)
    integrand=vec(FDA.evaluate(differences,integration_grid,FDA.PiecewiseLinear()))
    raw=2dot(quadrature_weights(integration_grid),integrand)/(last(grid)-first(grid))
    all(isfinite,surface) && isfinite(raw) || error("nonfinite common product covariance")
    details=(;bandwidth=h,degree=1,minobs=0,folds=0,selection=nothing,
        diagonal_source=:covariance,raw_variance=raw,truncated=raw<0,
        bounds=(lower,upper),integration_grid,diagonal,aggregation=:exact_common_product_geometry)
    FDA.SmoothedCovariance(collect(grid),surface,quadrature_weights(grid),h,h,kernel,
        :subj,nothing,1,FDA.sampling_design(fd),selection===nothing ? 0 : n,
        selection,max(0.,raw),details)
end
