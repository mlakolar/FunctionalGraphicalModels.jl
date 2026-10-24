"""The twelve Section 4 settings of Qiao et al. (2020), Tables 2 and 3.
The published experiment uses n=200 and 100 replicates per setting.
"""
dfgm2020_settings() = [(;design,p,schedule,n=200,replicates=100)
    for design in (:common,:irregular) for p in (50,100)
    for schedule in (design===:common ? (:T10,:T25,:T50) : (:T6_9,:T20_30,:T50))]

_dfgm_fourier(t) = [1.0,sqrt(2)*sinpi(2t),sqrt(2)*cospi(2t),
    sqrt(2)*sinpi(4t),sqrt(2)*cospi(4t)]

function _dfgm_dynamic(rng,p)
    p>=10 || throw(ArgumentError("dynamic simulation requires p>=10"))
    edges = [(j,k) for j in 1:p for k in j+1:p]
    # Published counts are 40/10 for p=50, 160/40 for p=100.
    initial_count = max(1,round(Int,.016p^2))
    change_count = max(1,round(Int,.004p^2))
    initial_count+change_count<=length(edges) || throw(ArgumentError("graph too small"))
    initial = zeros(length(edges))
    for e in Random.randperm(rng,length(edges))[1:initial_count]
        initial[e] = .1+.2Random.rand(rng)
    end
    events = NamedTuple[]
    decaying_until = fill(-Inf,length(edges))
    function edgeweight(e,t)
        idx = findlast(ev->ev.edge==e && ev.start<=t,events)
        idx===nothing && return initial[e]
        ev = events[idx]
        ev.value+(ev.target-ev.value)*clamp((t-ev.start)/.5,0,1)
    end
    weights(t) = [edgeweight(e,t) for e in eachindex(edges)]
    for start in 0.0:.1:.9
        w = weights(start)
        growing = findall(e->w[e]<=1e-10,eachindex(edges))
        decaying = findall(e->w[e]>1e-10 && decaying_until[e]<=start+1e-10,eachindex(edges))
        ng,nd = change_count,change_count
        min(length(growing),length(decaying))>=change_count || error("not enough eligible dynamic edges")
        for e in growing[Random.randperm(rng,length(growing))[1:ng]]
            push!(events,(;start,edge=e,value=0.,target=.1+.2Random.rand(rng)))
        end
        for e in decaying[Random.randperm(rng,length(decaying))[1:nd]]
            push!(events,(;start,edge=e,value=w[e],target=0.))
            decaying_until[e] = start+.5
        end
    end
    function precision(t)
        S = Matrix(.25I,p,p)
        for (e,(j,k)) in enumerate(edges)
            w = edgeweight(e,t)
            w = max(0,w)
            S[j,j] += w; S[k,k] += w
            S[j,k] = S[k,j] = -w
        end
        S
    end
    precision,(;edges,initial,events,initial_count,change_count,
        overlap_policy=:replace_growing_ramp_when_selected_to_decay)
end

"""
    simulate_dfgm2020(rng; p=50,n=200,design=:irregular,schedule=:T6_9,
                      integration_points=101,noise_variance=.5,
                      evaluation_grid=nothing)

Biometrika Section 4/Appendix A design: evolving Laplacian precision with .25
baseline, transitions of length .5, Gaussian coregionalization correlation
exp(-(u-v)^2/2), and projection onto five orthonormal Fourier functions.
Numerical integration uses the trapezoidal grid of `integration_points` points.
All random draws use the supplied RNG. Common grids are stratified into five
equal intervals, with T/5 uniform points per interval. Irregular schedules use
independent uniform times and the specified discrete-uniform visit counts.

`nominal_precision` is the dynamic precision before Fourier projection, as in
the paper's construction. `precision` is the actual pointwise precision of the
generated finite-rank process. Both are returned: they need not agree. The
paper does not specify the integration resolution, Fourier ordering, or how
to resolve overlapping edge transitions; these conventions are recorded.
Use integration-grid sensitivity checks before claiming a table replication.
"""
function simulate_dfgm2020(rng::Random.AbstractRNG; p::Integer=50,n::Integer=200,
        design::Symbol=:irregular,schedule::Symbol=:T6_9,
        integration_points::Integer=101,noise_variance::Real=.5,evaluation_grid=nothing)
    n>=2 && integration_points>=21 || throw(ArgumentError("invalid sample/integration size"))
    isfinite(noise_variance) && noise_variance>=0 || throw(ArgumentError("invalid noise variance"))
    design in (:common,:irregular) || throw(ArgumentError("unknown design"))
    schedule in (:T10,:T25,:T50,:T6_9,:T20_30) || throw(ArgumentError("unknown visit schedule"))
    design===:common && !(schedule in (:T10,:T25,:T50)) && throw(ArgumentError("invalid common schedule"))
    nominal,graph = _dfgm_dynamic(rng,p)
    ig = collect(range(0,1;length=integration_points))
    qw = quadrature_weights(ig)
    phi = reduce(vcat,transpose.(_dfgm_fourier.(ig)))
    loadings = [Matrix(cholesky(Symmetric(inv(nominal(t)))).L) for t in ig]
    R = exp.(-((ig .- transpose(ig)).^2)./2)
    E = eigen(Symmetric(R))
    keep = findall(>(1e-12*maximum(E.values)),E.values)
    coefficient_covariance = zeros(5p,5p)
    for r in keep
        B = zeros(5p,p)
        for u in eachindex(ig),j in 1:p,l in 1:5
            B[(j-1)*5+l,:] .+= qw[u]*phi[u,l]*sqrt(E.values[r])*E.vectors[u,r].*loadings[u][j,:]
        end
        coefficient_covariance .+= B*transpose(B)
    end
    ce = eigen(Symmetric(coefficient_covariance))
    minimum(ce.values)>=-1e-8*maximum(ce.values) || error("invalid projected covariance")
    root = ce.vectors*Diagonal(sqrt.(max.(ce.values,0)))
    coefficients = Random.randn(rng,n,5p)*transpose(root)
    common = if design===:common
        T = schedule===:T10 ? 10 : schedule===:T25 ? 25 : 50
        sort!([(b+Random.rand(rng))/5 for b in 0:4 for _ in 1:div(T,5)])
    else
        nothing
    end
    nodes = FDA.IrregularFunctionalData[]
    for j in 1:p
        times,values = Vector{Float64}[],Vector{Float64}[]
        for i in 1:n
            T = schedule===:T6_9 ? Random.rand(rng,6:9) :
                schedule===:T20_30 ? Random.rand(rng,20:30) :
                schedule===:T10 ? 10 : schedule===:T25 ? 25 : 50
            ti = common===nothing ? sort!(Random.rand(rng,T)) : copy(common)
            yi = [dot(_dfgm_fourier(t),coefficients[i,(j-1)*5+1:j*5])+
                sqrt(noise_variance)*Random.randn(rng) for t in ti]
            push!(times,ti); push!(values,yi)
        end
        push!(nodes,FDA.IrregularFunctionalData(times,values;domain=(0.,1.)))
    end
    grid = evaluation_grid===nothing ? (common===nothing ? collect(0:.05:1) : common) : Float64.(evaluation_grid)
    covariance,precision,nominal_precision = zeros(p,p,length(grid)),zeros(p,p,length(grid)),zeros(p,p,length(grid))
    for (v,t) in enumerate(grid)
        F = zeros(p,5p)
        for j in 1:p
            F[j,(j-1)*5+1:j*5] .= _dfgm_fourier(t)
        end
        covariance[:,:,v] .= F*coefficient_covariance*transpose(F)
        precision[:,:,v] .= inv(Symmetric(covariance[:,:,v]))
        nominal_precision[:,:,v] .= nominal(t)
    end
    (;data=FDA.MultivariateFunctionalData(nodes),grid,covariance,precision,
        nominal_precision,coefficient_covariance,graph,
        settings=(;p,n,design,schedule,integration_points,noise_variance,
            fourier=:constant_sin1_cos1_sin2_cos2,rho_sigma=1.,
            integration=:trapezoidal,projection_rank=5))
end
