# A seeded integration experiment, not a published simulation design.
using FunctionalGraphicalModels, LinearAlgebra, Random, TOML

function main(args=ARGS)
    length(args) == 1 || error("usage: julia --project experiments/other_methods_smoke.jl NEW_OUTPUT_DIRECTORY")
    out = abspath(only(args))
    ispath(out) && error("output already exists: $out")
    mkpath(out)
    seed = 20261024
    rng = MersenneTwister(seed)
    n,p,m,q = 80,4,2,13
    grid = collect(range(0,1; length=q))
    Q1 = Matrix{Float64}(I,p,p); Q1[1,2]=Q1[2,1]=-0.35
    Q2 = Matrix{Float64}(I,p,p); Q2[3,4]=Q2[4,3]=-0.4
    covariance1,covariance2 = inv(Q1),0.3*inv(Q2)
    a = randn(rng,n,p)*cholesky(Symmetric(covariance1)).U
    b = randn(rng,n,p)*cholesky(Symmetric(covariance2)).U
    phi(t) = sqrt(2)*cos(2pi*t)
    values = zeros(n,p,q)
    truth = falses(p,p)
    truth[1,2]=truth[2,1]=truth[3,4]=truth[4,3]=true
    for t in 1:q
        values[:,:,t] .= a .+ phi(grid[t]).*b
    end
    data = DenseFunctionalData(values,grid)
    ps = fit_separable(data; lambda=0.1,alpha=0.5,ncomponents=m)
    ps.diagnostics.converged || error("partial-separability solver did not converge")
    dense = fit_doubly(data;lambda=0.1,ncomponents=m)
    times = Matrix{Vector{Float64}}(undef,n,p)
    observed = similar(times)
    for j in 1:p, i in 1:n
        times[i,j] = sort!(rand(rng,6+mod(i+j,4)))
        observed[i,j] = a[i,j] .+ b[i,j].*phi.(times[i,j]) .+ 0.1.*randn(rng,length(times[i,j]))
    end
    sparse = fit_doubly(SparseFunctionalData(times,observed);lambda=0.1,
                       grid,bandwidth=0.25,ncomponents=m,noise_variance=0.01)
    all(isfinite,sparse.precisions) || error("nonfinite sparse DFGM output")
    open(joinpath(out,"metrics.csv"),"w") do io
        println(io,"method,time,target,tp,fp,tn,fn,tpr,fpr,f1")
        v=graph_metrics(ps.adjacency,truth)
        println(io,join(("partial_separability","all","latent_functional_union",v.tp,v.fp,v.tn,v.fn,v.tpr,v.fpr,v.f1),','))
        for t in 1:q
            # This latent comparison is descriptive for the bounded sparse
            # schedules. It is not equality of the DFGM reconstruction target.
            exact = inv(covariance1 + phi(grid[t])^2*covariance2)
            target = abs.(exact) .> 1e-9
            for j in 1:p; target[j,j]=false; end
            for (name,fit) in (("dfgm_dense",dense),("dfgm_sparse",sparse))
                v=graph_metrics(fit.adjacency[:,:,t],target)
                label = name == "dfgm_sparse" ? "latent_comparison_only" : "latent_pointwise"
                println(io,join((name,grid[t],label,v.tp,v.fp,v.tn,v.fn,v.tpr,v.fpr,v.f1),','))
            end
        end
    end
    settings = Dict("seed"=>seed,"n"=>n,"p"=>p,"components"=>m,"grid_points"=>q,
        "julia_version"=>string(VERSION),"source_revision"=>get(ENV,"FGM_SOURCE_REV","not_recorded"),
        "status"=>"integration_smoke_not_published_design","lambda"=>0.1,
        "sparse_bandwidth"=>0.25,"sparse_noise_variance"=>0.01,
        "ps_converged"=>ps.diagnostics.converged,"ps_primal"=>ps.diagnostics.primal_residual,
        "ps_dual"=>ps.diagnostics.dual_residual,
        "dense_clime_positive_definite_count"=>count(dense.positive_definite),
        "sparse_clime_positive_definite_count"=>count(sparse.positive_definite),
        "dfgm_sparse_target"=>"conditional-score reconstructed process, not generally latent pointwise precision")
    open(joinpath(out,"settings.toml"),"w") do io; TOML.print(io,settings;sorted=true); end
    println("other methods integration smoke passed; outputs: $out")
end
main()
