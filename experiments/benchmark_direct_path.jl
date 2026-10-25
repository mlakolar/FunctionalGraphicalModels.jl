# Validate and time reused models on one observed time slice of the saved
# full-size production canaries. No covariance reconstruction is repeated.
include("dfgm2020.jl")

function benchmark_direct_path(case,out)
    mkpath(out)
    cache=deserialize(joinpath(case,"covariances.bin")); C=cache.direct
    p,_,q,H=size(C.matrices)
    time_index=argmin(abs.(C.grid.-.5))
    ell=(3log(max(p,2))+log(H))/C.n
    domain=get(C.settings,:domain,(0.,1.))
    rates=(domain[2]-domain[1])./(C.min_visits.*C.bandwidths)
    u=sqrt(ell).*(1 .+rates).+ell.*(1 .+rates.^2)
    v=sqrt(ell).*rates.+ell.*(1 .+rates.^2)
    ks=exp.(range(log(.03),log(2.);length=7))
    open(joinpath(out,"benchmark.csv"),"w") do io
        row(io,"p","method","repeat","candidates","cold_seconds","path_seconds",
            "max_budget_difference","max_column_difference","max_violation","errors")
        for method in (:aclime,:l2_aclime),repeat in 1:2
            selected=method===:aclime ? [2] : collect(1:H)
            A=C.matrices[:,:,time_index,selected]
            normkind=method===:aclime ? :linf : :l2
            pilot=findfirst(==(2),selected)
            tick=time()
            cold=[bias_aclime(A;u=u[selected],v=v[selected],bias=C.bandwidths[selected].^2,
                pilot_index=pilot,K2,noise_variance=.5,bias_norm=normkind,
                diagonal_floor=C.n^-2) for K2 in ks]
            cold_seconds=time()-tick
            GC.gc()
            tick=time()
            path=bias_aclime_path(A,ks;u=u[selected],v=v[selected],bias=C.bandwidths[selected].^2,
                pilot_index=pilot,noise_variance=.5,bias_norm=normkind,diagonal_floor=C.n^-2)
            path_seconds=time()-tick
            failed=count(!isempty,path.errors)
            budget=maximum(maximum(abs.(cold[k].budgets-path.budgets[k])) for k in eachindex(ks))
            columns=maximum(maximum(abs.(cold[k].columns-path.columns[k])) for k in eachindex(ks))
            violation=maximum(maximum(v) for v in path.violations)
            row(io,p,method,repeat,length(ks),cold_seconds,path_seconds,budget,columns,violation,failed)
            @assert failed==0
            @assert budget<=1e-4*(1+maximum(maximum(c.budgets) for c in cold))
            @assert violation<=1e-5*(1+maximum(maximum(v) for v in path.budgets))
            println("p=$p method=$method repeat=$repeat cold=$cold_seconds path=$path_seconds budgetgap=$budget columngap=$columns")
            flush(stdout)
            GC.gc()
        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run under Slurm")
    benchmark_direct_path(option("case",""),option("out","results/path-benchmark"))
end
