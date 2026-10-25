# Independent replay of reporting quantities from saved arrays. No refits,
# simulation, driver inclusion or calls to its metrics/roc/BIC functions.
using FunctionalGraphicalModels,LinearAlgebra,Serialization,TOML
BLAS.set_num_threads(1)
function audit_saved(folder,out)
    settings=TOML.parsefile(joinpath(folder,"settings.toml"))
    paths=deserialize(joinpath(folder,"paths.bin"))
    covariances=deserialize(joinpath(folder,"covariances.bin"))
    targets=deserialize(joinpath(folder,"targets.bin"))
    n,p=settings["n"],settings["p"];cutoff=settings["threshold"]
    mkpath(out)
    open(joinpath(out,"independent-metrics.csv"),"w") do io
        println(io,"method,target,time_index,candidate,bic,operator,frobenius,normalized_squared_frobenius,fpr,tpr,edges,nonedges")
        for method in sort!(collect(keys(paths));by=string)
            S=method in (:DFGM_C,:DFGM_G) ? covariances.reconstructed :
              method in (:DGM_C,:DGM_G) ? covariances.dynamic :
              method in (:GGL,:Naive_C,:Naive_G) ? covariances.raw : covariances.direct.matrices[:,:,:,2]
            for (candidate,precision) in paths[method],t in axes(precision,3)
                P=precision[:,:,t]
                @assert size(P)==(p,p) && all(isfinite,P)
                @assert maximum(abs.(P-transpose(P)))<=1e-10*(1+maximum(abs,P))
                factor=cholesky(Symmetric(P);check=false)
                bic=if issuccess(factor)
                    dimension=sum(abs(P[i,j])>cutoff for i in 1:p,j in 1:p)
                    n*(tr(P*S[:,:,t])-2sum(log,diag(factor.L)))+log(n)*dimension
                else
                    Inf
                end
                for (target,truth) in ((:actual,targets.actual),(:nominal,targets.nominal))
                    D=P-truth[:,:,t]
                    operator=maximum(abs,eigvals(Symmetric(D)))
                    squared=sum(abs2,D)
                    positives=falsepositives=truepositives=negatives=0
                    for i in 1:p-1,j in i+1:p
                        truth_edge=abs(truth[i,j,t])>cutoff
                        fitted_edge=abs(P[i,j])>cutoff
                        positives+=truth_edge;negatives+=!truth_edge
                        truepositives+=truth_edge&&fitted_edge
                        falsepositives+=!truth_edge&&fitted_edge
                    end
                    fpr=negatives==0 ? NaN : falsepositives/negatives
                    tpr=positives==0 ? NaN : truepositives/positives
                    println(io,join((method,target,t,candidate,bic,operator,sqrt(squared),squared/p,
                        fpr,tpr,positives,negatives),","))
                end
            end
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run under Slurm")
    length(ARGS)==2 || error("expected saved dataset and fresh audit output folder")
    ispath(ARGS[2]) && error("audit output exists; preserve it and choose a new attempt")
    audit_saved(ARGS...)
end
