using LinearAlgebra, Serialization, TOML, SHA
function audit_spectral_saved(folder)
    settings=TOML.parsefile(joinpath(folder,"settings.toml"));cal=settings
    covariance=deserialize(joinpath(folder,"covariances.bin")).direct
    witnesses=deserialize(joinpath(folder,"audit-witnesses.bin"))
    stages=deserialize(joinpath(folder,"diagnostics.bin"))
    N=settings["n"];p=settings["p"];hs=settings["bandwidths"];h_d=settings["diagonal_bandwidth"]
    d=max.(diag(covariance.matrices[:,:,1,findfirst(==(h_d),hs)]),N^-2)
    B=2settings["kappa_bound"]./d
    bias=settings["C_K"]*settings["L"].*hs.^2
    ell=settings["ell"];q=settings["domain_length"]./(covariance.min_visits.*hs)
    u=sqrt(ell).*(1 .+q).+ell.*(1 .+q.^2);v=sqrt(ell).*q.+ell.*(1 .+q.^2)
    maximum_relative_violation=0.;checked_columns=0;skipped_failure_columns=0
    for (index,r) in enumerate(witnesses)
        r.d==d && r.B==B || error("diagonal calibration differs from saved covariance")
        for j in 1:p
            beta=r.columns[:,j];t=r.budgets[j];b=r.bias_vectors[:,:,j]
            record=stages[1].stages[index].column_records[j]
            if record.status in (:infeasible_zero_fallback,:numerical_failure)
                skipped_failure_columns+=1;continue
            end
            all(isfinite,beta) || error("nonfinite accepted witness column")
            rel=max(0.,norm(beta,1)-t)/(1+t)
            rel=max(rel,max(0.,norm(beta)-B[j])/(1+B[j]))
            for h in eachindex(hs)
                lambda=r.Kscore.*sqrt.(d.*B[j]).*(u[h]+settings["noise_variance"]*B[j]*v[h])
                residual=covariance.matrices[:,:,1,h]*beta-Float64.(1:p .== j)-b[:,h]
                rel=max(rel,maximum(max.(0.,abs.(residual).-lambda)./(1 .+lambda)),
                    max(0.,norm(b[:,h])-bias[h]*t)/(1+bias[h]*t),max(0.,norm(b[:,h])-bias[h]*B[j])/(1+bias[h]*B[j]))
            end
            maximum_relative_violation=max(maximum_relative_violation,rel);checked_columns+=1
        end
    end
    maximum_relative_violation<=1e-5 || error("independent primal replay failed: $maximum_relative_violation")
    # Check every time/candidate label and cap/diagonal from its recorded stage.
    output=deserialize(joinpath(folder,"paths.bin"))
    for t in eachindex(stages),k in eachindex(output.paths)
        s=stages[t].stages[k];records=s.column_records
        for x in records
            x.status===:infeasible_zero_fallback && !(x.termination=="INFEASIBLE" && x.dual_status=="INFEASIBILITY_CERTIFICATE") && error("invalid fallback certificate label")
        end
        P=output.paths[k][:,:,t];W=output.columns[k][:,:,t]
        if s.status===:optimal
            all(isfinite,P) && issymmetric(P) || error("invalid accepted precision")
            maximum(abs.(diag(P).-s.cap.*diag(W)))<=1e-7*(1+maximum(abs,P)) || error("symmetrization diagonal changed")
            opnorm(P,1)<=N+1e-5*(1+N) || error("radial cap failed")
        end
    end
    receipt=Dict("status"=>"passed","scope"=>"Independent first-time primal witness replay; all-time column labels, diagonal and cap checks; no fits or risk theorem validation",
        "checked_columns"=>checked_columns,"skipped_failure_columns"=>skipped_failure_columns,"max_relative_violation"=>maximum_relative_violation)
    open(joinpath(folder,"independent-check.toml"),"w") do io;TOML.print(io,receipt);end
    receipt
end
if abspath(PROGRAM_FILE)==@__FILE__;println(audit_spectral_saved(ARGS[1]));end
