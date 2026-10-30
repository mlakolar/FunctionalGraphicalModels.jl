# Solver-only deterministic canary; no data generation, targets, loss tuning or RNG.
# julia --project experiments/benchmark_multiscale_path.jl [p=12] [repeats=3]
using FunctionalGraphicalModels, LinearAlgebra

function benchmark_multiscale_path(p=12,repeats=3)
    S=[0.3^abs(j-k) for j in 1:p,k in 1:p]
    specs=[(;matrices=cat([S+(.005h+.01a)*I for h in 1:3]...;dims=3),
        u=[.02,.03,.04],v=[.01,.015,.02],bias=[.01,.02,.03]) for a in 1:3]
    penalties=[.1,.2,.3,.5,.75,1.]
    native(reuse)=multiscale_l2_aclime_path(specs,penalties;noise_variance=.5,reuse)
    old()=[bias_aclime_path(s.matrices,penalties;u=s.u,v=s.v,bias=s.bias,bias_norm=:l2,noise_variance=.5) for s in specs]
    # Compile every branch and check equivalence before collecting timings.
    path=native(true);cold=native(false);jump=old()
    for a in eachindex(specs),b in eachindex(penalties)
        @assert path.results[a][b].status===cold.results[a][b].status===:optimal
        @assert isempty(jump[a].errors[b])
        @assert isapprox(path.results[a][b].budgets,cold.results[a][b].budgets;atol=1e-4,rtol=1e-4)
        @assert isapprox(path.results[a][b].budgets,jump[a].budgets[b];atol=1e-4,rtol=1e-4)
    end
    println("p,repeat,implementation,seconds,bytes,workspace_builds,matrix_updates,rhs_updates")
    methods=[:native_reuse,:native_cold,:jump_path]
    for repeat in 1:repeats
        for method in (iseven(repeat) ? reverse(methods) : methods)
            GC.gc()
            timing=@timed (method===:jump_path ? old() : native(method===:native_reuse))
            counts=method===:jump_path ? (-1,-1,-1) :
                (timing.value.diagnostics.workspace_builds,timing.value.diagnostics.matrix_updates,timing.value.diagnostics.rhs_updates)
            println(join((p,repeat,method,timing.time,timing.bytes,counts...),','))
        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    benchmark_multiscale_path(isempty(ARGS) ? 12 : parse(Int,ARGS[1]),length(ARGS)<2 ? 3 : parse(Int,ARGS[2]))
end
