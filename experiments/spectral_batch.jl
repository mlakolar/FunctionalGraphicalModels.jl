# One persistent Julia process per GNU Parallel/srun worker; source and inputs
# are private to the job. Completed scientific outputs are never fitted again.
include("spectral_comparison.jl")
include("audit_spectral_saved.jl")
tasks=ARGS[1];output=ARGS[2];ntasks=parse(Int,ARGS[3]);failed=false
for (sequence,line) in enumerate(eachline(tasks))
    isempty(line) && continue
    kind,index,rep=split(line,'\t');s=parse(Int,index);r=parse(Int,rep)
    out=joinpath(output,"$kind-$s-$r");ENV["DFGM_CASE_SEQUENCE"]=string(sequence)
    try
        if isfile(joinpath(out,"complete.toml"))
            settings=TOML.parsefile(joinpath(out,"settings.toml"))
            settings["kind"]==kind && settings["setting_index"]==s && settings["replicate"]==r && settings["seed"]==2026102500+10000s+r || error("completed output identity differs")
            settings["source_sha256"]==source_hashes() || error("completed output belongs to a different source; reconcile journal")
            if !isfile(joinpath(out,"independent-check.toml"));audit_spectral_saved(out);end
            println("retained completed spectral dataset $kind/$s/$r");flush(stdout)
            continue
        end
        ok=spectral_dataset(kind,s,r,out;ntasks)
        audit_spectral_saved(out)
        failed |= !ok
    catch e
        failed=true
        println(stderr,"spectral dataset failed $kind/$s/$r: ",sprint(showerror,e));flush(stderr)
        if isdir(out)
            open(joinpath(out,"pipeline-failure.txt"),"w") do io;showerror(io,e,catch_backtrace());end
        end
    end
end
exit(failed ? 1 : 0)
