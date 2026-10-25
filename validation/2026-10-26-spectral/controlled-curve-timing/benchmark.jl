# Purposeful timing repair only. Reuse saved covariances; never replace primary fits.
using FunctionalGraphicalModels, LinearAlgebra, Serialization, TOML, Sockets, SHA
const F=FunctionalGraphicalModels
const R="/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct"
const SOURCE=joinpath(R,"spectral-20261026-v5-source")
Base.include(F,joinpath(SOURCE,"src/spectral_direct.jl"))
function slice_work(C,t,cal)
 calibration_start=time_ns();d=max.(diag(C.matrices[:,:,t,cal.diagonal_index]),C.n^-2)
 diagonal_seconds=Float64(time_ns()-calibration_start)/1e9
 core_start=time_ns();results=F.spectral_l2_aclime_path(C.matrices[:,:,t,:],[.01,.03,.1,.3,1.,2.];
  d,u=cal.u,v=cal.v,bias=cal.bias,kappa_bound=10.,noise_variance=.5)
 core_seconds=Float64(time_ns()-core_start)/1e9
 symmetry_seconds=0.;cap_seconds=0.;accepted=0;infeasible=0;numerical=0
 for x in results
  if all(isfinite,x.columns)
   try
    symmetry_start=time_ns();P=F.interval_symmetrize(x.columns);symmetry_seconds+=Float64(time_ns()-symmetry_start)/1e9
    cap_start=time_ns();norm1=opnorm(P,1);P.*=norm1==0 ? 1. : min(1.,C.n/norm1);cap_seconds+=Float64(time_ns()-cap_start)/1e9
   catch
    numerical+=1;continue
   end
  end
  accepted+=x.status===:optimal;infeasible+=x.status===:infeasible_zero_fallback;numerical+=x.status===:numerical_failure
 end
 (;diagonal_seconds,core_seconds,symmetry_seconds,cap_seconds,accepted,infeasible,numerical)
end
function curve_work(C,cal)
 q=length(C.grid);result=Vector{Any}(undef,q)
 @sync for worker in 1:min(4,q)
  Threads.@spawn for t in worker:4:q
   result[t]=slice_work(C,t,cal)
   println("timed slice=$t/$q core_seconds=$(result[t].core_seconds)");flush(stdout)
  end
 end
 result
end
function timed_curve(C,cal)
 # This start is local to the timing wrapper and cannot be reset by slice tasks.
 phase_start=time_ns();r=curve_work(C,cal);wall=Float64(time_ns()-phase_start)/1e9
 sum(x.core_seconds for x in r)<=4wall+1. || error("private clock/concurrency inconsistency")
 maximum(x.core_seconds for x in r)<=wall+1. || error("slice exceeds whole phase")
 (;wall,results=r)
end
out=ARGS[1];mkpath(out);cases=Dict[]
for setting in (7,6)
 p=joinpath(R,"outputs/spectral-v5/paper-$setting-1");C=deserialize(joinpath(p,"covariances.bin")).direct
 before=time_ns();cal=F.spectral_calibration(C.n,size(C.matrices,1),C.min_visits;L=1.,kappa_bound=10.,h_max=.25)
 calibration_seconds=Float64(time_ns()-before)/1e9
 C.bandwidths==cal.bandwidths || error("grid differs")
 result=timed_curve(C,cal);r=result.results
 receipt=Dict("setting"=>setting,"replicate"=>1,"n"=>C.n,"p"=>size(C.matrices,1),"minimum_visits"=>C.min_visits,
 "times"=>length(C.grid),"widths"=>length(C.bandwidths),"candidates_per_time"=>6,"threads"=>Threads.nthreads(),"solver_threads"=>1,
 "precision_phase_wall_seconds"=>result.wall,"column_model_solve_sum_seconds"=>sum(x.core_seconds for x in r),
 "symmetrization_sum_seconds"=>sum(x.symmetry_seconds for x in r),"cap_sum_seconds"=>sum(x.cap_seconds for x in r),
 "calibration_seconds"=>calibration_seconds,"diagonal_calibration_sum_seconds"=>sum(x.diagonal_seconds for x in r),
 "accepted_candidate_times"=>sum(x.accepted for x in r),"infeasible_candidate_times"=>sum(x.infeasible for x in r),"numerical_candidate_times"=>sum(x.numerical for x in r))
 push!(cases,receipt)
 open(joinpath(out,"case-$setting.toml"),"w") do io;TOML.print(io,receipt);end
end
open(joinpath(out,"receipt.toml"),"w") do io
 TOML.print(io,Dict("status"=>"completed","cpu_model"=>Sys.cpu_info()[1].model,"hostname"=>gethostname(),
 "benchmark_sha256"=>bytes2hex(sha256(read(@__FILE__))),"source"=>SOURCE,"cases"=>cases,
 "scope"=>"Separate two-case full-curve precision timing repair, all exact-grid six-score programs, four private time tasks. Cached covariances; no resimulation or primary-estimate replacement. Excludes covariance construction, checkpoint I/O, likelihood scoring and statistical evaluation."))
end
