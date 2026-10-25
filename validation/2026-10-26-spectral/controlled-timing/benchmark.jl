# Controlled timing only. No regenerated observations or primary baseline outputs.
using FunctionalGraphicalModels, LinearAlgebra, Serialization, TOML, SHA, Sockets
const F=FunctionalGraphicalModels
const R="/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct"
const S=joinpath(R,"spectral-20261026-v5-source")
Base.include(F,joinpath(S,"src/spectral_direct.jl"))
quote_csv(x)="\""*replace(string(x),'"'=>"\"\"")*"\""
row(io,x...)=(println(io,join(quote_csv.(x),','));flush(io))
function one_time(C,index)
 F.DirectCovariance(C.grid[index:index],C.bandwidths,C.matrices[:,:,index:index,:],C.n,C.min_visits,C.pair_subjects[:,index:index,:],C.settings)
end
# The clock belongs to this call. No shared closure variable can reset it.
function timed_call(f)
 before=time_ns()
 try
  result=f();return (;seconds=Float64(time_ns()-before)/1e9,result,error="")
 catch e
  return (;seconds=Float64(time_ns()-before)/1e9,result=nothing,error=sprint(showerror,e))
 end
end
function spectral_path(C,ks)
 calibration_start=time_ns();cal=F.spectral_calibration(C.n,size(C.matrices,1),C.min_visits;L=1.,kappa_bound=10.,h_max=.25)
 diagonal=max.(diag(C.matrices[:,:,1,cal.diagonal_index]),C.n^-2)
 calibration_seconds=Float64(time_ns()-calibration_start)/1e9
 path_start=time_ns();r=F.spectral_l2_aclime_path(C.matrices[:,:,1,:],ks;d=diagonal,u=cal.u,v=cal.v,bias=cal.bias,kappa_bound=10.,noise_variance=.5)
 column_path_seconds=Float64(time_ns()-path_start)/1e9
 symmetry_start=time_ns();matrices=[];errors=String[]
 for x in r
  if all(isfinite,x.columns)
   try
    P=F.interval_symmetrize(x.columns);norm1=opnorm(P,1);P.*=norm1==0 ? 1. : min(1.,C.n/norm1);push!(matrices,P)
   catch e;push!(errors,sprint(showerror,e));end
  end
 end
 (;candidates=length(ks),accepted=count(x->x.status===:optimal,r),infeasible=count(x->x.status===:infeasible_zero_fallback,r),
 numerical=count(x->x.status===:numerical_failure,r)+length(errors),zero_matrices=count(P->all(iszero,P),matrices),
 calibration_seconds,column_path_seconds,symmetry_cap_seconds=Float64(time_ns()-symmetry_start)/1e9,errors=join(errors,"; "))
end
function legacy_path(C,method,ks)
 r=F.fit_direct_path(C,ks;method,pilot_index=2,bias_constant=1.,noise_variance=.5,ntasks=1)
 (;candidates=length(ks),accepted=count(isempty,r.errors),infeasible=0,numerical=count(!isempty,r.errors),zero_matrices=count(f->f!==nothing&&all(iszero,f.precisions),r.fits),
 calibration_seconds=NaN,column_path_seconds=NaN,symmetry_cap_seconds=NaN,errors=join(filter(!isempty,r.errors),"; "))
end
function scaled(C,ls)
 accepted=0;errors=String[];zero=0
 for lambda0 in ls
  try
   r=F.fit_direct(C;method=:scaled_lasso,pilot_index=2,bias_constant=1.,noise_variance=.5,lambda0)
   accepted+=1;zero+=all(iszero,r.precisions)
  catch e;push!(errors,sprint(showerror,e));end
 end
 (;candidates=length(ls),accepted,infeasible=0,numerical=length(errors),zero_matrices=zero,calibration_seconds=NaN,column_path_seconds=NaN,symmetry_cap_seconds=NaN,errors=join(errors,"; "))
end
out=ARGS[1];mkpath(out)
open(joinpath(out,"settings.toml"),"w") do io
 TOML.print(io,Dict("source"=>S,"runtime"=>joinpath(R,"common-fast-v2-source"),"benchmark_sha256"=>bytes2hex(sha256(read(@__FILE__))),
 "cpu_model"=>Sys.cpu_info()[1].model,"hostname"=>gethostname(),"threads"=>Threads.nthreads(),"solver_threads"=>1,
 "scope"=>"Separate mechanical one-time-point benchmark on existing matrices; not additional independent simulations or changes to primary estimates",
 "menus"=>"six original K2 grid positions1,7,13,19,25,31 at fixed bias1; original sixSZ penalties; exact-grid six spectral scores. Different constraints/calibration retained; equal candidate counts do not imply equal statistical tuning.",
 "failure_classification"=>"spectral numerical failures and certified infeasibility are separate; legacy failures remain unclassified with raw errors",
 "warmup"=>"sample0 first call in this process, sample1 warm repeat; no discarded calls; same allocation sequential methods",
 "reason"=>"Independent controlled wall clocks after detecting shared aggregate timing variable in original parallel spectral worker"))
end
open(joinpath(out,"timings.csv"),"w") do io
 row(io,"setting","replicate","time_index","target_time","method","sample","seconds","candidates","accepted","infeasible","failed_or_numerical","zero_matrices","calibration_seconds","column_path_seconds","symmetry_cap_seconds","error")
 for setting in (7,6)
  original=joinpath(R,"spectral-input-view-v1/outputs/paper-v4/setting-$setting-rep-1")
  C0=deserialize(joinpath(original,"covariances.bin")).direct
  folder=setting==7 ? joinpath(R,"outputs/spectral-v5/paper-7-1") : joinpath(R,"outputs/spectral-v4-canary/paper-6-1")
  C1=deserialize(joinpath(folder,"covariances.bin")).direct
  index=cld(length(C0.grid),2);A=one_time(C0,index);B=one_time(C1,index)
  A.grid==B.grid || error("input grids differ")
  A.matrices[:,:,:,2]==B.matrices[:,:,:,findfirst(==(.25),B.bandwidths)] || error("shared input covariance differs")
  ks=collect(exp.(range(log(.03),log(2.);length=31)))[[1,7,13,19,25,31]]
  ls=Any[:quantile,.2,.35,.5,.75,sqrt(2log(size(A.matrices,1))/A.n)*(1+1/(A.min_visits*.25))]
  for sample in 0:1
   for method in (:scaled_lasso,:scaled_lasso_tuned,:aclime,:l2_aclime,:spectral_l2_aclime,:spectral_l2_aclime_tuned)
    call=method===:scaled_lasso ? ()->scaled(A,ls[1:1]) : method===:scaled_lasso_tuned ? ()->scaled(A,ls) :
         method in (:aclime,:l2_aclime) ? ()->legacy_path(A,method,ks) : method===:spectral_l2_aclime ? ()->spectral_path(B,[2.]) : ()->spectral_path(B,[.01,.03,.1,.3,1.,2.])
    result=timed_call(call);x=result.result
    if x===nothing
     row(io,setting,1,index,A.grid[1],method,sample,result.seconds,method in (:scaled_lasso,:spectral_l2_aclime) ? 1 : 6,0,0,NaN,NaN,NaN,NaN,NaN,result.error)
    else
     row(io,setting,1,index,A.grid[1],method,sample,result.seconds,x.candidates,x.accepted,x.infeasible,x.numerical,x.zero_matrices,x.calibration_seconds,x.column_path_seconds,x.symmetry_cap_seconds,x.errors)
    end
    println("timing setting=$setting method=$method sample=$sample seconds=$(result.seconds)");flush(stdout);GC.gc(false)
   end
  end
 end
end
open(joinpath(out,"complete.toml"),"w") do io;TOML.print(io,Dict("status"=>"completed","scope"=>"Timing benchmark only; errors/outcomes retained in every row"));end
