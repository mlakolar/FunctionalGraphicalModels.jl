using FunctionalGraphicalModels, LinearAlgebra, Random, TOML, Test
const F=FunctionalGraphicalModels
const FDA=F.FDA
Base.include(F,joinpath(ARGS[1],"src/spectral_direct.jl"))
rng=MersenneTwister(2026102601);N=120;m=8;p=2
S=[1. .2;.2 1.5];truth=inv(S);Z=randn(rng,N,p)*transpose(cholesky(Symmetric(S)).L)
nodes=FDA.IrregularFunctionalData[]
for j in 1:p
 times=[sort!(rand(rng,m)) for _ in 1:N]
 vals=[Z[i,j].+randn(rng,m) for i in 1:N]
 push!(nodes,FDA.IrregularFunctionalData(times,vals;domain=(0.,1.)))
end
cal=F.spectral_calibration(N,p,m;L=1.,kappa_bound=2.,h_max=.25)
C=F.spectral_covariance(FDA.MultivariateFunctionalData(nodes),[.5],cal).covariance
d=max.(diag(C.matrices[:,:,1,cal.diagonal_index]),N^-2);B=4 ./d
max_score_ratio=0.
@testset "Actual independent-schedule unit-noise truth feasibility" begin
 for j in 1:p
  @test norm(truth[:,j])<=B[j]
  for h in eachindex(cal.bandwidths)
   lambda=2 .*sqrt.(d.*B[j]).*(cal.u[h]+B[j]*cal.v[h])
   residual=abs.(C.matrices[:,:,1,h]*truth[:,j]-Float64.(1:p .==j))
   @test all(residual.<=lambda)
  end
 end
end
open(ARGS[2],"w") do io
 TOML.print(io,Dict("status"=>"passed","N"=>N,"p"=>p,"m"=>m,"noise_variance"=>1.,"centering"=>"known_zero",
 "design"=>"independent_uniform","target_time"=>.5,"curvature"=>0.,"supplied_L"=>1.,"condition_bound"=>2.,"seed"=>2026102601,
 "scope"=>"Finite-sample true columns satisfy all Program(P) constraints with b=0 and t=column l1; implementation check, not a risk theorem"))
end
