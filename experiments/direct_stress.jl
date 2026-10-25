# Additional controlled experiments: exact sparse target at the center, known
# operator bias, stationary control, degree and visit/noise sensitivity.
include("dfgm2020.jl")

function stress_data(rng;p=50,n=200,m=8,degree=10,curvature=2.,noise_variance=.5)
    1<=degree<p || error("invalid degree")
    O = Matrix{Float64}(I,p,p)
    for j in 2:degree+1
        O[1,1] += .3/degree; O[j,j] += .3/degree
        O[1,j] = O[j,1] = -.3/degree
    end
    S = inv(O)
    loading = zeros(p); loading[1] = .4
    rest = S-loading*transpose(loading)
    Z = randn(rng,n,p)*transpose(cholesky(Symmetric(rest)).L)
    U = randn(rng,n)
    nodes = FDA.IrregularFunctionalData[]
    for j in 1:p
        times,values = Vector{Float64}[],Vector{Float64}[]
        for i in 1:n
            ti = sort!(rand(rng,m))
            yi = [Z[i,j]+(1+curvature*(t-.5)^2)*loading[j]*U[i]+
                sqrt(noise_variance)*randn(rng) for t in ti]
            push!(times,ti);push!(values,yi)
        end
        push!(nodes,FDA.IrregularFunctionalData(times,values;domain=(0.,1.)))
    end
    # At .5, covariance is exactly S. For Epanechnikov, second moment=1/5.
    # Exact operator bias is ((1+c*h²/5)^2-1)*norm(loading)^2.
    exact_bias(h) = ((1+curvature*h^2/5)^2-1)*sum(abs2,loading)
    (;data=FDA.MultivariateFunctionalData(nodes),precision=O,exact_bias)
end

function stress_worker(index,out;replicates=3)
    configurations = [(;degree,m,curvature,noise_variance)
        for degree in (2,10) for m in (8,30) for curvature in (0.,4.) for noise_variance in (.5,2.)]
    config = configurations[index]
    mkpath(out)
    isfile(joinpath(out,"stress.csv")) && error("existing stress results; choose a new attempt directory")
    open(joinpath(out,"settings.toml"),"w") do io
        TOML.print(io,Dict("index"=>index,"replicates"=>replicates,
            "degree"=>config.degree,"m"=>config.m,"curvature"=>config.curvature,
            "noise_variance"=>config.noise_variance,"p"=>50,"n"=>200,
            "seed_recipe"=>"2026102500 + 10000*index + replicate",
            "selection"=>"heldout subjects: covariance precision residual; no latent truth"))
    end
    hs = [.1,.2,.3]
    failures = 0
    open(joinpath(out,"stress.csv"),"w") do io
        row(io,"replicate","method","candidate","status","validation_loss","selected","operator","normalized_squared_frobenius","error")
        for rep in 1:replicates
            rng = MersenneTwister(2026102500+10000index+rep)
            sim = stress_data(rng;config...)
            train = sim.data[1:160]
            validation = sim.data[161:200]
            Ct = direct_covariance(train,[.5];bandwidths=hs)
            Cv = direct_covariance(validation,[.5];bandwidths=[.2])
            # Strictly positive envelope also covers stationary zero-bias case.
            bconstant = max(.01,maximum(sim.exact_bias(h)/h^2 for h in hs))
            for method in (:aclime,:l2_aclime,:scaled_lasso)
                candidates = method===:scaled_lasso ? [1.] : [.1,.3,1.]
                fits = []
                for (k,c) in enumerate(candidates)
                    try
                        f = fit_direct(Ct;method,pilot_index=2,bias_constant=bconstant,
                            noise_variance=config.noise_variance,K2=c)
                        loss = norm(Cv.matrices[:,:,1,1]*f.precisions[:,:,1]-I)
                        push!(fits,(k,f,loss))
                    catch e
                        failures += 1
                        row(io,rep,method,k,"failed",NaN,false,NaN,NaN,sprint(showerror,e))
                    end
                end
                isempty(fits) && continue
                best = argmin(last.(fits))
                for (r,(k,f,loss)) in enumerate(fits)
                    D = f.precisions[:,:,1]-sim.precision
                    row(io,rep,method,k,"ok",loss,r==best,opnorm(D),sum(abs2,D)/50,"")
                end
            end
        end
    end
    open(joinpath(out,"complete.toml"),"w") do io
        TOML.print(io,Dict("fit_failures"=>failures,"expected_fits"=>7replicates,
            "status"=>failures==0 ? "completed" : "completed_with_failures"))
    end
    println("stress index=$index complete")
    failures==0
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run scientific work inside Slurm")
    ok = stress_worker(parse(Int,option("setting","1")),option("out","results/direct-stress");
        replicates=parse(Int,option("replicates","3")))
    exit(ok ? 0 : 1)
end
