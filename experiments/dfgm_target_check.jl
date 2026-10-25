# Deterministic target/integration/bias diagnostics. This does not tune fits.
include("dfgm2020.jl")

function smoothed_fourier(t,h;points=4097)
    weight,_=FunctionalGraphicalModels._equivalent_kernel(t,h,(0.,1.),:local_linear)
    lo,hi=max(-1.,-t/h),min(1.,(1-t)/h)
    u=collect(range(lo,hi;length=points))
    w=quadrature_weights(u)
    sum((w[r]*weight(u[r])/(hi-lo).*FunctionalGraphicalModels._dfgm_fourier(t+h*u[r])
        for r in eachindex(u)))
end

function target_check(out)
    mkpath(out)
    open(joinpath(out,"target-check.csv"),"w") do io
        row(io,"p","quantity","value","status")
        for p in (50,100)
            a=simulate_dfgm2020(MersenneTwister(1229);p,n=2,design=:irregular,
                schedule=:T6_9,integration_points=101)
            b=simulate_dfgm2020(MersenneTwister(1229);p,n=2,design=:irregular,
                schedule=:T6_9,integration_points=201)
            discrepancy=maximum(opnorm(a.precision[:,:,t]-a.nominal_precision[:,:,t]) for t in eachindex(a.grid))
            resolution=maximum(opnorm(a.precision[:,:,t]-b.precision[:,:,t]) for t in eachindex(a.grid))
            row(io,p,"max_nominal_actual_operator_gap",discrepancy,"numerical")
            row(io,p,"max_101_201_operator_gap",resolution,"numerical")
            row(io,p,"nominal_graph_resolution_gap",maximum(abs.(a.nominal_precision-b.nominal_precision)),"numerical")
            biasratio=0.; entryratio=0.
            for h in (.15,.25,.4),(t,u) in enumerate(a.grid)
                phi=smoothed_fourier(u,h)
                F=zeros(p,5p)
                for j in 1:p
                    F[j,(j-1)*5+1:j*5].=phi
                end
                J=F*a.coefficient_covariance*transpose(F)-a.covariance[:,:,t]
                biasratio=max(biasratio,opnorm(J)/h^2)
                entryratio=max(entryratio,maximum(abs,J)/h^2)
            end
            row(io,p,"max_operator_bias_divided_h_squared",biasratio,"quadrature_diagnostic_not_certificate")
            row(io,p,"max_entry_bias_divided_h_squared",entryratio,"quadrature_diagnostic_not_certificate")
        end
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    get(ENV,"SLURM_JOB_ID","")!="" || error("run scientific work inside Slurm")
    target_check(option("out","results/dfgm-target-check"))
end
