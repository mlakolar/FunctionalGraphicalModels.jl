# Exact sufficient-statistic aggregation for a single shared observation grid.
# It preserves the pooled local-linear objective and FDA's fixed-mean whole-
# curve CV criterion. This helper is an experiment optimization, not a new FDA
# covariance estimator. Golden tests compare it with FDA on the raw products.
import KernelSmoothing

function _dfgm_common_covariance(fd,meanfit,grid;bandwidth,cache=Dict())
    n = FDA.nobs(fd)
    times = fd[1][1]
    all(fd[i][1]==times for i in 1:n) || error("common aggregation requires one shared grid")
    T = length(times)
    pairs = [(r,s) for r in 1:T for s in 1:T if r!=s]
    x,y = first.(pairs),last.(pairs)
    coordinates = (times[x],times[y])
    mu = FDA.evaluate(meanfit,times)
    R = reduce(vcat,[transpose(fd[i][2]-mu) for i in 1:n])
    products = [R[i,r]*R[i,s] for (r,s) in pairs,i in 1:n]
    average = vec(sum(products;dims=2)./n)
    kernel = KernelSmoothing.GaussianKernel()
    selection = nothing
    h = if bandwidth isa FDA.CurveCV
        bandwidth.folds==n || error("common CV optimization requires leave-one-curve-out folds")
        bandwidth.pilot===nothing || error("explicit common covariance bandwidth candidates required")
        errors = Float64[]
        for candidate in bandwidth.bandwidths
            key = (Tuple(times),candidate)
            H = get!(cache,key) do
                KernelSmoothing.smoothermatrix(coordinates,collect(zip(coordinates...));
                    bandwidth=candidate,kernel,degree=1)
            end
            predictions = (n.*(H*average) .- H*products)./(n-1)
            push!(errors,sum(abs2,products-predictions)/length(products))
        end
        selection = (;method=:exact_common_curve_cv,bandwidths=copy(bandwidth.bandwidths),
            errors,assignments=collect(1:n))
        bandwidth.bandwidths[argmin(errors)]
    else
        Float64(bandwidth)
    end
    surface = KernelSmoothing.localpolynomial(coordinates,average,(grid,grid);
        bandwidth=h,kernel,degree=1)
    surface = (surface+transpose(surface))/2
    squares = vec(sum(abs2,R;dims=1)./n)
    diagonal = KernelSmoothing.localpolynomial(times,squares,grid;bandwidth=h,kernel,degree=1)
    lower,upper = first(grid)+(last(grid)-first(grid))/4,last(grid)-(last(grid)-first(grid))/4
    integration_grid = vcat(lower,grid[lower .< grid .< upper],upper)
    differences = FDA.GridFunctionalData(reshape(diagonal-diag(surface),1,:),grid)
    integrand = vec(FDA.evaluate(differences,integration_grid,FDA.PiecewiseLinear()))
    raw = 2dot(quadrature_weights(integration_grid),integrand)/(last(grid)-first(grid))
    details = (;bandwidth=h,degree=1,minobs=0,folds=0,selection=nothing,
        diagonal_source=:covariance,raw_variance=raw,truncated=raw<0,
        bounds=(lower,upper),integration_grid,diagonal,aggregation=:exact_shared_grid)
    FDA.SmoothedCovariance(collect(grid),surface,quadrature_weights(grid),h,h,kernel,
        :subj,nothing,1,FDA.sampling_design(fd),selection===nothing ? 0 : n,
        selection,max(0.,raw),details)
end
