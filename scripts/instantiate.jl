# Explicit URL dependencies work on Julia 1.10 as well as newer Pkg versions.
# Run on the user's compute host; Manifest.toml records exact resolved commits.
using Pkg, TOML
root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(root)
deps = TOML.parsefile(joinpath(root, "Project.toml"))["deps"]
urls = [
    ("ProximalBase", "https://github.com/mlakolar/ProximalBase.jl", "master"),
    ("ProximalSolvers", "https://github.com/mlakolar/ProximalSolvers.jl", "master"),
    ("CovarianceSelection", "https://github.com/mlakolar/CovarianceSelection.jl", "master"),
    ("KernelSmoothing", "https://github.com/mlakolar/KernelSmoothing.jl", "master"),
    ("FunctionalDataAnalysis", "https://github.com/mlakolar/FunctionalDataAnalysis.jl", "main"),
]
specs = Pkg.PackageSpec[]
for (name, url, rev) in urls
    if haskey(deps, name) || (name == "ProximalSolvers" && haskey(deps, "CovarianceSelection"))
        push!(specs, Pkg.PackageSpec(; name, url, rev))
    end
end
# ProximalSolvers must be resolved alongside its unregistered downstream
# package, but must not become a stale direct dependency of this package.
project_path = joinpath(root, "Project.toml")
original_project = read(project_path, String)
try
    isempty(specs) || Pkg.add(specs)
    Pkg.resolve()
    Pkg.instantiate()
finally
    # Newer Pkg versions can also rewrite sources during resolve/instantiate.
    write(project_path, original_project)
end
