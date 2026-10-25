# Approved reference driver at published commit 578fc1d. Its covariance,
# reconstruction and generating functions are unchanged by path reuse.
const DFGM_REFERENCE_COVARIANCE_DRIVER =
    "d6e83e097b6c6bdacdd53dc5930b3e3315ef6f7108cdb4677f179c4a5174e99c"
const DFGM_APPROVED_COVARIANCE_DRIVERS = (
    DFGM_REFERENCE_COVARIANCE_DRIVER,
    "b167369d4a39c2a34a671416ed34f9b0068464ec8d5d4ff05ced60126a31c4a8",
    "fb3b0bc7930e4924474e9bb8d38a794aa4d3b8493bc47d6c422f6f349e777147")

function load_dfgm_covariance_cache(folder;setting_index,replicate,seed,n,p,
        design,schedule,mode,integration_points,work_points)
    prior=TOML.parsefile(joinpath(folder,"settings.toml"))
    expected=Dict("setting_index"=>setting_index,"replicate"=>replicate,"seed"=>seed,
        "n"=>n,"p"=>p,"design"=>string(design),"schedule"=>string(schedule),
        "mode"=>string(mode),"integration_points"=>integration_points,"work_points"=>work_points,
        "noise_variance"=>.5)
    for (key,value) in expected
        get(prior,key,nothing)==value || error("cache setting mismatch: $key")
    end
    root=normpath(joinpath(@__DIR__,".."))
    for file in ("src/direct.jl","src/dfgm_simulation.jl",
            "experiments/dfgm_common.jl","experiments/dfgm_factorized.jl")
        prior["source_sha256"][file]==bytes2hex(sha256(read(joinpath(root,file)))) ||
            error("cache covariance source mismatch: $file")
    end
    driver=bytes2hex(sha256(read(joinpath(@__DIR__,"dfgm2020.jl"))))
    prior["source_sha256"]["experiments/dfgm2020.jl"] in
        (driver,DFGM_APPROVED_COVARIANCE_DRIVERS...) || error("unverified covariance driver in cache")
    prior["dependency_manifest_sha256"]==bytes2hex(sha256(read(joinpath(root,"Manifest.toml")))) ||
        error("cache dependency lock differs")
    prior["julia_version"]==string(VERSION) || error("cache Julia version differs")
    covariance=deserialize(joinpath(folder,"covariances.bin"))
    targets=deserialize(joinpath(folder,"targets.bin"))
    C=covariance.direct
    C.n==n && size(C.matrices,1)==p || error("cache dimensions differ")
    C.grid==covariance.grid==targets.grid==prior["evaluation_grid"] || error("cache grids differ")
    C.bandwidths==[.15,.25,.4] || error("cache bandwidth menu differs")
    covariance.reconstructed===nothing && error("cache reconstruction failed; regenerate rather than reuse")
    receipt=(;folder=abspath(folder),covariance_sha256=bytes2hex(sha256(read(joinpath(folder,"covariances.bin")))),
        targets_sha256=bytes2hex(sha256(read(joinpath(folder,"targets.bin")))),
        original_driver_sha256=prior["source_sha256"]["experiments/dfgm2020.jl"])
    (;covariance,targets,receipt)
end
