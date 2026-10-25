# Frozen campaign runtime

These Project/Manifest files belong to common-fast-v2-source, Julia1.13.2.
Manifest SHA256: b3c2e4c3d606a301e65429c19c00a6b41ae664a12f41626c59e8b92e27a4de5e.
Its git-tree hashes pin the dependency contents despite historical repo-rev
labels main/master. Do not run Pkg.update for campaign recreation. Current
package CI intentionally tests newer dependency masters and is a separate
environment; FDA noise defaults and API names have since changed.

A bounded source comparison verified that EVERY src/*.jl file in the frozen
runtime and Project.toml match published commit
d284642d20e8b8d03d9f3b89e29d101be6834c58 byte for byte. source-hashes.json
records those file hashes. All seven entries in spectral-v5/source.json
match published commit8c5110542531af74c65b8ae33242389f36ed740d. The actual
corrected persistent batch file matches commitcaa3a84a754ccd8a9bcf1d2c1d8765b21b908d6c.
The original v5 batch file was not the dispatcher used for production.

To recreate the numerical runtime in a separate checkout, check out d284642,
copy this Manifest.toml into its root, and instantiate that locked environment
with Julia1.13.2. The new spectral program/worker/audit files can be extracted
from8c51105 into a separate source folder, then loaded by the committed
corrected batch dispatcher. Follow the dispatch shell scripts' DFGM_ROOT,
DFGM_RUNTIME_PROJECT and SPECTRAL_SOURCE paths; do not modify a running source.
The recipes use the original seed IDs and recorded cached-data paths.

The v5 worker has the documented aggregate precision-clock bug. This affects
only those timing fields. Use the separate controlled timing scripts for
wall-time evidence; the current public worker also repairs those clocks.
Preserved v4 canaries and original baseline variants remain separately
identified by their source manifests and per-dataset settings. A contemporary
working checkout is not a substitute for one of these frozen variants.

This receipt checks source/environment identity only. It runs no statistical
fit and does not repeat passed scientific tests.
