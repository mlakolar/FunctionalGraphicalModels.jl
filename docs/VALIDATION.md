# Initial validation receipt — 24 October 2026

All numerical execution was on USC Discovery through Slurm, using Julia 1.13.2,
2 CPUs, 8 GB, two Julia threads and one BLAS thread per allocation.

| Job | Source at execution | Outcome |
|---|---|---|
| 13350744 | 24b3342 | Foundation suite: 72/72 checks; 25 seconds |
| 13350766 | 388a1f8 | 108 passed; one unqualified `checksquare` error, corrected in the estimator branch |
| 13350788 | 589ad03 | Cancelled after 4m34s stuck in filesystem I/O during setup; tests did not start |
| 13350865 | bace77c | All four model families: 316/316 checks; 2m16s |
| 13350955 | 49c6b8d | Complete suite: 4055/4055 checks, including graph-generator/metric cases; both experiment drivers passed; 1m42s; peak RSS about 1.2 GB |

The CI path repair rebased the stack after source staging. The tested `src/`
tree for job 13350955 and the resulting reproduction branch are identical:
`39f24da31577e3201d2dcbe57e72cc7ddd3ac161`. The rebase changed only the
GitHub workflow's bootstrap path. Logs preserve the original source revisions.

## Numerical scope

- Qiao smoke: n=40, p=6, T=41, true coefficient dimension 5, retained M=3,
  fixed spline L=7, seed 20261024. Three generating models, six penalties and
  two estimators: **36 fits converged; zero failures**. These are integration
  experiments, not the paper's dimensions or data-driven tuning.
- Other-methods smoke: n=80, p=4, two latent components, 13 evaluation points.
  Partially separable joint graphical lasso converged (primal residual
  9.44e-9; dual residual 6.60e-6). Dense and irregular DFGM each produced 13
  finite pointwise CLIME estimates; all happened to be positive definite.
  This does not imply that CLIME always has that property.
- The sparse DFGM latent-graph comparison is descriptive. With bounded sparse
  measurements the target is the conditional-score reconstructed process.
- Small-sample graph performance is not evidence that one method dominates
  the others. No published AUC table or real-data result is claimed here.

Selected numerical outputs and settings are committed under
[`validation/2026-10-24`](../validation/2026-10-24). Full small outputs, saved
supports, the test log and dependency Manifest are in the local ignored
`results/discovery-13350955/` directory. Remote receipts remain at
`/home1/mkolar/juliaPkg/runs/20261024-fgm-development/`.
The parent workspace's `.slurm/runs/20261024-fgm-development.json` records all
attempts, commands and resources. Failed and cancelled attempts are retained.
