# Spectral SOCP implementation checks and pilot

Program(P) is a distinct modification of l2-ACLIME. These are implementation
checks and cost/identity evidence, not a full paired comparison or a theorem.

- Job13446815 passed72 focused checks, including both bias constraints, the
  column ball, diagonal calibration, raw indefinite inputs, certified-infeasible
  fallback and unchanged existing-method outputs.
- Job13446494 checked the existing paper setting7/replicate1 observations and
  actual target, then evaluated one time with13 exact-grid widths. Additional
  covariance cost0.853seconds; six-score path21.102seconds including first JIT.
- Full control cases1:1 and16:1 retain genuine infeasible candidates at small
  score constants. Independent residual checks passed. Failed candidate states
  remain visible; these outputs are not universally converged paths.
- The first test job13446350 failed two assertions: the intended indefinite
  feasibility example used incompatible matrices at two widths. The corrected
  compatible indefinite example passed in13446367. Original logs are retained.

Full paired paper/control comparisons are still running. See
[protocol](../../docs/SPECTRAL_L2_ACLIME.md) for prespecified constants and scopes.
Source manifests identify every frozen variant; original baseline sources and
fits are preserved. Paper calibration constants are empirical sensitivity inputs.
