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

All240 new paper and1600 new control datasets are complete and independently
checked, with infeasible and numerical candidate outcomes preserved. The
[dated comparison report](../../docs/SPECTRAL_COMPARISON_20261026.md) includes
paired MCSE, norm/graph/PD/selection/ROC coverage, failures, controlled timing
and inspected figures. It still lacks17 baseline outputs in setting6; the
original recovery jobs remain active. This is not final completion. See
[protocol](../../docs/SPECTRAL_L2_ACLIME.md) for prespecified constants and scopes.
Source manifests identify every frozen variant; original baseline sources and
fits are preserved. Paper calibration constants are empirical sensitivity inputs.

The read-only selection diagnosis13448886 separates3120/5920 valid pointwise
selections from7/240 complete curves and records the zero-column mechanism.
Full-curve timing13448656 completed on two saved cases,93.00s and2255.36s
with four workers, preserving all numerical/infeasible outcomes. CI7b46b44
passed all four Julia platforms; currentLinux passed4602/4602 checks. These
receipts supplement the dated report without refitting baseline simulations.
