# CI compatibility repair

Run40232051532 failed at the base branch naming check, before dependency
resolution or Julia tests. Its merge base was40c80c9. The existing nested
covariance moment functions M0/M1/M2 violate that new rule. They are renamed
to moment0/moment1/moment2; reversing just those identifiers reproduces the
entire original direct.jl byte for byte. The current base naming checker
passes locally. No expression, kernel weight, solver or calibration changes.
All frozen simulations retain their original source files/hashes.

The batch utility now uses a function scope and explicitly selects the frozen
numerical source. Persistent execution was checked in13446911, with both
cases independently verified. Its nonzero allocation status is from preserved
small-score infeasible candidates, not a pipeline exception.

The additional finite-sample feasibility job13446925 passed26 checks on actual
independent uniform schedules, unit noise, known centering and a constant
Gaussian covariance kernel. It is an implementation check, not a risk theorem.

Run40233802138 then passed4570 checks (including all72 spectral tests),
but three covariance tests errored after KernelSmoothing master renamed
localpolynomial/smoothermatrix to local_polynomial/smoother_matrix. The
experiment adapter accepts both APIs with identical arguments. Focused
old-lock golden checks13447971 passed23/23. The first focused bundle13447958
omitted its raw-reference test helper and failed; that receipt is retained.

The frozen v4/v5 experiment driver also captured a shared fit_tick variable
in parallel time workers. Its aggregate precision wall/path times are invalid.
The current driver uses separate point and phase clocks. Numerical outputs,
individual column timers, covariance and scoring timers were unaffected.
No primary fits were repeated to repair timing. A separate controlled
benchmark13447505 uses the same two saved covariance inputs for all methods.

Report13447385 mixed baseline methods sharing candidate indices in its paper
selection filter. Its paper comparisons are superseded by corrected
report13447865, which filters the method explicitly and checks selected
metric coverage. The original report and its qualification remain preserved.
Saved-matrix audit13448207 independently recomputed996 metric/score/selection
rows across five cases: maximum norm error8.17e-14 and score error7.28e-12.
