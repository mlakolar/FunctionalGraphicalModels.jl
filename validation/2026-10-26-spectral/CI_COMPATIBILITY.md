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
