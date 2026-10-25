# Development conventions

Types and modules use UpperCamelCase. Descriptive functions, fields and keyword
arguments use snake_case; mutating functions end in `!`. Retain standard numerical
abbreviations, mathematical identifiers and names required by upstream interfaces.
Run `python3 scripts/check_naming.py` to check new declarations.

Run `julia scripts/instantiate.jl` before tests. Workspace dependencies use sibling
master checkouts, including transitive dependencies. Missing checkouts are cloned
from master. Existing checkouts are never reset or silently updated. Package
version numbering and registry changes are deferred.

Run the complete tests with `julia scripts/test_master.jl`. This uses an isolated
master-based environment on Julia 1.10 and newer without changing package deps.
