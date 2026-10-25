"""Read the explicit reporting scope; preserve legacy behavior when absent."""
import json
import pathlib

OWN_METHODS = frozenset(("aclime", "l2_aclime"))


def load_scope(root):
    path = pathlib.Path(root) / "comparison-scope.json"
    if not path.exists():
        return None
    scope = json.loads(path.read_text())
    assert scope["own_method_replicates"] == 20
    assert scope["other_method_replicates"] == 100
    assert scope["selection"] == "first_fixed_replicate_ids"
    return scope


def included(scope, method, replicate):
    return scope is None or method not in OWN_METHODS or int(replicate) <= scope["own_method_replicates"]
