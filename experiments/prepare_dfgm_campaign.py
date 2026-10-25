"""Render balanced per-pool task tables. Never submit or run scientific work.

Completed production canaries may later be omitted only after source/settings
and output validation, recorded in the campaign registry. This initial manifest
contains every intended dataset so missing work cannot disappear from coverage.
"""
import argparse
import json
import pathlib


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--pools", type=int, default=64)
    parser.add_argument("--nlambda", type=int, default=31)
    args = parser.parse_args()
    if args.replicates < 1 or args.pools < 1 or args.nlambda < 2:
        parser.error("positive replicate/pool counts and at least two penalties required")
    args.output.mkdir(parents=True)  # Refuse an existing manifest; never overwrite.
    rows = [f"{setting}\t{replicate}\tpaper\t0\t0\t101\t{args.nlambda}\t31\n"
            for setting in range(1, 13) for replicate in range(1, args.replicates + 1)]
    (args.output / "tasks.tsv").write_text("".join(rows))
    pools = []
    for index in range(args.pools):
        name = f"pool-{index + 1:03d}.tsv"
        subset = rows[index::args.pools]
        (args.output / name).write_text("".join(subset))
        pools.append(dict(pool=index + 1, task_file=name, datasets=len(subset)))
    manifest = dict(datasets=len(rows), settings=12, replicates=args.replicates,
                    pools=pools, nlambda=args.nlambda, bias_constants=[.01, .1, 1., 10., 20.],
                    submitted=False, source="Resolve immutable source from campaign registry",
                    note="All intended datasets retained; validated reuse must be recorded separately")
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Rendered {len(rows)} datasets in {args.pools} pools; no jobs submitted")


if __name__ == "__main__":
    main()
