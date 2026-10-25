"""Render an explicit task table; this command never submits or computes fits."""
import argparse
import pathlib

p = argparse.ArgumentParser()
p.add_argument("output", type=pathlib.Path)
p.add_argument("--mode", choices=("canary", "pilot", "paper"), default="pilot")
p.add_argument("--replicates", type=int)
p.add_argument("--nlambda", type=int)
args = p.parse_args()
reps = args.replicates if args.replicates is not None else (100 if args.mode == "paper" else 1)
nlambda = args.nlambda if args.nlambda is not None else (31 if args.mode == "paper" else 3)
if reps < 1 or nlambda < 2:
    p.error("replicates must be positive and nlambda at least two")
settings = (1, 7) if args.mode == "canary" else range(1, 13)
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open("x") as f:
    for setting in settings:
        for rep in range(1, reps + 1):
            # Canary exercises dimensions only; pilot uses published n,p.
            n, dimension = (30, 10) if args.mode == "canary" else (0, 0)
            integration = 31 if args.mode == "canary" else 101
            work = 15 if args.mode == "canary" else 31
            mode = "paper" if args.mode == "paper" else "pilot"
            f.write(f"{setting}\t{rep}\t{mode}\t{n}\t{dimension}\t{integration}\t{nlambda}\t{work}\n")
print(f"Wrote {len(settings) * reps} explicit tasks to {args.output}; no jobs submitted")
