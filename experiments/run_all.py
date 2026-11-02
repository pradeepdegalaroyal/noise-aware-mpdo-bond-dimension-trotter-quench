"""One-command reproduction: simulation stages -> results/*.csv|json -> figures/*.png.

    python experiments/run_all.py            # full study (about 10 minutes on a laptop)
    python experiments/run_all.py --quick    # smoke run (< 1 minute), written to results_quick/
    python experiments/run_all.py --stage sweep|size|ablation|validation|figures
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "6")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("MKL_NUM_THREADS", "6")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

import pandas as pd  # noqa: E402

import stages  # noqa: E402

FULL = dict(
    sweep=dict(ns=(6, 8), seeds=(0, 1, 2), p2_grid=(0.01, 0.02, 0.035, 0.05, 0.07, 0.1), steps=12, max_cap=64),
    size=dict(ns=(10, 12), seeds=(0, 1, 2), p2=0.05, steps=12, max_cap=32),
    ablation=dict(n=8, seeds=(0, 1, 2), p2=0.05, steps=12, max_cap=64),
)
QUICK = dict(
    sweep=dict(ns=(6,), seeds=(0, 1), p2_grid=(0.03, 0.1), steps=6, max_cap=16),
    size=dict(ns=(8,), seeds=(0,), p2=0.05, steps=6, max_cap=16),
    ablation=dict(n=6, seeds=(0, 1), p2=0.05, steps=6, max_cap=16),
)


def _save(df_rows: list[dict], path: Path, append: bool = False) -> None:
    df = pd.DataFrame(df_rows)
    if append and path.exists():
        df = pd.concat([pd.read_csv(path), df], ignore_index=True)
    df.to_csv(path, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--stage", default="all", choices=["all", "sweep", "size", "ablation", "validation", "figures"])
    args = ap.parse_args()
    cfg = QUICK if args.quick else FULL
    out = ROOT / ("results_quick" if args.quick else "results")
    fig = ROOT / ("figures_quick" if args.quick else "figures")
    out.mkdir(exist_ok=True)
    fig.mkdir(exist_ok=True)
    t0 = time.time()
    want = lambda s: args.stage in ("all", s)  # noqa: E731

    if want("sweep"):
        print("[sweep] noise sweep, exact-referenced", flush=True)
        steps_rows, cap_rows, gap_rows = stages.stage_noise_sweep(**cfg["sweep"])
        _save(steps_rows, out / "chi_star_by_step.csv")
        _save(cap_rows, out / "cap_errors.csv")
        _save(gap_rows, out / "ideal_gap.csv")
    if want("size"):
        print("[size] size scaling", flush=True)
        steps_rows, cap_rows, gap_rows = stages.stage_size_scaling(**cfg["size"])
        _save(steps_rows, out / "chi_star_by_step.csv", append=True)
        _save(cap_rows, out / "cap_errors.csv", append=True)
        _save(gap_rows, out / "ideal_gap.csv", append=True)
    if want("ablation"):
        print("[ablation] noise-source ablation", flush=True)
        _save(stages.stage_ablation(**cfg["ablation"]), out / "ablation.csv")
    if want("validation"):
        print("[validation] Qiskit Aer cross-check", flush=True)
        try:
            (out / "validation.json").write_text(json.dumps(stages.stage_validation(), indent=2))
        except ImportError as exc:  # qiskit optional
            print("  skipped:", exc)
    if want("figures"):
        print("[figures]", flush=True)
        import plots

        summary = plots.make_all(out, fig, cfg["sweep"]["steps"])
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
