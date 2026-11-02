"""Experiment stages. Each stage returns plain row dictionaries that run_all.py writes to results/."""

from __future__ import annotations

import time

import numpy as np

from noisy_mpdo import NoiseSpec, QuenchSpec, run_exact, run_ideal_mps, run_mpdo
from noisy_mpdo.search import chi_star_curve, obs_error, search_caps

DELTAS = (0.01, 0.03)  # tolerance on observable error: 0.01 equals the shot noise of 10^4 shots


def _caps_pow2(lo: int, hi: int) -> list[int]:
    out, c = [], lo
    while c <= hi:
        out.append(c)
        c *= 2
    return out


def _ideal_search(spec: QuenchSpec, delta: float):
    full = 2 ** (spec.n // 2)
    ref = run_ideal_mps(spec, full).obs
    return ref, search_caps(lambda c: run_ideal_mps(spec, c), ref, _caps_pow2(2, full), delta)


def _reference(spec: QuenchSpec, noise: NoiseSpec, exact_max_n: int, ref_cap: int) -> np.ndarray:
    if spec.n <= exact_max_n:
        return run_exact(spec, noise).obs
    return run_mpdo(spec, noise, ref_cap).obs


def _rows(kind, spec, p2, noise_tag, caps, errs, secs, reached, step_rows, cap_rows):
    for k, c in enumerate(caps):
        cap_rows.append(dict(kind=kind, n=spec.n, seed=spec.seed, p2=p2, noise=noise_tag, cap=c,
                             final_err=float(errs[k].max()), seconds=float(secs[k]), bond_reached=int(reached[k])))
    for dl in DELTAS:
        star = chi_star_curve(caps, errs, dl)
        for t, v in enumerate(star):
            step_rows.append(dict(kind=kind, n=spec.n, seed=spec.seed, p2=p2, noise=noise_tag, delta=dl,
                                  step=t, chi_star=float(v), max_cap=int(caps[-1])))


def stage_noise_sweep(ns, seeds, p2_grid, steps, max_cap, exact_max_n=10, ref_cap=64):
    """chi*(t) of the noise-aware MPDO and of the ideal MPS, plus the ideal-vs-noisy gap."""
    step_rows, cap_rows, gap_rows = [], [], []
    for n in ns:
        for seed in seeds:
            spec = QuenchSpec(n=n, steps=steps, seed=seed)
            ideal_ref, (caps, errs, secs, reached) = _ideal_search(spec, DELTAS[0])
            _rows("ideal_mps", spec, 0.0, "none", caps, errs, secs, reached, step_rows, cap_rows)
            for p2 in p2_grid:
                noise = NoiseSpec.from_p2(p2)
                ref = _reference(spec, noise, exact_max_n, ref_cap)
                gap_rows.append(dict(n=n, seed=seed, p2=p2, ideal_vs_noisy_gap=float(obs_error(ideal_ref, ref).max()),
                                     signal=float(np.abs(ref[1:]).max())))
                caps, errs, secs, reached = search_caps(lambda c: run_mpdo(spec, noise, c), ref,
                                                        _caps_pow2(2, max_cap), DELTAS[0])
                _rows("mpdo", spec, p2, "full", caps, errs, secs, reached, step_rows, cap_rows)
            print(f"  sweep n={n} seed={seed} done", flush=True)
    return step_rows, cap_rows, gap_rows


def stage_ablation(n, seeds, p2, steps, max_cap, fixed_cap=16):
    """Cumulative noise-source ablation: chi* (coarse) and truncation error at a fixed cap (fine).

    Order: no noise -> + two-qubit depolarising -> + idle dephasing -> + single-qubit depolarising.
    """
    full = NoiseSpec.from_p2(p2)
    configs = [("no noise", NoiseSpec()),
               ("+ two-qubit depolarising", NoiseSpec(p2=full.p2)),
               ("+ idle dephasing", NoiseSpec(p2=full.p2, pd=full.pd)),
               ("+ single-qubit depolarising", full)]
    rows = []
    for seed in seeds:
        spec = QuenchSpec(n=n, steps=steps, seed=seed)
        _, (caps, errs, _, _) = _ideal_search(spec, DELTAS[0])
        psi_sq = float(chi_star_curve(caps, errs, DELTAS[0])[-1]) ** 2
        for name, noise in configs:
            ref = run_exact(spec, noise).obs
            err_fixed = float(obs_error(run_mpdo(spec, noise, fixed_cap).obs, ref).max())
            if noise.is_noiseless():
                chi = psi_sq  # exact bond dimension of |psi><psi| is chi_psi squared
            else:
                caps, errs, _, _ = search_caps(lambda c: run_mpdo(spec, noise, c), ref, _caps_pow2(2, max_cap), DELTAS[0])
                chi = float(chi_star_curve(caps, errs, DELTAS[0])[-1])
            rows.append(dict(seed=seed, config=name, chi_star=chi, fixed_cap=fixed_cap, err_fixed_cap=err_fixed))
        print(f"  ablation seed={seed} done", flush=True)
    return rows


def stage_size_scaling(ns, seeds, p2, steps, max_cap, exact_max_n=10, ref_cap=64):
    """Final-step chi* versus qubit number n at fixed p2 (MPDO and ideal MPS)."""
    return stage_noise_sweep(ns, seeds, [p2], steps, max_cap, exact_max_n, ref_cap)


def stage_validation():
    """Aer density-matrix agreement and transpiled two-qubit gate count."""
    from noisy_mpdo.qiskit_check import aer_observables, transpiled_cx_count

    out = {"aer_max_abs_deviation": 0.0, "cases": []}
    for n, steps, p2 in [(4, 3, 0.03), (6, 4, 0.05), (6, 4, 0.0)]:
        spec = QuenchSpec(n=n, steps=steps, seed=2)
        noise = NoiseSpec.from_p2(p2)
        dev = float(np.abs(aer_observables(spec, noise, steps) - run_exact(spec, noise).obs[-1]).max())
        out["cases"].append(dict(n=n, steps=steps, p2=p2, deviation=dev))
        out["aer_max_abs_deviation"] = max(out["aer_max_abs_deviation"], dev)
    out["transpile_n8_one_step"] = transpiled_cx_count(QuenchSpec(n=8), 1)
    return out
