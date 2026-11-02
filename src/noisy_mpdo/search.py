"""Minimal bond dimension needed to meet an observable-error tolerance, with censoring."""

from __future__ import annotations

import numpy as np


def chi_star_curve(caps: list[int], errors: np.ndarray, delta: float) -> np.ndarray:
    """Smallest bond cap whose running-maximum error stays below ``delta``, for every step.

    ``errors[k, t]`` is the largest observable error at step ``t`` for cap ``caps[k]``
    (caps ascending). Steps no tested cap resolves are returned as ``inf`` (censored).
    """
    running = np.maximum.accumulate(np.asarray(errors, dtype=float), axis=1)
    star = np.full(running.shape[1], np.inf)
    for k in range(len(caps) - 1, -1, -1):
        star[running[k] <= delta] = caps[k]
    return star


def obs_error(obs: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Largest absolute error over the tracked observables at every step."""
    return np.abs(obs - ref).max(axis=1)


def search_caps(run, ref_obs: np.ndarray, caps: list[int], delta: float):
    """Run ``run(cap)`` for ascending caps until the final running error is <= ``delta``.

    Returns (caps_used, errors[k, step], wall_seconds[k], max_bond_reached[k]).
    """
    import time

    used, errs, secs, reached = [], [], [], []
    for cap in caps:
        t0 = time.perf_counter()
        tr = run(cap)
        secs.append(time.perf_counter() - t0)
        used.append(cap)
        errs.append(obs_error(tr.obs, ref_obs))
        reached.append(int(tr.bond.max()))
        if np.max(errs[-1]) <= delta:
            break
    return used, np.array(errs), np.array(secs), np.array(reached)
