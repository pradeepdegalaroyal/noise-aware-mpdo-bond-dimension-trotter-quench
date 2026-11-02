"""Three simulators of the same Trotter circuit: exact density operator, ideal MPS, noisy MPDO."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .circuits import QuenchSpec, rx, rz, rzz
from .mps import MPS
from .noise import NoiseSpec, super_1q, super_2q
from .observables import mpdo_observables, pure_observables, tensor_observables


@dataclass
class Trace:
    """Observables and bond dimension after every Trotter step (row 0 is the initial state)."""

    obs: np.ndarray  # (steps + 1, 3): ms, zz, mx
    bond: np.ndarray  # (steps + 1,) maximal bond dimension
    discarded: float = 0.0


def _vec_basis(bit: int) -> np.ndarray:
    v = np.zeros(4, dtype=complex)
    v[3 * bit] = 1.0  # |bit><bit| -> index 2*ket + bra
    return v


def _bonds(n: int):
    """Brickwork order: even bonds, then odd bonds."""
    return list(range(0, n - 1, 2)), list(range(1, n - 1, 2))


def run_exact(spec: QuenchSpec, noise: NoiseSpec) -> Trace:
    """Exact vectorised density-operator evolution (reference, n <= ~10)."""
    n = spec.n
    zz, ax, az = spec.angles()
    vec = np.ones((1,), dtype=complex)
    for b in spec.neel_bits():
        vec = np.multiply.outer(vec, _vec_basis(b))
    vec = vec.reshape((4,) * n)
    obs = [tensor_observables(vec)]
    for t in range(spec.steps):
        for parity in _bonds(n):
            for i in parity:
                s = super_2q(rzz(zz[t, i]), noise.p2).reshape(4, 4, 4, 4)
                vec = np.moveaxis(np.tensordot(s, vec, axes=([2, 3], [i, i + 1])), [0, 1], [i, i + 1])
        for q in range(n):
            s = super_1q(rz(az[t, q]) @ rx(ax[t, q]), noise.p1, noise.pd)
            vec = np.moveaxis(np.tensordot(s, vec, axes=(1, q)), 0, q)
        obs.append(tensor_observables(vec))
    return Trace(np.array(obs), np.full(spec.steps + 1, np.nan))


def run_mpdo(spec: QuenchSpec, noise: NoiseSpec, max_bond: int, cutoff: float = 1e-12) -> Trace:
    """Noise-aware matrix-product density-operator TEBD with bond cap ``max_bond``."""
    n = spec.n
    zz, ax, az = spec.angles()
    m = MPS.product([_vec_basis(b) for b in spec.neel_bits()])
    obs = [mpdo_observables(m)]
    bond = [1]
    even, odd = _bonds(n)
    for t in range(spec.steps):
        m.canonicalise()
        for i in reversed(even):
            m.apply_two(i, super_2q(rzz(zz[t, i]), noise.p2), max_bond, cutoff, "left")
        for i in odd:
            m.apply_two(i, super_2q(rzz(zz[t, i]), noise.p2), max_bond, cutoff, "right")
        for q in range(n):
            m.apply_one(q, super_1q(rz(az[t, q]) @ rx(ax[t, q]), noise.p1, noise.pd))
        obs.append(mpdo_observables(m))
        bond.append(m.max_bond())
    return Trace(np.array(obs), np.array(bond), m.discarded)


def run_ideal_mps(spec: QuenchSpec, max_bond: int, cutoff: float = 1e-12) -> Trace:
    """Noiseless pure-state MPS TEBD (the standard classical-verification baseline)."""
    n = spec.n
    zz, ax, az = spec.angles()
    m = MPS.product([np.eye(2)[b] for b in spec.neel_bits()])
    obs = [pure_observables(m)]
    bond = [1]
    even, odd = _bonds(n)
    for t in range(spec.steps):
        for i in reversed(even):
            m.apply_two(i, rzz(zz[t, i]), max_bond, cutoff, "left", normalise=True)
        for i in odd:
            m.apply_two(i, rzz(zz[t, i]), max_bond, cutoff, "right", normalise=True)
        for q in range(n):
            m.apply_one(q, rz(az[t, q]) @ rx(ax[t, q]))
        obs.append(pure_observables(m))
        bond.append(m.max_bond())
    return Trace(np.array(obs), np.array(bond), m.discarded)
