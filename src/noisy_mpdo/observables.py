"""Local observables tracked along the quench: staggered magnetisation, ZZ correlator, X mean."""

from __future__ import annotations

import numpy as np

from .circuits import X, Z
from .mps import MPS

OBS_NAMES = ("ms", "zz", "mx")

# weight vectors w[2*ket + bra] such that Tr(O rho) = sum w * vec(rho)
W_I = np.array([1, 0, 0, 1], dtype=complex)
W_Z = np.array([1, 0, 0, -1], dtype=complex)
W_X = np.array([0, 1, 1, 0], dtype=complex)


def _combine(z: np.ndarray, zz: np.ndarray, x: np.ndarray) -> np.ndarray:
    n = len(z)
    sign = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    return np.array([np.dot(sign, z) / n, zz.mean(), x.mean()])


def mpdo_observables(m: MPS) -> np.ndarray:
    """(ms, zz, mx) of a vectorised density operator, divided by its (truncation-drifted) trace."""
    n = m.n
    left = [np.ones(1, dtype=complex)]
    for a in m.t:
        left.append(left[-1] @ np.tensordot(a, W_I, axes=(1, 0)))
    right = [np.ones(1, dtype=complex)]
    for a in reversed(m.t):
        right.append(np.tensordot(a, W_I, axes=(1, 0)) @ right[-1])
    right = right[::-1]  # right[i] = contraction of sites i..n-1
    trace = left[n][0]
    z = np.empty(n)
    x = np.empty(n)
    for i, a in enumerate(m.t):
        z[i] = (left[i] @ np.tensordot(a, W_Z, axes=(1, 0)) @ right[i + 1]).real
        x[i] = (left[i] @ np.tensordot(a, W_X, axes=(1, 0)) @ right[i + 1]).real
    zz = np.empty(n - 1)
    for i in range(n - 1):
        v = left[i] @ np.tensordot(m.t[i], W_Z, axes=(1, 0))
        v = v @ np.tensordot(m.t[i + 1], W_Z, axes=(1, 0))
        zz[i] = (v @ right[i + 2]).real
    return _combine(z, zz, x) / trace.real


def pure_observables(m: MPS) -> np.ndarray:
    """(ms, zz, mx) of a pure-state MPS using left and right environments."""
    n = m.n
    left = [np.ones((1, 1), dtype=complex)]
    for a in m.t:
        left.append(np.einsum("xy,xsr,ysq->rq", left[-1], a.conj(), a))
    right = [np.ones((1, 1), dtype=complex)]
    for a in reversed(m.t):
        right.append(np.einsum("xsr,ysq,rq->xy", a.conj(), a, right[-1]))
    right = right[::-1]
    norm = left[n][0, 0].real
    z = np.empty(n)
    x = np.empty(n)
    for i, a in enumerate(m.t):
        for op, out in ((Z, z), (X, x)):
            oa = np.einsum("ps,lsr->lpr", op, a)
            out[i] = np.einsum("xy,xsr,ysq,rq->", left[i], a.conj(), oa, right[i + 1]).real
    zz = np.empty(n - 1)
    for i in range(n - 1):
        a, b = m.t[i], m.t[i + 1]
        za = np.einsum("ps,lsr->lpr", Z, a)
        zb = np.einsum("ps,lsr->lpr", Z, b)
        e = np.einsum("xy,xsr,ysq->rq", left[i], a.conj(), za)
        zz[i] = np.einsum("xy,xsr,ysq,rq->", e, b.conj(), zb, right[i + 2]).real
    return _combine(z, zz, x) / norm


def tensor_observables(vec: np.ndarray) -> np.ndarray:
    """(ms, zz, mx) of a full vectorised density operator stored as an (4,)*n tensor."""
    n = vec.ndim

    def contract(ws: dict[int, np.ndarray]) -> float:
        t = vec
        for site in range(n - 1, -1, -1):
            t = np.tensordot(t, ws.get(site, W_I), axes=(site, 0))
        return float(t.real)

    tr = contract({})
    z = np.array([contract({i: W_Z}) for i in range(n)])
    x = np.array([contract({i: W_X}) for i in range(n)])
    zz = np.array([contract({i: W_Z, i + 1: W_Z}) for i in range(n - 1)])
    return _combine(z, zz, x) / tr
