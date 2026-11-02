"""Matrix-product state with time-evolving block decimation (TEBD) updates.

The same class stores a pure state (local dimension 2) or a vectorised density operator
(local dimension 4, a matrix-product density operator, MPDO). Tensors have shape
(left bond, physical, right bond) and the orthogonality centre is tracked explicitly so that
every truncation is the optimal two-norm truncation.
"""

from __future__ import annotations

import numpy as np


class MPS:
    """Open-boundary matrix-product state with truncated two-site updates."""

    def __init__(self, tensors: list[np.ndarray], center: int = 0) -> None:
        self.t = tensors
        self.n = len(tensors)
        self.center = center
        self.discarded = 0.0  # accumulated discarded weight (sum of per-gate truncation errors)

    @classmethod
    def product(cls, local_states: list[np.ndarray]) -> "MPS":
        return cls([np.asarray(v, dtype=complex).reshape(1, -1, 1) for v in local_states], 0)

    def bond_dims(self) -> list[int]:
        return [self.t[i].shape[2] for i in range(self.n - 1)]

    def max_bond(self) -> int:
        return max(self.bond_dims()) if self.n > 1 else 1

    # ---- canonical form -------------------------------------------------------------
    def _shift_right(self, i: int) -> None:
        a = self.t[i]
        l, d, r = a.shape
        q, rr = np.linalg.qr(a.reshape(l * d, r))
        self.t[i] = q.reshape(l, d, -1)
        self.t[i + 1] = np.tensordot(rr, self.t[i + 1], axes=(1, 0))

    def _shift_left(self, i: int) -> None:
        a = self.t[i]
        l, d, r = a.shape
        q, rr = np.linalg.qr(a.reshape(l, d * r).T)
        self.t[i] = q.T.reshape(-1, d, r)
        self.t[i - 1] = np.tensordot(self.t[i - 1], rr.T, axes=(2, 0))

    def move_center(self, j: int) -> None:
        while self.center < j:
            self._shift_right(self.center)
            self.center += 1
        while self.center > j:
            self._shift_left(self.center)
            self.center -= 1

    def canonicalise(self) -> None:
        """Full left-to-right sweep; afterwards the centre sits on the last site."""
        for i in range(self.n - 1):
            self._shift_right(i)
        self.center = self.n - 1

    # ---- gates ------------------------------------------------------------------------
    def apply_one(self, i: int, m: np.ndarray) -> None:
        self.t[i] = np.einsum("ps,lsr->lpr", m, self.t[i])

    def apply_two(
        self,
        i: int,
        m: np.ndarray,
        max_bond: int,
        cutoff: float = 1e-12,
        direction: str = "right",
        normalise: bool = False,
    ) -> None:
        """Apply a (d^2 x d^2) matrix on sites (i, i+1), then truncate the bond by SVD."""
        self.move_center(i)
        a, b = self.t[i], self.t[i + 1]
        l, d, _ = a.shape
        r = b.shape[2]
        theta = np.tensordot(a, b, axes=(2, 0))  # l, d, d, r
        theta = np.einsum("qpsu,lsur->lqpr", m.reshape(d, d, d, d), theta)
        u, s, vh = np.linalg.svd(theta.reshape(l * d, d * r), full_matrices=False)
        tot = float(np.sum(s**2))
        k = int(min(max_bond, np.sum(s > cutoff * s[0]))) if s[0] > 0 else 1
        k = max(k, 1)
        if tot > 0:
            self.discarded += float(np.sum(s[k:] ** 2)) / tot
        s = s[:k]
        if normalise:
            s = s / np.linalg.norm(s)
        u, vh = u[:, :k], vh[:k]
        if direction == "right":
            self.t[i] = u.reshape(l, d, k)
            self.t[i + 1] = (s[:, None] * vh).reshape(k, d, r)
            self.center = i + 1
        else:
            self.t[i] = (u * s[None, :]).reshape(l, d, k)
            self.t[i + 1] = vh.reshape(k, d, r)
            self.center = i

    # ---- contractions -----------------------------------------------------------------
    def expect_vec(self, weights: list[np.ndarray]) -> complex:
        """Contract a vectorised density operator with one weight vector per site."""
        env = np.ones(1, dtype=complex)
        for a, w in zip(self.t, weights):
            env = env @ np.tensordot(a, w, axes=(1, 0))
        return complex(env[0])

    def expect_pure(self, ops: dict[int, np.ndarray]) -> complex:
        """<psi| prod_i O_i |psi> for a pure state (not necessarily normalised)."""
        env = np.ones((1, 1), dtype=complex)
        for i, a in enumerate(self.t):
            o = ops.get(i)
            b = a if o is None else np.einsum("ps,lsr->lpr", o, a)
            env = np.einsum("xy,xsr,ysq->rq", env, a.conj(), b)
        return complex(env[0, 0])
