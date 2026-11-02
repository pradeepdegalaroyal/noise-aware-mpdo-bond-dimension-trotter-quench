"""Pauli noise model and the superoperators used by the density-operator simulators.

A density operator is stored vectorised: each qubit carries a four-level index
``s = 2 * ket + bra``, so an operator ``rho -> K rho K^dagger`` acts as ``kron(K, conj(K))``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np

from .circuits import PAULIS, Z


@dataclass(frozen=True)
class NoiseSpec:
    """Per-gate Pauli noise.

    p2: total probability of a non-identity two-qubit Pauli error after each ZZ gate.
    p1: same for a single-qubit depolarising error after each RX/RZ layer.
    pd: Z-flip probability on every qubit once per Trotter step (idle dephasing proxy).
    """

    p2: float = 0.0
    p1: float = 0.0
    pd: float = 0.0

    @classmethod
    def from_p2(cls, p2: float, p1_ratio: float = 0.1, pd_ratio: float = 0.2) -> "NoiseSpec":
        """Hardware-style ratios: p1 = 0.1 p2 and idle dephasing pd = 0.2 p2."""
        return cls(p2=p2, p1=p1_ratio * p2, pd=pd_ratio * p2)

    def is_noiseless(self) -> bool:
        return self.p2 == 0 and self.p1 == 0 and self.pd == 0


def _kron_super(k: np.ndarray) -> np.ndarray:
    return np.kron(k, k.conj())


def super_1q(unitary: np.ndarray | None = None, p1: float = 0.0, pd: float = 0.0) -> np.ndarray:
    """4x4 superoperator: unitary, then 1q depolarising(p1), then dephasing(pd)."""
    s = np.eye(4, dtype=complex) if unitary is None else _kron_super(unitary)
    if p1 > 0:
        dep = (1 - p1) * np.eye(4, dtype=complex) + (p1 / 3) * sum(_kron_super(P) for P in PAULIS[1:])
        s = dep @ s
    if pd > 0:
        s = ((1 - pd) * np.eye(4, dtype=complex) + pd * _kron_super(Z)) @ s
    return s


def _pair_order(s: np.ndarray) -> np.ndarray:
    """Permute a 16x16 superoperator from ((a_i a_j),(b_i b_j)) to ((a_i b_i),(a_j b_j)) ordering."""
    t = s.reshape(2, 2, 2, 2, 2, 2, 2, 2).transpose(0, 2, 1, 3, 4, 6, 5, 7)
    return t.reshape(16, 16)


def super_2q(unitary: np.ndarray, p2: float = 0.0) -> np.ndarray:
    """16x16 superoperator on two vectorised qubits: unitary followed by 2q depolarising(p2)."""
    s = _kron_super(unitary)
    if p2 > 0:
        dep = (1 - p2) * np.eye(16, dtype=complex)
        for a, b in product(range(4), repeat=2):
            if a == 0 and b == 0:
                continue
            dep = dep + (p2 / 15) * _kron_super(np.kron(PAULIS[a], PAULIS[b]))
        s = dep @ s
    return _pair_order(s)


def pure_pair(unitary: np.ndarray) -> np.ndarray:
    """Plain 4x4 two-qubit unitary for the pure-state MPS."""
    return unitary
