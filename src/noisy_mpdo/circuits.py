"""Trotterised mixed-field Ising quench circuit (Rydberg-type chain) and its Qiskit export."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = (I2, X, Y, Z)


def rx(theta: float) -> np.ndarray:
    """Single-qubit rotation exp(-i theta X / 2)."""
    return np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * X


def rz(theta: float) -> np.ndarray:
    """Single-qubit rotation exp(-i theta Z / 2)."""
    return np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * Z


def rzz(theta: float) -> np.ndarray:
    """Two-qubit rotation exp(-i theta Z Z / 2)."""
    return np.diag(np.exp(-0.5j * theta * np.array([1, -1, -1, 1]))).astype(complex)


@dataclass(frozen=True)
class QuenchSpec:
    """One brickwork Trotter step is ZZ(even bonds), ZZ(odd bonds), then RX and RZ on every qubit.

    The circuit starts from the Neel product state |0101...>. Angles carry a Gaussian
    calibration jitter (relative standard deviation ``jitter``) that is fixed by ``seed``.
    The default fields are a generic non-integrable (fast-scrambling) choice.
    """

    n: int
    steps: int = 16
    dt: float = 0.4
    J: float = 1.0
    hx: float = 0.9045
    hz: float = 0.809
    jitter: float = 0.1
    seed: int = 0

    def angles(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return rotation angles (zz[step, bond], rx[step, qubit], rz[step, qubit])."""
        rng = np.random.default_rng(self.seed)

        def jit(shape: tuple[int, int]) -> np.ndarray:
            return 1.0 + self.jitter * rng.standard_normal(shape)

        zz = 2 * self.J * self.dt * jit((self.steps, self.n - 1))
        ax = 2 * self.hx * self.dt * jit((self.steps, self.n))
        az = 2 * self.hz * self.dt * jit((self.steps, self.n))
        return zz, ax, az

    def neel_bits(self) -> list[int]:
        """Initial computational-basis bit string 0101..."""
        return [i % 2 for i in range(self.n)]


def to_qiskit(spec: QuenchSpec, steps: int | None = None):
    """Build the same circuit as a Qiskit ``QuantumCircuit`` (``id`` gates mark the idle-dephasing slot)."""
    from qiskit import QuantumCircuit

    steps = spec.steps if steps is None else steps
    zz, ax, az = spec.angles()
    qc = QuantumCircuit(spec.n)
    for q, b in enumerate(spec.neel_bits()):
        if b:
            qc.x(q)
    for t in range(steps):
        for parity in (0, 1):
            for i in range(parity, spec.n - 1, 2):
                qc.rzz(float(zz[t, i]), i, i + 1)
        for q in range(spec.n):
            qc.rx(float(ax[t, q]), q)
            qc.rz(float(az[t, q]), q)
        for q in range(spec.n):
            qc.id(q)
    return qc
