"""Cross-check against Qiskit Aer's density-matrix simulator and transpilation-depth accounting."""

from __future__ import annotations

import numpy as np

from .circuits import QuenchSpec, X, Z, to_qiskit
from .noise import NoiseSpec
from .observables import _combine


def aer_noise_model(noise: NoiseSpec):
    """Aer noise model equivalent to ``NoiseSpec`` (Pauli-error probabilities -> depolarising parameters)."""
    from qiskit_aer.noise import NoiseModel, depolarizing_error, pauli_error

    nm = NoiseModel(basis_gates=["rzz", "rx", "rz", "id", "x"])
    if noise.p2 > 0:
        nm.add_all_qubit_quantum_error(depolarizing_error(16 * noise.p2 / 15, 2), ["rzz"])
    if noise.p1 > 0:
        nm.add_all_qubit_quantum_error(depolarizing_error(4 * noise.p1 / 3, 1), ["rz"])
    if noise.pd > 0:
        nm.add_all_qubit_quantum_error(pauli_error([("Z", noise.pd), ("I", 1 - noise.pd)]), ["id"])
    return nm


def aer_observables(spec: QuenchSpec, noise: NoiseSpec, steps: int) -> np.ndarray:
    """(ms, zz, mx) after ``steps`` Trotter steps from Aer density-matrix simulation."""
    from qiskit import transpile
    from qiskit_aer import AerSimulator

    qc = to_qiskit(spec, steps)
    qc.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=aer_noise_model(noise))
    rho = np.asarray(sim.run(transpile(qc, sim, optimization_level=0)).result().data()["density_matrix"])
    n = spec.n

    def ev(ops: dict[int, np.ndarray]) -> float:
        # Qiskit is little-endian: qubit 0 is the least-significant tensor factor.
        full = np.array([[1.0]], dtype=complex)
        for q in range(n - 1, -1, -1):
            full = np.kron(full, ops.get(q, np.eye(2)))
        return float(np.trace(full @ rho).real)

    z = np.array([ev({i: Z}) for i in range(n)])
    x = np.array([ev({i: X}) for i in range(n)])
    zz = np.array([ev({i: Z, i + 1: Z}) for i in range(n - 1)])
    return _combine(z, zz, x)


def transpiled_cx_count(spec: QuenchSpec, steps: int = 1) -> dict[str, int]:
    """Two-qubit gate count and depth of ``steps`` Trotter steps on a line after transpilation to {cx, rz, sx, x}."""
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap

    qc = to_qiskit(spec, steps)
    tq = transpile(
        qc,
        basis_gates=["cx", "rz", "sx", "x"],
        coupling_map=CouplingMap.from_line(spec.n),
        optimization_level=1,
        seed_transpiler=0,
    )
    ops = tq.count_ops()
    return {"cx": int(ops.get("cx", 0)), "rzz_in_circuit": steps * (spec.n - 1), "depth": int(tq.depth())}
