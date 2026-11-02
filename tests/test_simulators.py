"""Fast correctness tests: superoperators, simulator agreement, truncation behaviour, Qiskit cross-check."""

from __future__ import annotations

import numpy as np
import pytest

from noisy_mpdo import NoiseSpec, QuenchSpec, run_exact, run_ideal_mps, run_mpdo
from noisy_mpdo.circuits import rx, rz, rzz
from noisy_mpdo.noise import super_1q, super_2q
from noisy_mpdo.search import chi_star_curve


def test_superoperators_are_trace_preserving():
    tr1 = np.array([1, 0, 0, 1])
    s1 = super_1q(rz(0.3) @ rx(0.7), p1=0.05, pd=0.02)
    assert np.allclose(tr1 @ s1, tr1)
    tr2 = np.kron(tr1, tr1)
    s2 = super_2q(rzz(0.9), p2=0.04)
    assert np.allclose(tr2 @ s2, tr2)


def test_depolarising_damps_single_qubit_z():
    # Z -> (1 - 16 p / 15) Z on each qubit of a two-qubit depolarising channel
    p = 0.06
    s = super_2q(np.eye(4), p2=p)
    z_vec = np.kron([1, 0, 0, -1], [1, 0, 0, 1]).astype(complex)  # Z (x) I as bra-weights
    rho = np.kron([1, 0, 0, 0], [1, 0, 0, 0]).astype(complex)  # |00><00|
    out = z_vec @ (s @ rho)
    assert out == pytest.approx(1 - 16 * p / 15)


@pytest.mark.parametrize("n", [4, 6])
def test_untruncated_mpdo_matches_exact_density_operator(n):
    spec = QuenchSpec(n=n, steps=4, seed=3)
    noise = NoiseSpec.from_p2(0.03)
    ref = run_exact(spec, noise)
    got = run_mpdo(spec, noise, max_bond=4 ** (n // 2))
    assert np.abs(ref.obs - got.obs).max() < 1e-10


def test_ideal_mps_matches_noiseless_exact_and_bond_growth():
    spec = QuenchSpec(n=6, steps=5, seed=1)
    ref = run_exact(spec, NoiseSpec())
    got = run_ideal_mps(spec, max_bond=8)
    assert np.abs(ref.obs - got.obs).max() < 1e-10
    assert got.bond.max() <= 8 and got.bond[0] == 1


def test_truncation_error_shrinks_with_bond_cap():
    spec = QuenchSpec(n=6, steps=6, seed=2)
    noise = NoiseSpec.from_p2(0.03)
    ref = run_exact(spec, noise)
    err = [np.abs(run_mpdo(spec, noise, c).obs - ref.obs).max() for c in (2, 4, 16)]
    assert err[0] > err[1] > err[2]


def test_chi_star_curve_prefix_logic():
    caps = [2, 4, 8]
    # error per step for each cap (row = cap); cap 4 fails at step 2, cap 8 is fine everywhere
    errs = np.array([[0, 0.5, 0.5], [0, 0.0, 0.2], [0, 0.0, 0.0]])
    star = chi_star_curve(caps, errs, delta=0.1)
    assert list(star) == [2, 4, 8]
    errs[2, 2] = 0.3
    assert np.isinf(chi_star_curve(caps, errs, delta=0.1)[2])


def test_qiskit_aer_density_matrix_agrees():
    pytest.importorskip("qiskit_aer")
    from noisy_mpdo.qiskit_check import aer_observables

    spec = QuenchSpec(n=4, steps=3, seed=2)
    noise = NoiseSpec.from_p2(0.03)
    ref = run_exact(spec, noise).obs[-1]
    assert np.abs(aer_observables(spec, noise, 3) - ref).max() < 1e-9
