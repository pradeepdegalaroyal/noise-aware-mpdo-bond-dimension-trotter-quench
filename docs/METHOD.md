# Method

## Circuit and noise model

One Trotter step of the mixed-field Ising chain (open boundary, n qubits, brickwork) is: a two-qubit ZZ rotation (RZZ) of angle 2 J dt on even bonds, RZZ on odd bonds, then an X rotation (RX) of angle 2 hx dt and a Z rotation (RZ) of angle 2 hz dt on every qubit, with coupling J = 1, transverse field hx = 0.9045, longitudinal field hz = 0.809 and Trotter step dt = 0.4 (a generic non-integrable choice that scrambles quickly). The initial state is the Néel product state |0101...>. Every angle carries a Gaussian calibration jitter of 10 % relative standard deviation drawn from the seed, so seeds change the circuit instance. The circuit itself, not the continuous-time Hamiltonian, is the target, so Trotter error is not part of the question.

Noise is Pauli noise attached to gates: a two-qubit depolarising error with total probability p2 after every RZZ, a single-qubit depolarising error p1 = 0.1 p2 after each RX/RZ layer, and a Z-flip probability pd = 0.2 p2 per qubit per step as an idle-dephasing proxy. Hardware reading: one RZZ transpiles to two controlled-NOT (CNOT) gates (checked with Qiskit on a line), so p2 is about twice the per-CNOT error.

## Simulators

1. **Exact density operator (reference, n <= 10).** The vectorised density operator is stored as a rank-n tensor with local dimension 4 and updated with 16 x 16 superoperators.
2. **Ideal matrix-product state (MPS) baseline (classical standard).** Noiseless pure-state MPS evolved with time-evolving block decimation (TEBD), bond cap chi, two-norm truncation in canonical form. This is how the classical cost of an ideal circuit is usually assessed.
3. **Noise-aware matrix-product density operator (MPDO, the method).** The vectorised density operator stored as an MPS with local dimension 4. Gate and noise are combined into one superoperator per bond, and the same TEBD engine truncates the bond to cap chi. Observables are traces against product weight vectors and are divided by the (truncation-drifted) trace.
4. **Standard quantum baseline.** The plain noisy Trotter circuit with no mitigation, i.e. what the device would return. Its output is the MPDO target; the noise-blind ideal MPS is also scored against it (figure `noise_blind_check_gap.png`).

The superoperator construction is validated against the Qiskit Aer density-matrix simulator (Aer is Qiskit's simulator package) with an equivalent noise model (`results/validation.json`).

## Metrics

- chi*(t): smallest bond cap in {2, 4, ..., cap_max} whose running-maximum observable error up to step t is at most delta = 0.01 (0.03 is also reported as a tolerance check). If no tested cap meets delta, the value is *censored* and reported as "> cap_max".
- Cost ratio: (4 chi*_MPDO)^3 / (2 chi*_MPS)^3, the ratio of per-update singular-value-decomposition costs.
- Noise-blind gap: maximum absolute difference between the ideal and the noisy observables.

## Ablation

Cumulative noise-source ablation at n = 8, p2 = 0.05, final step: ideal circuit (chi*_MPS squared, the exact bond dimension of |psi><psi|), then add two-qubit depolarising noise, then idle dephasing, then single-qubit depolarising noise. Each bar is the change in log2 chi*.

## Sweeps and seeds

- Noise sweep: p2 in {0.01, 0.02, 0.035, 0.05, 0.07, 0.1}, n = 6 and 8, exact reference, 12 steps, caps up to 64.
- Size sweep: n = 6, 8, 10, 12 at p2 = 0.05 (n = 12 uses an MPDO reference with cap 64, since an exact 12-qubit density operator is too slow for the time box; caps up to 32 are tested).
- Three seeds (0, 1, 2) everywhere; tables give mean +- sample standard deviation over seeds.

## Reproduction

`python experiments/run_all.py` runs all stages and writes `results/` and `figures/`; `--quick` is a smoke test.
