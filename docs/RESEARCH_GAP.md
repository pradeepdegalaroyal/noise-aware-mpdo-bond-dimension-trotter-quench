# Research gap

## State of the art (from the provided literature)

The literature supplied for this study holds five papers in two roles (three landmark papers and two reviews). No recent-role or adjacent-role entries were available, so the table below lists all five and nothing else is cited.

Arute et al. (2019, doi:10.1038/s41586-019-1666-5) sampled a 53-qubit superconducting processor and judged it against classical simulation of the same circuits, so the classical cost of a circuit family is the yardstick of the claim. Havlíček et al. (2019, doi:10.1038/s41586-019-0980-2) use a quantum state space as a feature space and note that kernels become expensive to estimate when the feature space grows; their argument that the feature space is only efficiently accessible on a quantum computer again rests on a classical cost. Bernien et al. (2017, doi:10.1038/nature24622) study quench dynamics of a programmable Ising-type spin model with up to 51 atoms, the physical setting that a digital Trotter circuit on a gate-based device imitates. Bauer et al. (2020, doi:10.1021/acs.chemrev.9b00829) state that arbitrary simulations cannot be expected to speed up and that areas of quantum advantage must be identified carefully for ground-state, dynamics and thermal-state problems. Gaita-Ariño et al. (2019, doi:10.1038/s41557-019-0232-y) name scalability as a central challenge of spin-based quantum computing.

### Related work table

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Arute et al., doi:10.1038/s41586-019-1666-5 | 2019 | landmark | 53-qubit random-circuit sampling, checked against classical simulation | How the classical cost of a structured dynamics circuit changes with the gate error of the device |
| Havlíček et al., doi:10.1038/s41586-019-0980-2 | 2019 | landmark | Variational classifier and kernel estimator on a superconducting processor | A noise-dependent classical-cost estimate for the circuits it runs |
| Bernien et al., doi:10.1038/nature24622 | 2017 | landmark | 51-atom programmable Ising simulator, quench dynamics and ordered phases | Digital (gate-based) emulation, noise accounting and bond-dimension cost |
| Bauer et al., doi:10.1021/acs.chemrev.9b00829 | 2020 | review | Maps quantum algorithms for ground-state, dynamics and thermal-state problems | A quantitative classical-simulability boundary for noisy dynamics circuits |
| Gaita-Ariño et al., doi:10.1038/s41557-019-0232-y | 2019 | review | Spin qubits and the scalability challenge of chemistry-based quantum hardware | Algorithm-level benchmarking against classical simulation |

## The gap

The landmark experiments are credible only relative to a classical cost, and the review asks that advantage regimes be identified carefully. In practice that cost is estimated for the *ideal* circuit with a noiseless matrix-product state (MPS). A device, however, runs a *noisy* circuit, and noise damps high-weight operator content, which should make the output cheaper to reproduce with a noise-aware tensor-network method. The quantitative question is open: at what two-qubit error rate, depth and size does a noise-aware matrix-product density operator (MPDO) become cheaper than the ideal-circuit MPS, and how much does each noise source contribute? We answer it for a Trotterised mixed-field Ising quench, the digital counterpart of the Rydberg-type dynamics of Bernien et al., at simulation scale (6-12 qubits, 12 Trotter steps). All statements are simulation results on a laptop and do not extrapolate to 50+ qubits.

## Hypotheses and success metric

The quantity measured is the bond dimension chi* that a truncated time-evolving block decimation (TEBD) simulation needs for every tracked local observable (staggered magnetisation, nearest-neighbour ZZ correlator, mean X) to stay within 0.01 of the exact answer at all steps up to the stated depth. The tolerance 0.01 equals the statistical error of 10^4 shots. The cost proxy of one two-site update is (d chi)^3, with local dimension d = 2 for the MPS and d = 4 for the MPDO.

- **H1 (hardware-relevant break-even).** For n = 8 the final-step cost ratio (4 chi*_MPDO)^3 / (2 chi*_MPS)^3 drops below 1 at a two-qubit error p2 <= 0.01, the range of present-day devices. Falsified if the break-even p2 is above 0.01 or absent.
- **H2 (power law).** chi*_MPDO falls as p2^(-1) on a log-log plot. Falsified if the fitted slope differs from -1 by more than the seed-to-seed spread.
- **H3 (size saturation).** At fixed p2 the MPDO bond dimension saturates with n while the ideal MPS needs 2^(n/2). Falsified if chi*_MPDO keeps growing with n over the tested sizes.

## Novelty

This study measures, with a validated noise-aware MPDO engine and fixed-tolerance bond-dimension searches, where the noisy-circuit and ideal-circuit classical costs cross for a Trotterised Ising quench, and attributes the saving to individual noise sources.
