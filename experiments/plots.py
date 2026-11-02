"""Figures and summary statistics (earth-tone palette) from results/*.csv."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.ticker  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

BROWN, OCHRE, TEAL, DTEAL, TERRA, SAND, CLAY = "#6b4423", "#c8962e", "#2a7f7f", "#1d4e4e", "#b5523b", "#e4d6b8", "#8c6d4f"
DELTA = 0.01
P_FIXED = 0.05  # gate error used for the size-scaling and ablation studies
N_COLOURS = (OCHRE, TEAL, TERRA, DTEAL)
N_MARKERS = "osD^"


def _final(df: pd.DataFrame, kind: str, steps: int, delta: float = DELTA) -> pd.DataFrame:
    d = df[(df.kind == kind) & (df.delta == delta) & (df.step == steps)].copy()
    d["censored"] = ~np.isfinite(d.chi_star)
    d["chi_eff"] = np.where(d.censored, d.max_cap, d.chi_star)
    return d


def _stats(vals: np.ndarray) -> tuple[float, float]:
    sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
    return float(np.mean(vals)), sd


def fit_slope(x: np.ndarray, y: np.ndarray) -> float | None:
    """Least-squares slope of log10 y against log10 x (None with fewer than 3 points)."""
    if len(x) < 3:
        return None
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])


def cost_ratio(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """log10 [(4 chi_rho*)^3 / (2 chi_psi*)^3] per (p2, step), mean over seeds (censored -> lower bound)."""
    rho = df[(df.kind == "mpdo") & (df.n == n) & (df.delta == DELTA)].copy()
    rho["chi"] = np.where(np.isfinite(rho.chi_star), rho.chi_star, rho.max_cap)
    psi = df[(df.kind == "ideal_mps") & (df.n == n) & (df.delta == DELTA)][["seed", "step", "chi_star"]]
    m = rho.merge(psi, on=["seed", "step"], suffixes=("_rho", "_psi"))
    m["log_ratio"] = 3 * np.log10(4 * m.chi / (2 * m.chi_star_psi))
    return m.groupby(["p2", "step"]).log_ratio.mean().reset_index()


def _fig_loglog(df: pd.DataFrame, steps: int, fig: Path, summary: dict) -> None:
    f, ax = plt.subplots(figsize=(6.2, 4.4))
    d_all = _final(df, "mpdo", steps)
    summary["final_chi"], summary["slopes"] = [], {}
    sweep_ns = [n for n in sorted(d_all.n.unique()) if d_all[d_all.n == n].p2.nunique() >= 2]
    for n, col, mk, off in zip(sweep_ns, N_COLOURS, N_MARKERS, ((10, 10), (10, -16))):
        xs, ms, ss, cens = [], [], [], []
        for p2, g in d_all[d_all.n == n].groupby("p2"):
            m, s = _stats(g.chi_eff.values)
            xs.append(p2), ms.append(m), ss.append(s), cens.append(bool(g.censored.any()))
            summary["final_chi"].append(dict(n=int(n), p2=float(p2), mean=m, std=s, n_seeds=int(len(g)),
                                             censored_seeds=int(g.censored.sum()), max_cap=int(g.max_cap.iloc[0])))
        xs, ms, ss, cens = map(np.array, (xs, ms, ss, cens))
        ok = ~cens
        ax.errorbar(xs[ok], ms[ok], ss[ok], color=col, marker=mk, lw=1.6, capsize=3, label=f"n = {n} qubits")
        for x, m in zip(xs[~ok], ms[~ok]):
            ax.plot(x, m, marker=mk, mfc="white", color=col, ls="none")
            ax.annotate("", (x, m * 1.5), (x, m), arrowprops=dict(arrowstyle="->", color=col))
        sl = fit_slope(xs[ok], ms[ok])
        summary["slopes"][str(int(n))] = sl
        if sl is not None:
            j = len(xs[ok]) - 2
            ax.annotate(f"fitted slope {sl:.2f}", (xs[ok][j], ms[ok][j]), textcoords="offset points", xytext=off,
                        color=col, fontsize=9, fontweight="bold")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlabel("two-qubit error probability p2 (per ZZ gate)")
    ax.set_ylabel(f"bond dimension needed, final step (error <= {DELTA})")
    ax.set_title("Noise-aware MPDO cost versus gate error")
    ax.legend(frameon=False)
    if (~d_all.chi_star.apply(np.isfinite)).any():
        ax.text(0.02, 0.03, "open marker + arrow: not resolved at the largest tested cap", transform=ax.transAxes,
                fontsize=7.5, color=CLAY)
    f.tight_layout(), f.savefig(fig / "loglog_bond_vs_gate_error.png", dpi=200), plt.close(f)


def _fig_phase(df: pd.DataFrame, fig: Path, summary: dict) -> None:
    mp = df[df.kind == "mpdo"]
    nmax = int(max(n for n in mp.n.unique() if mp[mp.n == n].p2.nunique() >= 2))
    cr = cost_ratio(df, nmax)
    pv, sv = np.sort(cr.p2.unique()), np.sort(cr.step.unique())
    z = cr.pivot(index="step", columns="p2", values="log_ratio").loc[sv, pv].values
    f, ax = plt.subplots(figsize=(6.2, 4.4))
    cmap = LinearSegmentedColormap.from_list("earth", [DTEAL, TEAL, SAND, OCHRE, TERRA, BROWN])
    cs = ax.contourf(pv, sv, np.clip(z, -3, 6), levels=np.linspace(-3, 6, 19), cmap=cmap, extend="both")
    if z.min() < 0 < z.max():
        ax.contour(pv, sv, z, levels=[0], colors="black", linewidths=2)
    ax.set_xscale("log")
    ax.set_xticks([0.01, 0.02, 0.05, 0.1])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    f.colorbar(cs, ax=ax).set_label("log10 of cost ratio, MPDO / ideal MPS")
    ax.set_xlabel("two-qubit error probability p2 (per ZZ gate)")
    ax.set_ylabel("Trotter step (circuit depth)")
    ax.set_title(f"Where simulating the noisy circuit is cheaper (n = {nmax})")
    ax.text(0.03, 0.05, "black line: equal cost; above/right of it the noisy MPDO is cheaper", transform=ax.transAxes,
            fontsize=7.5, color="#2b2118")
    f.tight_layout(), f.savefig(fig / "phase_cost_ratio.png", dpi=200), plt.close(f)
    last, cross = z[-1], None
    for a in range(len(pv) - 1):
        if last[a] > 0 >= last[a + 1]:
            w = last[a] / (last[a] - last[a + 1])
            cross = float(10 ** (np.log10(pv[a]) + w * (np.log10(pv[a + 1]) - np.log10(pv[a]))))
    summary["breakeven_p2_final_step"] = dict(n=nmax, p2=cross,
                                              log10_ratio_by_p2={f"{p:g}": float(v) for p, v in zip(pv, last)})
    first_cheap = {f"{p:g}": (int(sv[np.argmax(z[:, k] <= 0)]) if (z[:, k] <= 0).any() else None) for k, p in enumerate(pv)}
    summary["first_step_mpdo_cheaper"] = first_cheap


def _fig_waterfall(out: Path, fig: Path, summary: dict) -> None:
    ab = pd.read_csv(out / "ablation.csv")
    order = list(dict.fromkeys(ab.config))
    cap = int(ab.fixed_cap.iloc[0])
    lg = ab.assign(l=np.log10(ab.err_fixed_cap)).groupby("config", sort=False).l.agg(["mean", "std"]).loc[order]
    lev = lg["mean"].values
    f, ax = plt.subplots(figsize=(7.4, 4.5))
    base = lev.min() - 0.6  # bars are drawn from a common floor so that decades are comparable
    ax.bar(0, lev[0] - base, bottom=base, color=BROWN, yerr=lg["std"].values[0], capsize=3)
    ax.text(0, lev[0] + 0.12, f"{10 ** lev[0]:.3f}", ha="center", fontsize=9)
    for i in range(1, len(lev)):
        d = lev[i] - lev[i - 1]
        ax.bar(i, d, bottom=lev[i - 1], color=TEAL if d < 0 else TERRA)
        ax.text(i, max(lev[i], lev[i - 1]) + 0.12, f"{d:+.2f} dex", ha="center", fontsize=9)
    ax.bar(len(lev), lev[-1] - base, bottom=base, color=OCHRE, yerr=lg["std"].values[-1], capsize=4)
    ax.text(len(lev), lev[-1] + 0.12, f"{10 ** lev[-1]:.4f}", ha="center", fontsize=9)
    ax.set_xticks(range(len(lev) + 1))
    labels = [order[0]] + [c.replace(" ", chr(10), 1) for c in order[1:]] + ["full noise" + chr(10) + "model"]
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(f"log10 of observable error at bond cap {cap}")
    ax.set_title(f"Cumulative effect of each noise source (n = 8, p2 = {P_FIXED})")
    ax.set_ylim(base, lev.max() + 0.6)
    f.tight_layout(), f.savefig(fig / "waterfall_noise_source_ablation.png", dpi=200), plt.close(f)
    chi = ab.groupby("config", sort=False).chi_star.agg(["mean", "std"]).loc[order]
    summary["ablation"] = {k: dict(err_fixed_cap_mean=float(ab[ab.config == k].err_fixed_cap.mean()),
                                   err_fixed_cap_std=float(ab[ab.config == k].err_fixed_cap.std(ddof=1)),
                                   chi_star_mean=float(chi.loc[k, "mean"]), chi_star_std=float(chi.loc[k, "std"]),
                                   step_dex=(float(lev[i] - lev[i - 1]) if i else None))
                           for i, k in enumerate(order)}
    summary["ablation_fixed_cap"] = cap


def _fig_size(df: pd.DataFrame, steps: int, fig: Path, summary: dict) -> None:
    f, ax = plt.subplots(figsize=(6.2, 4.4))
    rho = _final(df, "mpdo", steps)
    rho = rho[np.isclose(rho.p2, P_FIXED)]
    psi = _final(df, "ideal_mps", steps)
    summary["size_scaling"] = []
    for kind, d, col, lab in (("ideal MPS", psi, BROWN, "ideal circuit, pure-state MPS"),
                              ("MPDO", rho, TEAL, f"noisy circuit (p2 = {P_FIXED}), MPDO")):
        ns, ms, ss = [], [], []
        for n, g in d.groupby("n"):
            m, s = _stats(g.chi_eff.values)
            ns.append(n), ms.append(m), ss.append(s)
            summary["size_scaling"].append(dict(kind=kind, n=int(n), mean=m, std=s, n_seeds=int(len(g)),
                                                censored_seeds=int(g.censored.sum()), max_cap=int(g.max_cap.iloc[0])))
        ax.plot(ns, ms, color=col, marker="o", lw=1.8, label=lab + " (mean)")
        ax.scatter(d.n, d.chi_eff, s=14, color=col, alpha=0.45, zorder=1)
    ax.set_yscale("log", base=2)
    ax.set_xlabel("number of qubits n"), ax.set_ylabel("bond dimension needed, final step")
    ax.set_title("Size scaling of the classical simulation cost (12 steps)")
    ax.legend(frameon=False)
    f.tight_layout(), f.savefig(fig / "size_scaling_bond_dimension.png", dpi=200), plt.close(f)


def _fig_gap(gap: pd.DataFrame, fig: Path, summary: dict) -> None:
    f, ax = plt.subplots(figsize=(6.0, 4.0))
    g = gap.groupby(["n", "p2"]).ideal_vs_noisy_gap.agg(["mean", "std"]).reset_index().fillna(0.0)
    for n, col in zip(sorted(g.n.unique()), N_COLOURS):
        d = g[g.n == n]
        ax.errorbar(d.p2, d["mean"], d["std"], color=col, marker="o", capsize=3, label=f"n = {n}")
    ax.axhline(DELTA, color=BROWN, ls="--", lw=1)
    ax.text(g.p2.min(), DELTA * 1.15, "shot-noise level, 10^4 shots", color=BROWN, fontsize=8)
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlabel("two-qubit error probability p2 (per ZZ gate)")
    ax.set_ylabel("max |ideal - noisy| over observables and steps")
    ax.set_title("Bias of a noise-blind classical check")
    ax.legend(frameon=False)
    f.tight_layout(), f.savefig(fig / "noise_blind_check_gap.png", dpi=200), plt.close(f)
    summary["ideal_gap"] = [dict(n=int(r.n), p2=float(r.p2), mean=float(r["mean"]), std=float(r["std"])) for _, r in g.iterrows()]


def _final_table(df: pd.DataFrame, steps: int, out: Path) -> None:
    """Final-step bond dimension (mean and std over seeds) and cost ratio for both tolerances."""
    rows = []
    for delta in (0.01, 0.03):
        d = df[(df.step == steps) & (df.delta == delta)].copy()
        d["chi"] = np.where(np.isfinite(d.chi_star), d.chi_star, d.max_cap)
        psi = d[d.kind == "ideal_mps"].set_index(["n", "seed"]).chi
        for (kind, n, p2), g in d.groupby(["kind", "n", "p2"]):
            ratio = np.nan
            if kind == "mpdo":
                ratio = float(np.mean([(4 * r.chi) ** 3 / (2 * psi[(n, r.seed)]) ** 3 for r in g.itertuples()]))
            rows.append(dict(delta=delta, kind=kind, n=n, p2=p2, chi_mean=g.chi.mean(),
                             chi_std=g.chi.std(ddof=1), n_seeds=len(g), cost_ratio_mean=ratio))
    pd.DataFrame(rows).to_csv(out / "final_step_table.csv", index=False)


def make_all(out: Path, fig: Path, steps: int) -> dict:
    """Write all figures and return the summary statistics quoted in the README."""
    df = pd.read_csv(out / "chi_star_by_step.csv")
    caps = pd.read_csv(out / "cap_errors.csv")
    gap = pd.read_csv(out / "ideal_gap.csv")
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": BROWN, "axes.labelcolor": BROWN})
    summary: dict = {"delta": DELTA, "steps": steps}
    _final_table(df, steps, out)
    _fig_loglog(df, steps, fig, summary)
    _fig_phase(df, fig, summary)
    if (out / "ablation.csv").exists():
        _fig_waterfall(out, fig, summary)
    _fig_size(df, steps, fig, summary)
    _fig_gap(gap, fig, summary)
    summary["total_search_compute_seconds"] = float(caps.seconds.sum())
    summary["seconds_by_kind_cap"] = {f"{k[0]}_cap{k[1]}": float(v) for k, v in caps.groupby(["kind", "cap"]).seconds.mean().items()}
    return summary
