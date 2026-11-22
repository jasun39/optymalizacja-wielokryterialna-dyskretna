"""Computes the experiments of the scalarization report and saves figures (fig/) and result tables (results.md).

Uses the same code as the Streamlit page (owd/scalarization.py) with the page's default settings.
Requires matplotlib (only for this script):  pip install matplotlib
Run from the repository root:  python reports/cw3_skalaryzacja/generate_figures.py
"""
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))  # repository root

from owd.scalarization import (PRESETS, Problem, distance_to_ideal, first_occurrence, front_spline_segments,  # noqa: E402
                           parse_function, payoff_table, run_method, spline_segments, surface_triangles)

FIG = ROOT / "fig"
FIG.mkdir(exist_ok=True)
COLORS = {"S1": "#1f77b4", "S2": "#2ca02c", "S3 (p=2)": "#d62728", "S3 (p=∞)": "#9467bd"}
SEED, STARTS, DELTA, EPS_CRITERION = 0, 3, 0.1, 0

EXPERIMENTS = [
    # (id, title, n, preset, domain, step)
    ("a_convex", "Wariant (a) – front wypukły, U = [-1, 1]²", 2, "Wypukły (domyślny)", "square", 0.05),
    ("a_fonseca", "Wariant (a) – front niewypukły (Fonseca–Fleming), U = [-1, 1]²", 2,
     "Niewypukły front (Fonseca–Fleming)", "square", 0.05),
    ("a_cubic_ball", "Wariant (a) – wielomiany 3. stopnia, U = koło", 2, "Wielomiany 3. stopnia", "ball", 0.05),
    ("b_mixed", "Wariant (b) – jedno kryterium liniowe, U = [-1, 1]³", 3, "Mieszany – jedno liniowe (domyślny)",
     "square", 0.1),
    ("b_mixed_ball", "Wariant (b) – jedno kryterium liniowe, U = kula", 3, "Mieszany – jedno liniowe (domyślny)",
     "ball", 0.1),
]
VARIANTS = [("S1", "S1", None), ("S2", "S2", None), ("S3 (p=2)", "S3", 2.0), ("S3 (p=∞)", "S3", math.inf)]


def solve(n, preset, domain, step):
    expressions = [parse_function(text, n) for text in PRESETS[n][preset]]
    problem = Problem(expressions, n, domain)
    rng = np.random.default_rng(SEED)
    starts = list(problem.random_points(STARTS, rng))
    table, ideal, nadir = payoff_table(problem, starts)
    reference = ideal - DELTA * (nadir - ideal)
    results = {}
    for label, method, p in VARIANTS:
        settings = {"step": step, "eps_criterion": EPS_CRITERION, "p": p or 2.0, "reference": reference}
        results[label] = run_method(method, problem, settings, ideal, nadir, starts)
    return problem, table, ideal, nadir, reference, results


def arrays(result):
    solved = result.solved
    F = np.array([run.F for run in solved]).reshape(-1, 0 if not solved else len(solved[0].F))
    X = np.array([run.x for run in solved])
    grid = np.array([run.grid for run in solved], dtype=float)
    return F, X, grid, result.nondominated


def compromise(F, mask, ideal, nadir):
    """Index (into F) of the nondominated point closest to the ideal point (normalized, p = 2)."""
    candidates = np.flatnonzero(mask)
    return candidates[np.argmin(distance_to_ideal(F[candidates], ideal, nadir, 2))]


def plot_2d(exp_id, title, problem, ideal, nadir, reference, results):
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.4))
    for column, (label, result) in enumerate(results.items()):
        F, X, _, mask = arrays(result)
        color = COLORS[label]
        ax_f, ax_x = axes[0, column], axes[1, column]
        if problem.domain == "square":
            ax_x.add_patch(plt.Rectangle((-1, -1), 2, 2, fill=False, ls="--", color="grey"))
        else:
            ax_x.add_patch(plt.Circle((0, 0), 1, fill=False, ls="--", color="grey"))
        for piece in front_spline_segments(F[mask]):
            ax_f.plot(piece[:, 0], piece[:, 1], color=color, lw=1.3)
        ax_f.scatter(F[mask, 0], F[mask, 1], s=18, color=color, zorder=3, label="niezdominowane")
        if (~mask).any():
            ax_f.scatter(F[~mask, 0], F[~mask, 1], s=30, color="black", zorder=3, marker="x", label="zdominowane")
        order = np.lexsort((F[mask, 1], F[mask, 0]))
        for piece in spline_segments(X[mask][order]):
            ax_x.plot(piece[:, 0], piece[:, 1], color=color, lw=1.3)
        ax_x.scatter(X[mask, 0], X[mask, 1], s=18, color=color, zorder=3)
        best = compromise(F, mask, ideal, nadir)
        ax_f.scatter(*F[best], s=140, facecolors="none", edgecolors="#e6a800", lw=2, zorder=5,
                     label="najbliższy idealnemu")
        ax_x.scatter(*X[best], s=140, facecolors="none", edgecolors="#e6a800", lw=2, zorder=5)
        ax_f.scatter(*ideal, marker="*", s=160, color="gold", edgecolors="black", label="punkt idealny", zorder=4)
        if label.startswith("S3"):
            ax_f.scatter(*reference, marker="+", s=110, color="#8c564b", label="punkt r", zorder=4)
        ax_f.set(xlabel="F1", ylabel="F2", title=f"{label}: $F_P(U)$")
        ax_x.set(xlabel="x1", ylabel="x2", title=f"{label}: $P(U,F)$", aspect="equal",
                 xlim=(-1.1, 1.1), ylim=(-1.1, 1.1))
        ax_f.legend(fontsize=7, loc="upper right")
        for ax in (ax_f, ax_x):
            ax.grid(alpha=0.3)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(FIG / f"{exp_id}.png", dpi=120)
    plt.close(fig)


def plot_3d(exp_id, title, problem, ideal, nadir, reference, results):
    fig = plt.figure(figsize=(16, 8.5))
    for column, (label, result) in enumerate(results.items()):
        F, X, grid, mask = arrays(result)
        color = COLORS[label]
        for row, (values, names, space) in enumerate(((F, ("G1", "G2", "G3"), "ocen"),
                                                       (X, ("x1", "x2", "x3"), "decyzji"))):
            ax = fig.add_subplot(2, 4, row * 4 + column + 1, projection="3d")
            nd = values[mask]
            if len(nd) >= 3:
                triangles = surface_triangles(grid[mask], nd)
                if len(triangles):
                    ax.plot_trisurf(nd[:, 0], nd[:, 1], nd[:, 2], triangles=triangles, color=color, alpha=0.25,
                                    edgecolor=color, linewidth=0.3)
            ax.scatter(nd[:, 0], nd[:, 1], nd[:, 2], s=8, color=color)
            if (~mask).any():
                ax.scatter(*values[~mask].T, s=14, color="black", marker="x")
            if row == 0:
                ax.scatter(*ideal, marker="*", s=120, color="gold", edgecolors="black")
            best = compromise(F, mask, ideal, nadir)
            ax.scatter(*values[best], s=90, facecolors="none", edgecolors="#e6a800", lw=2, depthshade=False)
            ax.set(xlabel=names[0], ylabel=names[1], zlabel=names[2])
            ax.set_box_aspect(None, zoom=0.85)
            ax.set_title(f"{label} – przestrzeń {space}", fontsize=10)
            ax.tick_params(labelsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(FIG / f"{exp_id}.png", dpi=110)
    plt.close(fig)


def fmt(vector):
    return "(" + ", ".join(f"{v:.3f}" for v in vector) + ")"


def summary(exp_id, title, problem, table, ideal, nadir, reference, results):
    lines = [f"### {title} [{exp_id}]", "", "Tablica wypłat (wiersz i = F w minimum kryterium i):", ""]
    m = problem.n_criteria
    lines += ["| min | " + " | ".join(f"F{j + 1}" for j in range(m)) + " |", "|---" * (m + 1) + "|"]
    for i in range(m):
        lines.append(f"| F{i + 1} | " + " | ".join(f"{v:.4f}" for v in table[i]) + " |")
    lines += ["", f"ideal = {fmt(ideal)}, nadir = {fmt(nadir)}, r = {fmt(reference)}", ""]
    lines += ["| metoda | zadań | rozwiązanych | niezdominowanych | unikalnych niezd. | czas [s] | "
              "najbliższy idealnemu: x | F(x) | d₂ |", "|---" * 9 + "|"]
    scale = float(np.max(nadir - ideal))
    for label, result in results.items():
        F, X, _, mask = arrays(result)
        first = first_occurrence(F, scale)
        unique_nd = int(np.sum(mask & (first == np.arange(len(F)))))
        d = distance_to_ideal(F[mask], ideal, nadir, 2)
        best = int(np.argmin(d))
        lines.append(f"| {label} | {len(result.runs)} | {len(F)} | {int(mask.sum())} | {unique_nd} | "
                     f"{result.time_s:.2f} | {fmt(X[mask][best])} | {fmt(F[mask][best])} | {d[best]:.4f} |")
    if exp_id == "a_convex":
        lines += ["", "Zgodność z rozwiązaniem analitycznym (P = odcinek x1 = x2 ∈ [0, 1], "
                      "front √F1 + √F2 = √2), punkty niezdominowane:", ""]
        for label, result in results.items():
            F, X, _, mask = arrays(result)
            front_err = np.max(np.abs(np.sqrt(F[mask, 0]) + np.sqrt(F[mask, 1]) - math.sqrt(2)))
            set_err = np.max(np.abs(X[mask, 0] - X[mask, 1]))
            lines.append(f"- {label}: max |√F1 + √F2 − √2| = {front_err:.2e}, max |x1 − x2| = {set_err:.2e}")
    if exp_id == "a_fonseca":
        lines += ["", "Rozpiętość F1 punktów niezdominowanych (pokrycie frontu):", ""]
        for label, result in results.items():
            F, _, _, mask = arrays(result)
            values = np.sort(np.unique(np.round(F[mask, 0], 4)))
            gap = np.max(np.diff(values)) if len(values) > 1 else float("nan")
            lines.append(f"- {label}: {len(values)} różnych wartości F1, największa luka ΔF1 = {gap:.3f}")
    return lines + [""]


def main():
    report = []
    for exp_id, title, n, preset, domain, step in EXPERIMENTS:
        problem, table, ideal, nadir, reference, results = solve(n, preset, domain, step)
        (plot_2d if n == 2 else plot_3d)(exp_id, title, problem, ideal, nadir, reference, results)
        report += summary(exp_id, title, problem, table, ideal, nadir, reference, results)
        print("done", exp_id, flush=True)
    (ROOT / "results.md").write_text("\n".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
