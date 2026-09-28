"""One combined figure: every pathway as a row, sorted by MDF.

Each row is one flux mode and carries three panels:

1. the cumulative driving-force plot at the reference ethanol titer (2 M), the same
   plot ``run_mdf_analysis.py`` writes per pathway;
2. the MDF as a function of ethanol titer;
3. the NADH/NAD+, NADPH/NADP+ and Fd(red)/Fd(ox) ratios at the optimum, overlaid on one
   log axis, against ethanol titer.

Rows are sorted by the MDF at the reference titer, highest first, so the table of
pathways and the figure tell the same story in the same order.

    python combined_figure.py
    python combined_figure.py --modes 5 4 1 2     # a subset, still sorted by MDF
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import mdf_pathways as mdf  # noqa: E402

logging.getLogger("fontTools").setLevel(logging.ERROR)

INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d8d7d2"

# Panel 3 overlays three series on one axis, so these need to be distinguishable: the
# reference palette's first three categorical slots, which validate on all pairs.
RATIOS = [
    ("nadh/nad", r"NADH/NAD$^+$", "#2a78d6"),    # blue
    ("nadph/nadp", r"NADPH/NADP$^+$", "#1baf7a"),  # aqua
    ("fdred/fdox", r"Fd$_{red}$/Fd$_{ox}$", "#8c5a2b"),  # brown
]
MDF_COLOR = "#4a3aa7"

RATIO_MIN, RATIO_MAX = 0.01, 100.0

# Portrait, narrower than the page it sits on: at ~6 in wide the panels are close to
# square rather than wide and flat, and scaling the figure up to the text width of a
# document enlarges every label with it.
PAGE = (6.0, 11.0)


def style(ax, labelsize=6):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=labelsize, length=3)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=Path, default=mdf.DEFAULT_MODEL)
    ap.add_argument("--modes", type=int, nargs="+", default=None)
    ap.add_argument("--max-gL", type=float, default=100.0)
    ap.add_argument("--points", type=int, default=40)
    ap.add_argument("--ylim", type=float, nargs=2, default=(-145, 5))
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    ap.add_argument("--name", default="combined_pathway_overview",
                    help="output file stem (default: combined_pathway_overview)")
    ap.add_argument("--split", type=int, default=7,
                    help="rows on the first page; 0 for a single figure (default: 7)")
    ap.add_argument("--reuse-titration", action="store_true",
                    help="read the curves from results/tables/cofactor_titration.csv "
                         "instead of re-solving them (they are the same calculation)")
    args = ap.parse_args(argv)

    modes = mdf.list_modes(args.model)
    mode_ids = args.modes if args.modes is not None else list(modes.index)

    print("loading eQuilibrator component-contribution data ...", flush=True)
    cc = mdf.ComponentContribution()

    titers = np.logspace(np.log10(0.1), np.log10(args.max_gL), args.points)

    # solve everything first so the rows can be sorted by MDF before anything is drawn
    cached = None
    cache_path = args.outdir / "tables" / "cofactor_titration.csv"
    if args.reuse_titration and cache_path.exists():
        cached = pd.read_csv(cache_path)
        print(f"  reusing titration curves from {cache_path}", flush=True)

    results, curves = {}, {}
    for mode_id in mode_ids:
        results[mode_id] = mdf.run_mdf(args.model, mode_id,
                                       out_dir=args.outdir / "sbtab", comp_contrib=cc)
        if cached is not None and (cached.mode_id == mode_id).any():
            curves[mode_id] = cached[cached.mode_id == mode_id].sort_values("etoh_gL")
        else:
            curves[mode_id] = mdf.titrate_product(args.model, mode_id, titers, cc,
                                                  out_dir=args.outdir / "sbtab")
        print(f"  M{mode_id:02d} {modes.loc[mode_id, 'name']:<40s} "
              f"MDF = {results[mode_id].mdf_kj_per_mol:5.2f}", flush=True)

    order = sorted(mode_ids, key=lambda m: results[m].mdf_kj_per_mol, reverse=True)

    # one shared MDF scale across both pages, so rows compare directly
    mdf_hi = max(c.mdf.max() for c in curves.values())
    mdf_ylim = (-0.4, mdf_hi * 1.06)

    blocks = ([order[:args.split], order[args.split:]] if args.split
              else [order])
    blocks = [b for b in blocks if b]

    written = []
    for part, block in enumerate(blocks, start=1):
        fig = draw_block(block, results, curves, modes, args, mdf_ylim)
        stem = args.outdir / "figures" / args.name
        if len(blocks) > 1:
            stem = stem.with_name(f"{stem.name}_part{part}")
        written += mdf.save_figure(fig, stem, dpi=200)
        plt.close(fig)

    summary = pd.DataFrame(
        [{"rank": i + 1, "mode_id": m, "name": modes.loc[m, "name"],
          "mdf_kJ_per_mol": round(results[m].mdf_kj_per_mol, 2),
          "part": 1 + (i >= args.split if args.split else 0)}
         for i, m in enumerate(order)])
    print()
    print(summary.to_string(index=False))
    for p in written:
        print("wrote", p)
    return 0


def draw_block(order, results, curves, modes, args, mdf_ylim):
    """Draw one page: one row per pathway, three panels across."""
    n = len(order)
    fig, axes = plt.subplots(n, 3, figsize=PAGE,
                             gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})
    if n == 1:
        axes = np.array([axes])

    for row, mode_id in enumerate(order):
        result = results[mode_id]
        curve = curves[mode_id]
        ax1, ax2, ax3 = axes[row]

        # --- panel 1: cumulative driving force at the reference titer ----------------
        result.solution.plot_driving_forces(ax=ax1)
        ax1.set_ylim(args.ylim)
        ax1.set_xticks(ax1.get_xticks(), ax1.get_xticklabels(), rotation=90, ha="right")
        ax1.set_ylabel(r"cum. $\Delta_r G'$ (kJ/mol)", fontsize=6.5)
        ax1.set_xlabel("")
        # matplotlib keeps the left/center/right titles as separate artists, so
        # plot_driving_forces' centred "MDF = ..." survives unless it is cleared first
        ax1.set_title("")
        ax1.set_title(f"M{mode_id:02d}   {modes.loc[mode_id, 'name']}   "
                      f"(MDF = {result.mdf_kj_per_mol:.2f} kJ/mol)",
                      fontsize=6.5, color=INK, loc="left", pad=3)
        legend = ax1.get_legend()
        if legend is not None:
            legend.remove()
        if row == 0:
            ax1.legend(loc="lower left", fontsize=4.8, labelcolor=INK,
                       frameon=True, facecolor="white", framealpha=0.85,
                       edgecolor="none", handlelength=1.2, borderpad=0.25)
        style(ax1, labelsize=5.5)

        # --- panel 2: MDF against titer ----------------------------------------------
        ax2.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax2.plot(curve.etoh_gL, curve.mdf, color=MDF_COLOR, lw=2, zorder=3)
        ax2.set_xlim(0, args.max_gL)
        ax2.set_xticks([0, args.max_gL / 2, args.max_gL])
        ax2.set_ylim(mdf_ylim)                      # shared across every row and page
        ax2.set_xlabel("ethanol (g/L)", fontsize=6.5)
        ax2.set_ylabel("MDF (kJ/mol)", fontsize=6.5)
        style(ax2)

        # --- panel 3: the three redox ratios on one log axis --------------------------
        for column, label, color in RATIOS:
            if column in curve.columns and curve[column].notna().any():
                ax3.plot(curve.etoh_gL, curve[column], color=color, lw=2,
                         label=label, zorder=3)
        for bound in (RATIO_MIN, RATIO_MAX):
            ax3.axhline(bound, color=INK_MUTED, lw=0.8, ls=(0, (1, 3)), zorder=1)
        ax3.axhline(1.0, color=GRID, lw=1.0, zorder=1)
        ax3.set_yscale("log")
        ax3.set_ylim(RATIO_MIN / 2, RATIO_MAX * 2)
        ax3.set_xlim(0, args.max_gL)
        ax3.set_xticks([0, args.max_gL / 2, args.max_gL])
        ax3.set_xlabel("ethanol (g/L)", fontsize=6.5)
        ax3.set_ylabel("cofactor ratio", fontsize=6.5)
        style(ax3)
        # a legend on every row, because which pools a pathway uses changes row to row
        # above the axes rather than inside it: with only ~1.3 in of panel height there
        # is no interior band a curve is guaranteed not to reach
        ax3.legend(fontsize=5, labelcolor=INK, ncol=3, frameon=False,
                   loc="lower left", bbox_to_anchor=(-0.02, 1.0),
                   handlelength=1.0, columnspacing=0.7, borderpad=0.1,
                   handletextpad=0.35)

    fig.tight_layout(h_pad=1.0, w_pad=0.9, pad=0.6)
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
