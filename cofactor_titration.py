"""How the driving force and the redox pools respond to ethanol titer, one figure per
pathway.

For each flux mode this re-solves the MDF problem across a range of ethanol
concentrations and plots four quantities against titer: the MDF itself, and the
NADH/NAD+, NADPH/NADP+ and Fd(red)/Fd(ox) ratios at the optimum. Together they show
where a pathway spends its driving force as product accumulates, and which pool it has
to re-poise to keep going.

    python cofactor_titration.py                # every pathway
    python cofactor_titration.py --modes 1 2    # just these

Figures are written to ``results/figures/titration/`` and the underlying curves to
``results/tables/cofactor_titration.csv``.

Not every pathway contains every pool -- pathways without a ferredoxin-linked or
NADP-linked step leave the corresponding panel empty rather than dropping it, so the
figures stay comparable across modes.
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

# One panel per quantity, so no series share an axis and colour is not carrying identity
# between them; these are the reference palette's slots taken in fixed order.
PANELS = [
    ("mdf", "MDF (kJ/mol)", "#2a78d6", False),
    ("nadh/nad", "NADH/NAD$^+$", "#eb6834", True),
    ("nadph/nadp", "NADPH/NADP$^+$", "#1baf7a", True),
    ("fdred/fdox", "Fd$_{red}$/Fd$_{ox}$", "#4a3aa7", True),
]

# The concentration bounds hold every redox ratio between 1:100 and 100:1.
RATIO_MIN, RATIO_MAX = 0.01, 100.0


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=8, length=3)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def plot_mode(curve: pd.DataFrame, mode_id: int, name: str, max_gL: float):
    fig, axes = plt.subplots(1, 4, figsize=(11.0, 2.9))

    for ax, (column, label, color, log_scale) in zip(axes, PANELS):
        present = column in curve.columns and curve[column].notna().any()
        if present:
            ax.plot(curve.etoh_gL, curve[column], color=color, lw=2, zorder=3)
        else:
            ax.text(0.5, 0.5, "not used by\nthis pathway", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8, color=INK_MUTED,
                    linespacing=1.4)

        if log_scale:
            ax.set_yscale("log")
            ax.set_ylim(RATIO_MIN / 2, RATIO_MAX * 2)
            # the concentration bounds cap every ratio at 1:100 and 100:1
            for bound in (RATIO_MIN, RATIO_MAX):
                ax.axhline(bound, color=INK_MUTED, lw=0.8, ls=(0, (1, 3)), zorder=1)
            ax.axhline(1.0, color=GRID, lw=1.0, zorder=1)
        else:
            ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)

        ax.set_xlim(0, max_gL)
        ax.set_xlabel("ethanol (g/L)", fontsize=9)
        ax.set_ylabel(label, fontsize=9)
        style(ax)

    fig.suptitle(f"M{mode_id:02d}   {name}", fontsize=10, color=INK, x=0.005, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=Path, default=mdf.DEFAULT_MODEL)
    ap.add_argument("--modes", type=int, nargs="+", default=None)
    ap.add_argument("--max-gL", type=float, default=100.0)
    ap.add_argument("--points", type=int, default=40)
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    args = ap.parse_args(argv)

    modes = mdf.list_modes(args.model)
    mode_ids = args.modes if args.modes is not None else list(modes.index)

    print("loading eQuilibrator component-contribution data ...", flush=True)
    cc = mdf.ComponentContribution()

    # log-spaced, since every quantity here responds to the log of the titer
    titers = np.logspace(np.log10(0.1), np.log10(args.max_gL), args.points)

    figure_dir = args.outdir / "figures" / "titration"
    all_curves = []

    for mode_id in mode_ids:
        name = modes.loc[mode_id, "name"]
        curve = mdf.titrate_product(args.model, mode_id, titers, cc,
                                    out_dir=args.outdir / "sbtab")
        all_curves.append(curve.assign(mode_id=mode_id, name=name))

        fig = plot_mode(curve, mode_id, name, args.max_gL)
        stem = figure_dir / f"M{mode_id:02d}_{name.replace(' ', '_')}_titration"
        mdf.save_figure(fig, stem)
        plt.close(fig)

        span = f"{curve.mdf.iloc[0]:5.2f} -> {curve.mdf.iloc[-1]:5.2f}"
        ratio = curve["nadh/nad"]
        print(f"  M{mode_id:02d} {name:<40s} MDF {span} | "
              f"NADH/NAD+ {ratio.iloc[0]:.3f} -> {ratio.iloc[-1]:.3f}", flush=True)

    table = pd.concat(all_curves)
    table_path = args.outdir / "tables" / "cofactor_titration.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(table_path, index=False)
    print(f"\nwrote {len(mode_ids)} figures to {figure_dir}")
    print(f"wrote {table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
