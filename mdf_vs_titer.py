"""Plot MDF against ethanol titer for two or more pathways.

The default pair is the two *T. saccharolyticum* routes that differ only in the cofactor
used by the alcohol dehydrogenase:

* mode 1, the NADH-linked ethanologen: hydrogen cycling supplies NADH, and both ALDH and
  ADH are NADH-linked.
* mode 2, the NADPH-linked ethanologen: NfnAB supplies NADPH to an NADPH-linked ADH,
  while ALDH stays NADH-linked.

Engineered *T. saccharolyticum* strains that acquired NADPH-linked ADH activity reach
higher ethanol titers, and this is the thermodynamic reason why: sharing the NADH pool
between GAPDH and the terminal reductions costs driving force that grows with titer.

    python mdf_vs_titer.py                  # modes 1 and 2
    python mdf_vs_titer.py --modes 5 4      # native Cth vs NfnAB + engineered AdhE
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import mdf_pathways as mdf  # noqa: E402

# Categorical slots of the reference palette, assigned in fixed order.
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d8d7d2"

# Short labels for the pathways this script is usually pointed at.
LABELS = {
    1: "NADH-linked ADH",
    2: "NADPH-linked ADH",
    4: "NfnAB + engineered AdhE",
    5: "native C. thermocellum",
}


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=3)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=Path, default=mdf.DEFAULT_MODEL)
    ap.add_argument("--modes", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--max-gL", type=float, default=100.0)
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    ap.add_argument("--name", default=None, help="output file stem")
    args = ap.parse_args(argv)

    modes = mdf.list_modes(args.model)
    print("loading eQuilibrator component-contribution data ...", flush=True)
    cc = mdf.ComponentContribution()

    titers = np.concatenate([np.linspace(0.1, 2, 12),
                             np.linspace(2.5, args.max_gL, 80)])

    curves = {}
    for mode_id in args.modes:
        curves[mode_id] = mdf.titrate_product(
            args.model, mode_id, titers, cc, out_dir=args.outdir / "sbtab"
        )
        print(f"  M{mode_id:02d} {modes.loc[mode_id, 'name']}", flush=True)

    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    ax.axhline(0, color=INK_MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)

    for i, mode_id in enumerate(args.modes):
        curve = curves[mode_id]
        label = LABELS.get(mode_id, modes.loc[mode_id, "name"])
        ax.plot(curve.etoh_gL, curve.mdf, color=SERIES_COLORS[i % len(SERIES_COLORS)],
                lw=2, label=f"M{mode_id:02d}  {label}", zorder=3)

    ax.set_xlabel("ethanol concentration (g/L)")
    ax.set_ylabel("MDF (kJ/mol)")
    ax.set_xlim(0, args.max_gL)
    style(ax)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="upper right")

    fig.tight_layout()
    stem = args.name or "mdf_vs_titer_" + "_".join(f"M{m:02d}" for m in args.modes)
    paths = mdf.save_figure(fig, args.outdir / "figures" / stem)
    plt.close(fig)

    table = pd.concat(
        [c.assign(mode_id=m, name=modes.loc[m, "name"]) for m, c in curves.items()]
    )
    table_path = args.outdir / "tables" / f"{stem}.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(table_path, index=False)

    # where each pathway runs out of driving force, and the gap between the pair
    for mode_id, curve in curves.items():
        below = curve[curve.mdf <= 0]
        limit = "beyond the range tested" if below.empty else f"{below.etoh_gL.iloc[0]:.0f} g/L"
        print(f"  M{mode_id:02d} MDF reaches zero at: {limit}")
    for p in paths:
        print("wrote", p)
    print("wrote", table_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
