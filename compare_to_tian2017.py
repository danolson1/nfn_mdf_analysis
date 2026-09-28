"""Compare the MDF prediction against the measured metabolome of Tian et al. (2017).

MDF analysis of the native *C. thermocellum* pathway predicts that the NADH/NAD+ ratio
must rise as ethanol accumulates: GAPDH and the terminal ALDH/ADH steps draw on the same
cofactor pool and pull it in opposite directions, so the only way to keep the terminal
reductions moving against rising product is to let NADH/NAD+ climb -- at GAPDH's expense.

Tian et al. (2017) measured exactly that, by adding ethanol to a growing culture at
~5 g/L/h and quantifying intracellular metabolites by LC-MS. They independently found the
metabolite signature of a GAPDH bottleneck (accumulation upstream, depletion downstream)
and showed in vitro that the *C. thermocellum* enzyme is inhibited by high NADH/NAD+.

This script overlays the two.

    python compare_to_tian2017.py

Reference: Tian L, Perot SJ, Stevenson D, Jacobson T, Lanahan AA, Amador-Noguez D,
Olson DG, Lynd LR (2017) Metabolome analysis reveals a role for glyceraldehyde
3-phosphate dehydrogenase in the inhibition of C. thermocellum by ethanol.
Biotechnol Biofuels 10:276. doi:10.1186/s13068-017-0961-3
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from equilibrator_api import Q_  # noqa: E402
from equilibrator_pathway import ThermodynamicModel  # noqa: E402

import mdf_pathways as mdf  # noqa: E402

MW_ETOH = 46.069  # g/mol -> 1 mM = 0.046069 g/L

WT_MODE = 5  # Cth ethanologen, wild-type NADH-linked AdhE

# Categorical slots 1 and 2 of the reference palette. Colour encodes the source of the
# number: measured or modelled.
C_MEASURED = "#eb6834"  # slot 2, orange -- Tian et al. 2017 measurement
C_MODEL_WT = "#2a78d6"  # slot 1, blue   -- MDF model, native pathway

INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d8d7d2"


def mdf_vs_titer(mode_id, titers_gL, cc, model_path, sbtab_dir):
    """MDF and optimized NADH/NAD+ for one pathway across a range of ethanol titers."""
    _, sbtab_path = mdf.build_sbtab(model_path, mode_id, out_dir=sbtab_dir)
    rows = []
    for gL in titers_gL:
        mM = gL / (MW_ETOH / 1000.0)
        model = ThermodynamicModel.from_sbtab(str(sbtab_path), comp_contrib=cc)
        model.set_bounds("etoh", Q_(mM, "mM"), Q_(mM, "mM"))
        sol = model.mdf_analysis()
        conc = mdf._tidy_compound_df(sol.compound_df).set_index("compound_id")
        conc = conc["concentration_in_mM"]
        rows.append(
            {
                "etoh_gL": gL,
                "etoh_mM": mM,
                "mdf": float(sol.score),
                "nadh_nad": conc["nadh"] / conc["nad"],
            }
        )
    return pd.DataFrame(rows)


def style(ax):
    """Recessive axes: no top/right spines, muted grid behind the data."""
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
    ap.add_argument("--measured", type=Path,
                    default=Path(__file__).parent / "data" / "tian2017_nadh_nad.csv")
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    args = ap.parse_args(argv)

    measured = pd.read_csv(args.measured, comment="#")

    print("loading eQuilibrator component-contribution data ...", flush=True)
    cc = mdf.ComponentContribution()

    titers = np.concatenate([np.linspace(0.2, 2, 10), np.linspace(2.5, 60, 60)])
    sbtab_dir = args.outdir / "sbtab"
    wt = mdf_vs_titer(WT_MODE, titers, cc, args.model, sbtab_dir)

    # Two panels on a shared scale, so the shape of the two curves can be compared
    # directly. One series each, so the panel title names it and no legend is needed.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.6), sharey=True)

    # --- Panel A: measured ------------------------------------------------------------
    ax1.plot(measured.etoh_gL, measured.nadh_nad, color=C_MEASURED, lw=2,
             marker="o", ms=6, mew=1.2, mec="white", zorder=3)
    ax1.set_ylabel(r"NADH/NAD$^+$ ratio")
    ax1.set_title("A   Measured", fontsize=10, color=INK, loc="left", pad=8)
    ax1.text(0.97, 0.06, "Tian et al. 2017", transform=ax1.transAxes, ha="right",
             va="bottom", fontsize=9, color=C_MEASURED)

    # --- Panel B: predicted -----------------------------------------------------------
    ax2.plot(wt.etoh_gL, wt.nadh_nad, color=C_MODEL_WT, lw=2, zorder=3)
    ax2.set_title("B   Predicted", fontsize=10, color=INK, loc="left", pad=8)
    ax2.text(0.97, 0.06, "MDF analysis,\nnative pathway", transform=ax2.transAxes,
             ha="right", va="bottom", fontsize=9, color=C_MODEL_WT, linespacing=1.3)

    for ax in (ax1, ax2):
        ax.set_xlabel("ethanol concentration (g/L)")
        ax.set_xlim(0, 60)
        ax.set_ylim(0, 2.0)
        style(ax)

    fig.tight_layout()
    figure_dir = args.outdir / "figures"
    paths = mdf.save_figure(fig, figure_dir / "tian2017_comparison")
    plt.close(fig)

    table = args.outdir / "tables" / "tian2017_comparison.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    merged = measured.copy()
    merged["predicted_nadh_nad"] = np.interp(
        merged.etoh_gL, wt.etoh_gL, wt.nadh_nad
    )
    merged["predicted_mdf"] = np.interp(merged.etoh_gL, wt.etoh_gL, wt.mdf)
    merged.round(3).to_csv(table, index=False)

    print(merged.round(3).to_string(index=False))
    for p in paths:
        print("wrote", p)
    print("wrote", table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
