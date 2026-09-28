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

WT_MODE = 5   # Cth ethanologen, wild-type NADH-linked AdhE
NFN_MODE = 4  # NfnAB with NADPH-linked AdhE

# Categorical slots 1-3 of the reference palette, which validate on all pairs in both
# modes. Colour encodes the source, consistently across both panels.
C_MODEL_WT = "#2a78d6"   # slot 1, blue   -- MDF model, native pathway
C_MEASURED = "#eb6834"   # slot 2, orange -- Tian et al. 2017 measurement
C_MODEL_NFN = "#1baf7a"  # slot 3, aqua   -- MDF model, NfnAB pathway

INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d8d7d2"

# From Tian et al. (2017) Fig. 4: specific activity of purified Gapdh versus NADH/NAD+.
# The C. thermocellum enzyme loses more than half its activity at 0.2 and essentially
# all of it at 1.0, where the T. saccharolyticum enzyme retains ~30%.
GAPDH_HALF = 0.2
GAPDH_OFF = 1.0


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
    nfn = mdf_vs_titer(NFN_MODE, titers, cc, args.model, sbtab_dir)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.0))

    # --- Panel A: predicted vs measured NADH/NAD+ -------------------------------------
    ax1.axhspan(GAPDH_OFF, 2.4, color=GRID, alpha=0.28, lw=0, zorder=0)
    ax1.axhline(GAPDH_OFF, color=INK_MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax1.axhline(GAPDH_HALF, color=INK_MUTED, lw=1.0, ls=(0, (1, 3)), zorder=1)
    ax1.text(59, GAPDH_OFF + 0.05, "Cth Gapdh inactive", ha="right", va="bottom",
             fontsize=8, color=INK_MUTED)
    ax1.text(59, GAPDH_HALF + 0.04, "Cth Gapdh 50% activity", ha="right", va="bottom",
             fontsize=8, color=INK_MUTED)

    ax1.plot(wt.etoh_gL, wt.nadh_nad, color=C_MODEL_WT, lw=2,
             label="MDF model, native pathway", zorder=3)
    ax1.plot(measured.etoh_gL, measured.nadh_nad, color=C_MEASURED, lw=2,
             marker="o", ms=5, mew=1.2, mec="white",
             label="Measured (Tian et al. 2017)", zorder=4)

    ax1.set_xlabel("ethanol concentration (g/L)")
    ax1.set_ylabel(r"NADH/NAD$^+$ ratio")
    ax1.set_xlim(0, 60)
    ax1.set_ylim(0, 2.4)
    style(ax1)
    ax1.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=INK)
    ax1.set_title("Cofactor ratio rises with product titer",
                  fontsize=10, color=INK, loc="left", pad=8)

    # --- Panel B: MDF versus titer ----------------------------------------------------
    ax2.axhline(0, color=INK_MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax2.plot(nfn.etoh_gL, nfn.mdf, color=C_MODEL_NFN, lw=2, zorder=3)
    ax2.plot(wt.etoh_gL, wt.mdf, color=C_MODEL_WT, lw=2, zorder=3)

    ax2.text(30, nfn.mdf.iloc[-1] - 0.35, "NfnAB + engineered AdhE", ha="center",
             va="top", fontsize=9, color=C_MODEL_NFN)
    ax2.text(30, np.interp(30, wt.etoh_gL, wt.mdf) - 0.35, "native pathway",
             ha="center", va="top", fontsize=9, color=C_MODEL_WT)

    ax2.set_xlabel("ethanol concentration (g/L)")
    ax2.set_ylabel("MDF (kJ/mol)")
    ax2.set_xlim(0, 60)
    ax2.set_ylim(-0.4, 7.4)
    style(ax2)
    ax2.set_title("Driving force lost to product accumulation",
                  fontsize=10, color=INK, loc="left", pad=8)

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
