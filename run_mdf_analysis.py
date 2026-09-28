"""Run MDF analysis for one or more pathways and write the SBtab models, figures
and a summary table.

Examples
--------
    python run_mdf_analysis.py                 # every pathway in the workbook
    python run_mdf_analysis.py --modes 4 5     # just the two paper figures
    python run_mdf_analysis.py --list          # show the pathways and exit
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")  # no interactive display needed
import matplotlib.pyplot as plt  # noqa: E402

import mdf_pathways as mdf  # noqa: E402

# Subsetting DejaVu for TrueType embedding makes fontTools complain about the font's
# head-table timestamps once per saved file. Harmless, and it drowns out real output.
logging.getLogger("fontTools").setLevel(logging.ERROR)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", type=Path, default=mdf.DEFAULT_MODEL,
                        help="Excel workbook defining the pathways")
    parser.add_argument("--modes", type=int, nargs="+", default=None,
                        help="mode IDs to analyse (default: all)")
    parser.add_argument("--outdir", type=Path, default=Path("results"),
                        help="output directory (default: results)")
    parser.add_argument("--ylim", type=float, nargs=2, default=(-145, 5),
                        help="y-axis range, kJ/mol (default: -145 5)")
    parser.add_argument("--list", action="store_true",
                        help="list the pathways in the workbook and exit")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    modes = mdf.list_modes(args.model)
    if args.list:
        print(modes.to_string())
        return 0

    mode_ids = args.modes if args.modes is not None else list(modes.index)

    # Building the ComponentContribution object takes several seconds, so do it once.
    print("loading eQuilibrator component-contribution data ...", flush=True)
    cc = mdf.ComponentContribution()

    sbtab_dir = args.outdir / "sbtab"
    figure_dir = args.outdir / "figures"

    rows = []
    for mode_id in mode_ids:
        result = mdf.run_mdf(args.model, mode_id, out_dir=sbtab_dir, comp_contrib=cc)

        fig = mdf.plot_driving_forces(result, ylim=tuple(args.ylim))
        stem = figure_dir / result.sbtab_path.stem
        mdf.save_figure(fig, stem)
        plt.close(fig)

        # per-pathway detail tables, for anyone who wants the underlying numbers
        detail_dir = args.outdir / "tables"
        detail_dir.mkdir(parents=True, exist_ok=True)
        result.compound_df.to_csv(detail_dir / f"{result.sbtab_path.stem}_compounds.csv",
                                  index=False)
        result.reaction_df.astype(str).to_csv(
            detail_dir / f"{result.sbtab_path.stem}_reactions.csv", index=False)

        row = {
            "mode_id": mode_id,
            "name": result.name,
            "mdf_kJ_per_mol": round(result.mdf_kj_per_mol, 2),
            "atp_per_glucose": modes.loc[mode_id, "atp_per_glucose"],
            "n_reactions": len(result.reaction_df),
            "bottleneck_reactions": ", ".join(result.bottleneck_reactions),
            "net_reaction": result.net_reaction,
        }
        row.update({k: round(v, 3) for k, v in result.ratios().items()})
        rows.append(row)

        print(f"  M{mode_id:02d} {result.name:<40s} MDF = "
              f"{result.mdf_kj_per_mol:6.2f} kJ/mol", flush=True)

    summary = pd.DataFrame(rows).sort_values("mdf_kJ_per_mol", ascending=False)
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "mdf_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\nwrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
