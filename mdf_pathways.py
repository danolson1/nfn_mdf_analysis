"""
Max-min Driving Force (MDF) analysis of alternative ethanol production pathways.

The pathway model lives in an Excel workbook (``data/ethanol_pathway_model.xlsx``) with
two sheets:

``Reaction``
    One row per reaction, giving the reaction ID and formula, plus a block of flux
    columns -- one column per pathway ("elementary flux mode"). A blank cell means the
    reaction carries no flux in that mode and is dropped from the model.
``Compound``
    One row per metabolite, giving the KEGG identifier and the lower/upper concentration
    bounds (mM) used as constraints in the MDF linear program.

The thermodynamic conditions that apply to every pathway (pH, ionic strength, pMg,
temperature) are kept in ``data/thermodynamic_config.tsv`` rather than in the workbook,
so that the workbook never has to be rewritten programmatically. (The flux block mirrors
the reaction IDs with Excel formulas, and rewriting the file with openpyxl would discard
their cached values.)

For each mode the workbook is converted to an SBtab document, which is what
``equilibrator_pathway.ThermodynamicModel`` consumes. Writing the intermediate SBtab
file out makes each analysis independently reproducible and human-readable.

Reference: Noor et al. (2014) PLoS Comput Biol 10:e1003483, doi:10.1371/journal.pcbi.1003483
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from equilibrator_api import ComponentContribution
from equilibrator_pathway import ThermodynamicModel
from sbtab import SBtab, validatorSBtab

# Embed real TrueType fonts in vector output. Matplotlib's default (type 3) writes each
# glyph as PDF drawing operations, which leaves the text unselectable and unsearchable
# and is rejected by several publishers (Elsevier, IEEE, ACS). Set on import so the CLI
# and the notebook both produce submission-ready files.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

logger = logging.getLogger(__name__)

DEFAULT_MODEL = Path(__file__).parent / "data" / "ethanol_pathway_model.xlsx"
DEFAULT_CONFIG = Path(__file__).parent / "data" / "thermodynamic_config.tsv"

SBTAB_VERSION = "1.0"

# Rows 0-4 of the Reaction sheet are header rows: the group name ('info'/'rxn'/'stoich'/
# 'flux'), then the Mode ID / ATP yield / MDF annotation rows, then the column names.
_REACTION_HEADER_ROWS = [0, 1, 2, 3, 4]
_ANNOTATION_LEVELS = [0, 1, 2]  # dropped after selecting a column group


# --------------------------------------------------------------------------------------
# Reading the Excel model
# --------------------------------------------------------------------------------------

def read_config(config_path: Path | str = DEFAULT_CONFIG) -> pd.DataFrame:
    """Return the thermodynamic conditions as an Option/Value/Comment table."""
    config = pd.read_csv(config_path, sep="\t", dtype=str)
    return config.fillna("")


def _read_reaction_sheet(xlsx_path: Path | str) -> pd.DataFrame:
    return pd.read_excel(xlsx_path, sheet_name="Reaction", header=_REACTION_HEADER_ROWS)


def list_modes(xlsx_path: Path | str = DEFAULT_MODEL) -> pd.DataFrame:
    """List the pathways (flux modes) defined in the workbook.

    Returns a frame indexed by ``mode_id`` with the pathway ``name``, the ATP yield per
    glucose recorded in the workbook, and the number of reactions carrying flux.
    """
    flux_block = _read_reaction_sheet(xlsx_path).loc[:, "flux"]

    records = []
    for column in flux_block.columns:
        mode_id, atp_yield, _mdf, name = column
        if mode_id == "Mode ID":  # this is the reaction-ID column, not a pathway
            continue
        records.append(
            {
                "mode_id": int(mode_id),
                "name": name,
                "atp_per_glucose": atp_yield,
                "n_reactions": int(flux_block[column].notna().sum()),
            }
        )

    return pd.DataFrame.from_records(records).set_index("mode_id").sort_index()


def _read_tables(xlsx_path: Path | str, mode_id: int) -> dict:
    """Build the four SBtab data tables for one flux mode."""
    sheet = _read_reaction_sheet(xlsx_path)

    # --- reactions: ID + formula, for every reaction in the workbook ---
    reactions = sheet.loc[:, "rxn"].copy()
    reactions.columns = reactions.columns.droplevel(_ANNOTATION_LEVELS)
    reactions = reactions.rename(columns={"ReactionID": "ID"})
    reactions = reactions.dropna(subset="ID")

    # --- fluxes: the single flux column for this mode, aligned to the reaction rows ---
    # The flux block repeats the reaction IDs, but only as Excel formulas pointing back
    # at the rxn block, so take the IDs from the rxn block instead and align by row.
    flux_block = sheet.loc[:, "flux"]
    mode_columns = [c for c in flux_block.columns if c[0] == mode_id]
    if len(mode_columns) != 1:
        raise KeyError(f"mode {mode_id} is not defined exactly once in {xlsx_path}")
    mode_column = mode_columns[0]
    mode_name = mode_column[-1]

    fluxes = pd.DataFrame(
        {
            "QuantityType": "rate of reaction",
            "Reaction": reactions.ID,
            "Value": flux_block.loc[reactions.index, mode_column],
        }
    )

    # --- drop reactions that carry no flux in this mode ---
    active = fluxes.Value.notna()
    fluxes = fluxes.loc[active]
    reactions = reactions.loc[active[active].index]
    # SBtab serialisation requires strings
    fluxes["Value"] = fluxes["Value"].astype("str")
    if len(fluxes) != len(reactions) or fluxes.empty:
        raise ValueError(
            f"mode {mode_id}: {len(reactions)} reactions but {len(fluxes)} fluxes"
        )

    # --- compounds and concentration bounds ---
    compound_sheet = pd.read_excel(xlsx_path, sheet_name="Compound", header=1)
    compounds = compound_sheet.loc[:, ["CompoundID", "Name", "Identifiers"]].copy()
    compounds = compounds.rename(columns={"CompoundID": "ID"})

    bounds = compound_sheet.loc[:, ["CompoundID", "Min", "Max"]].copy()
    bounds = bounds.rename(columns={"CompoundID": "Compound"})
    bounds.insert(0, "QuantityType", "concentration")
    bounds["Min"] = bounds["Min"].astype("str")
    bounds["Max"] = bounds["Max"].astype("str")

    return {
        "name": mode_name,
        "reactions": reactions,
        "fluxes": fluxes,
        "compounds": compounds,
        "conc_constraints": bounds,
    }


# --------------------------------------------------------------------------------------
# Writing SBtab
# --------------------------------------------------------------------------------------

def _bang(df: pd.DataFrame) -> pd.DataFrame:
    """SBtab column headers are prefixed with an exclamation point."""
    out = df.copy()
    out.columns = [f"!{c}" for c in out.columns]
    return out


def build_sbtab(
    xlsx_path: Path | str,
    mode_id: int,
    out_dir: Path | str | None = None,
    config_path: Path | str = DEFAULT_CONFIG,
    validate: bool = False,
):
    """Convert one flux mode of the workbook into an SBtab document.

    Returns ``(document, path)``; ``path`` is ``None`` when ``out_dir`` is ``None``.
    """
    tables = _read_tables(xlsx_path, mode_id)
    mode_name = tables["name"]

    doc = SBtab.SBtabDocument()
    doc.set_name(mode_name)
    doc.set_filename(f"M{mode_id:02d}_{mode_name.replace(' ', '_')}.tsv")
    doc.set_version(SBTAB_VERSION)

    specs = [
        (_bang(read_config(config_path)), "Configuration", "Config", "Configuration", None),
        (_bang(tables["reactions"]), "Reaction", "Reaction", "Reaction", None),
        (_bang(tables["compounds"]), "Compound", "Compound", "Compound", None),
        (_bang(tables["fluxes"]), "Flux", "Quantity", "Flux", "dimensionless"),
        (
            _bang(tables["conc_constraints"]),
            "ConcentrationConstraint",
            "Quantity",
            "ConcentrationConstraint",
            "mM",
        ),
    ]

    for frame, table_id, table_type, table_name, unit in specs:
        kwargs = {} if unit is None else {"unit": unit}
        doc.add_sbtab(
            SBtab.SBtabTable.from_data_frame(
                frame,
                table_id=table_id,
                table_type=table_type,
                table_name=table_name,
                sbtab_version=doc.version,
                **kwargs,
            )
        )

    if validate:
        for message in validatorSBtab.ValidateDocument(doc).validate_document():
            logger.warning("SBtab validation: %s", message)

    path = None
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / doc.filename
        path.write_text(doc.to_str(), encoding="utf-8")

    return doc, path


# --------------------------------------------------------------------------------------
# Running MDF
# --------------------------------------------------------------------------------------

@dataclass
class MDFResult:
    """Everything needed to report or plot one pathway."""

    mode_id: int
    name: str
    mdf_kj_per_mol: float
    net_reaction: str
    reaction_df: pd.DataFrame
    compound_df: pd.DataFrame
    solution: object  # equilibrator_pathway PathwayMdfSolution
    sbtab_path: Path | None = None

    @property
    def bottleneck_reactions(self) -> list:
        """Reactions operating at the MDF, i.e. the thermodynamic bottleneck.

        These are the steps highlighted in red on the driving-force plot: at the
        optimum their driving force (-dG') equals the MDF to within rounding, and
        they are the steps whose dG'0 uncertainty or concentration bounds set the
        MDF for the whole pathway.
        """
        driving_force = self.reaction_df["optimized_dg_prime"].apply(
            lambda q: -q.m_as("kJ/mol")
        )
        at_optimum = (driving_force - self.mdf_kj_per_mol).abs() < 1e-3
        return list(self.reaction_df.reaction_id[at_optimum])

    def ratios(self) -> dict:
        """Cofactor ratios and key metabolite concentrations at the MDF optimum."""
        conc = self.compound_df.set_index("compound_id")["concentration_in_mM"]
        out = {}
        for numerator, denominator in [
            ("nadh", "nad"),
            ("nadph", "nadp"),
            ("fdred", "fdox"),
            ("atp", "adp"),
        ]:
            if numerator in conc.index and denominator in conc.index:
                out[f"{numerator}/{denominator}"] = conc[numerator] / conc[denominator]
        for metabolite in ("acald", "accoa"):
            if metabolite in conc.index:
                out[f"[{metabolite}] mM"] = conc[metabolite]
        return out


def _tidy_compound_df(compound_df: pd.DataFrame) -> pd.DataFrame:
    """Add mM columns alongside eQuilibrator's molar concentration columns."""
    df = compound_df.copy()
    for molar_col, mM_col in [
        ("concentration_in_molar", "concentration_in_mM"),
        ("lower_bound_in_molar", "lower_bound_in_mM"),
        ("upper_bound_in_molar", "upper_bound_in_mM"),
    ]:
        if molar_col in df.columns:
            df[mM_col] = df[molar_col] * 1000.0
    return df


def run_mdf(
    xlsx_path: Path | str = DEFAULT_MODEL,
    mode_id: int = 5,
    out_dir: Path | str = "results/sbtab",
    config_path: Path | str = DEFAULT_CONFIG,
    comp_contrib: ComponentContribution | None = None,
) -> MDFResult:
    """Build the SBtab model for one flux mode and solve the MDF problem.

    ``ComponentContribution`` is expensive to construct; pass one in when looping over
    several modes.
    """
    cc = comp_contrib or ComponentContribution()
    doc, sbtab_path = build_sbtab(
        xlsx_path, mode_id, out_dir=out_dir, config_path=config_path
    )

    model = ThermodynamicModel.from_sbtab(str(sbtab_path), comp_contrib=cc)
    solution = model.mdf_analysis()

    return MDFResult(
        mode_id=mode_id,
        name=doc.name,
        mdf_kj_per_mol=float(solution.score),
        net_reaction=model.net_reaction_formula,
        reaction_df=solution.reaction_df,
        compound_df=_tidy_compound_df(solution.compound_df),
        solution=solution,
        sbtab_path=sbtab_path,
    )


# --------------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------------

def plot_driving_forces(
    result: MDFResult,
    ax: plt.Axes | None = None,
    ylim=(-145, 5),
    figsize=(3, 3),
) -> plt.Figure:
    """Cumulative driving-force plot, in the style used for the published figures.

    The same y-range is used for every pathway so that panels are directly comparable.
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)
    fig = ax.get_figure()

    result.solution.plot_driving_forces(ax=ax)
    ax.set_ylim(ylim)
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=90, ha="right")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.3),
        fancybox=False,
        shadow=False,
        ncol=1,
    )
    return fig


def save_figure(fig: plt.Figure, stem: Path | str, formats=("png", "pdf"), dpi: int = 200):
    """Write a figure to one file per format, with the axes and legend box included."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = stem.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    return paths
