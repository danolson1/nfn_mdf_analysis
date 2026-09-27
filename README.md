# Thermodynamic (MDF) analysis of ethanol production pathways

Max-min Driving Force (MDF) analysis of alternative routes from glucose to ethanol in
the thermophilic anaerobes *Clostridium thermocellum* and *Thermoanaerobacterium
saccharolyticum*, evaluated at a high ethanol titer.

This repository contains the code, the pathway model, and the generated figures and
tables for the MDF analysis reported in [CITATION TO BE ADDED].

## What the analysis does

A pathway can only carry flux in the forward direction if every one of its reactions has
a negative ΔrG′. At a high product titer this becomes restrictive: mass action pushes the
terminal steps toward equilibrium. The MDF (Noor et al. 2014) is the largest value of
−ΔrG′ that can be achieved simultaneously by *every* reaction in a pathway, optimizing
over all metabolite concentrations allowed by a set of bounds. It is a single number
summarizing how much thermodynamic room a pathway has, and it identifies which reactions
are the bottleneck.

Here it is used to ask which ethanol pathway — and in particular which route for
delivering electrons from reduced ferredoxin to the alcohol dehydrogenase — retains the
most driving force at 2 M (≈9% w/v) ethanol.

## Repository layout

```
data/
  ethanol_pathway_model.xlsx    pathway model: reactions, flux modes, concentration bounds
  thermodynamic_config.tsv      pH, ionic strength, pMg, temperature
mdf_pathways.py                 library: Excel -> SBtab -> MDF -> figure
run_mdf_analysis.py             command-line driver
notebooks/
  mdf_analysis.ipynb            worked example / interactive exploration
results/
  mdf_summary.csv               MDF, ATP yield, bottleneck reactions, cofactor ratios
  sbtab/                        the SBtab model actually solved, one file per pathway
  figures/                      driving-force plots (PNG + PDF)
  tables/                       per-pathway reaction and compound detail
```

## The model

`data/ethanol_pathway_model.xlsx` holds the whole model in two sheets.

**`Reaction`** lists 30 reactions with their formulas: glycolysis (10), the malate shunt
(3), three pyruvate-oxidizing options (PFOR, PDC, PDH), seven ferredoxin/NAD(P)
electron-transfer options (RNF, bifurcating hydrogenase, ferredoxin-only hydrogenase,
NfnAB, NADH- and NADPH-linked ferredoxin:NAD(P) oxidoreductases, and a pyridine
nucleotide transhydrogenase), and seven aldehyde/alcohol dehydrogenase options differing
in cofactor specificity and in whether acetaldehyde is channeled. To the right is a block of flux columns, one
per pathway. Each column is an elementary flux mode: a stoichiometrically balanced route
from 1 glucose to 2 ethanol. A blank cell means the reaction carries no flux in that
pathway and is dropped before solving. The header rows of each column record the mode ID,
the ATP yield per glucose, and the pathway name.

**`Compound`** lists the 30 metabolites with KEGG identifiers and the lower/upper
concentration bounds used in the linear program.

Fifteen pathways are defined:

| Mode | Pathway | ATP/glucose |
|-----:|---------|------------:|
|  1 | Tsac NADH ethanologen | 2 |
|  2 | Tsac NADH+NADPH ethanologen | 2 |
|  3 | Tsac NADH+NADPH ethanologen fnorp | 2 |
|  4 | NFN-only pathway with engineered AdhE | 2 |
|  5 | Cth ethanologen WT adhE | 3 |
|  6 | Cth malate shunt adhp | 3 |
|  7 | PDC ethanol pathway | 2 |
|  8 | PDC Xhyd ethanol pathway | 2 |
|  9 | Tsac fd adh | 2 |
| 10 | Tsac NADH substrate-channeling | 2 |
| 11 | Tsac NADPH substrate-channeling | 2 |
| 12 | PDH ethanol pathway | 2 |
| 13 | PFOR FNOR ethanol pathway | 2 |
| 14 | PFOR FNOR Xhyd pathway | 2 |
| 15 | Cth malate shunt FNORp adhp | 2 |

### Hydrogen cycling

Modes **1, 2 and 10** use **hydrogen cycling** to move electrons from ferredoxin to NAD,
as *T. saccharolyticum* does: the ferredoxin-only hydrogenase (`hyd`) evolves H₂, and the
electron-bifurcating hydrogenase (`bif-hyd`) runs in the reverse, H₂-consuming direction
to reduce NAD⁺, giving a net Fd:NAD transhydrogenation with H₂ as a freely diffusing
intermediate. The reverse direction is entered in the workbook as a negative flux on
`bif-hyd`. For mode 1, for example, PFOR produces 4 reduced ferredoxin, `hyd` (flux 4)
consumes 8, and `bif-hyd` (flux −2) regenerates 4 while reducing 2 NAD⁺ — balanced in
both ferredoxin and H₂.

Because the intermediate H₂ pool is explicit, the assumed H₂ concentration (0.01 mM)
directly controls how much driving force this route costs; above roughly 2 mM H₂ it
becomes thermodynamically equivalent to a direct NADH-linked ferredoxin:NAD
oxidoreductase.

Modes 9 and 11 are named for *T. saccharolyticum* but do **not** use hydrogen cycling —
mode 9 uses a ferredoxin-linked ADH directly, and mode 11 uses NfnAB with a channeled
NADPH-linked ALDH/ADH.

### Conditions

All pathways are evaluated under one set of conditions (`data/thermodynamic_config.tsv`):

| Option | Value |
|---|---|
| pH | 7.0 |
| ionic strength | 0.25 M |
| pMg | 3.0 |
| temperature | 328 K (55 °C) |
| dg_confidence | 0.95 |
| ln_conc_confidence | 0.95 |

ΔrG′° values are not treated as point estimates: the linear program lets them vary within
a 95% chi-squared confidence ellipsoid of the component-contribution covariance matrix,
and the reported MDF is the worst case over that ellipsoid. The `stdev_factor = 1.96` row
inherited from the original SBtab files is retained for provenance but is **inert** —
equilibrator-pathway 0.8.1 does not read it in the MDF code path. It encodes the same
95% interval, so the intent is unchanged, and setting `dg_confidence` explicitly gives
the same numbers as the defaults.

Concentrations are free to vary between 0.001 and 10 mM except where fixed to represent a
boundary condition or a buffered pool: glucose 10 mM, **ethanol 2000 mM**, phosphate
10 mM, CoA 1 mM, NAD⁺ / NADP⁺ / ADP / oxidized ferredoxin 0.1 mM, H₂ 0.01 mM, CO₂
0.001 mM, water 55 M. Holding the oxidized partner of each redox pair fixed while letting
the reduced partner range over 0.001–10 mM constrains each cofactor ratio to between
1:100 and 100:1; ATP is bounded at 0.1–10 mM, limiting ATP/ADP to 1:1–100:1.

The H₂ and CO₂ values are deliberately low (see the Notes column of the `Compound`
sheet): both are freely exchanged with the gas phase, and the chosen values illustrate
how gas removal affects pathways that produce or consume them.

## Reproducing the analysis

```bash
conda env create -f environment.yml
conda activate mdf-ethanol

python run_mdf_analysis.py --list          # show the pathways
python run_mdf_analysis.py                 # all 15 pathways
python run_mdf_analysis.py --modes 4 5     # just the two paper figures
```

Everything under `results/` is regenerated from scratch by that command; it is committed
so the repository can be read without running anything.

Pinned versions are in `requirements.txt`. The results were produced with
equilibrator-api 0.8.1 and Python 3.12.14 on Windows 11. `equilibrator-cache-data` ships
the component-contribution training data as a Python package, so no download is needed on
first use. The MDF linear program is solved through CVXPY; `equilibrator-pathway`
defaults to the CLARABEL solver, which can be overridden by adding a `solver` row to
`thermodynamic_config.tsv`.

### A note on editing the workbook

The flux block mirrors the reaction IDs from column C using Excel formulas. Saving the
workbook with a Python library such as `openpyxl` discards the cached values of those
formulas and leaves the column blank. Edit the workbook in Excel, not from Python. (The
parser reads reaction IDs from the `rxn` block and aligns fluxes by row, so it is
unaffected either way — but other tools may not be.)

## Results

| Mode | Pathway | MDF (kJ/mol) | Bottleneck |
|-----:|---------|-------------:|------------|
|  4 | NFN-only pathway with engineered AdhE | 7.26 | upper glycolysis |
| 11 | Tsac NADPH substrate-channeling | 7.24 | upper glycolysis |
|  9 | Tsac fd adh | 7.22 | upper glycolysis, aldh |
|  2 | Tsac NADH+NADPH ethanologen | 4.80 | upper glycolysis, bif-hyd, nfn, aldh, adhp |
|  3 | Tsac NADH+NADPH ethanologen fnorp | 4.79 | upper glycolysis, aldh, adhp |
| 15 | Cth malate shunt FNORp adhp | 4.08 | upper glycolysis, pepck |
|  7 | PDC ethanol pathway | 3.25 | upper glycolysis, adh |
|  6 | Cth malate shunt adhp | 2.97 | upper glycolysis, pepck, aldh, adhp |
| 10 | Tsac NADH substrate-channeling | 2.59 | upper glycolysis, aldh_adh |
|  8 | PDC Xhyd ethanol pathway | 2.46 | upper glycolysis, xhyd, adhp |
|  1 | Tsac NADH ethanologen | 2.20 | upper glycolysis, aldh, adh |
|  5 | Cth ethanologen WT adhE | 2.16 | upper glycolysis, aldh, adh |
| 13 | PFOR FNOR ethanol pathway | 2.16 | upper glycolysis, aldh, adh |
| 12 | PDH ethanol pathway | 2.12 | upper glycolysis, aldh, adh |
| 14 | PFOR FNOR Xhyd pathway | 1.88 | upper glycolysis, xhyd, aldh |

"Upper glycolysis" is the fba/tpi/gap segment, which is the thermodynamic bottleneck in
every pathway: it is limiting regardless of the ethanol route, and sets a ceiling of
about 7.3 kJ/mol under these conditions.

The pathways that reach that ceiling are those that avoid spending driving force on the
terminal aldehyde and alcohol dehydrogenase steps, either by coupling them to NADPH
generated by the electron-bifurcating transhydrogenase NfnAB (mode 4), by channeling the
acetaldehyde intermediate (mode 11), or by using a ferredoxin-linked ADH (mode 9). The
native *C. thermocellum* route with wild-type NADH-linked AdhE (mode 5, 2.16 kJ/mol) and
the native *T. saccharolyticum* route with hydrogen cycling (mode 1, 2.20 kJ/mol) are
both limited at their `aldh` and `adh` steps.

### Sensitivity to ethanol titer

Because the bottleneck differs, the two pathways respond very differently as ethanol
accumulates (last cell of the notebook):

| ethanol (mM) | M05 Cth WT adhE | M04 NFN + engineered AdhE |
|---:|---:|---:|
|   10 | 4.35 | 7.26 |
|  100 | 3.39 | 7.26 |
|  500 | 2.73 | 7.26 |
| 1000 | 2.44 | 7.26 |
| 1500 | 2.28 | 7.26 |
| 2000 | 2.16 | 7.26 |
| 2500 | 2.07 | 7.26 |

The MDF of the NFN pathway is entirely independent of the ethanol concentration: its
bottleneck is upper glycolysis, and the terminal steps retain enough driving force that
raising the titer does not make them limiting. In the wild-type pathway the ALDH and ADH
steps *are* the bottleneck, so every increase in ethanol is paid for directly out of the
pathway's driving force. Thermodynamically, the engineered pathway removes ethanol titer
as a constraint on the ethanol pathway itself.

## References

- Noor E, Bar-Even A, Flamholz A, Reznik E, Liebermeister W, Milo R (2014). Pathway
  thermodynamics highlights kinetic obstacles in central metabolism.
  *PLoS Comput Biol* 10:e1003483. doi:10.1371/journal.pcbi.1003483
- Beber ME, Gollub MG, Mozaffari D, Shebek KM, Flamholz AI, Milo R, Noor E (2022).
  eQuilibrator 3.0: a database solution for thermodynamic constant estimation.
  *Nucleic Acids Res* 50:D603–D609. doi:10.1093/nar/gkab1106
- Lubitz T, Hahn J, Bergmann FT, Noor E, Klipp E, Liebermeister W (2016). SBtab: a flexible
  table format for data exchange in systems biology. *Bioinformatics* 32:2559–2561.
  doi:10.1093/bioinformatics/btw179

## License

MIT — see `LICENSE`.
