# Thermodynamic (MDF) analysis of ethanol production pathways

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23005004.svg)](https://doi.org/10.5281/zenodo.23005004)

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
  tian2017_nadh_nad.csv         measured NADH/NAD+ vs ethanol, Tian et al. 2017
mdf_pathways.py                 library: Excel -> SBtab -> MDF -> figure
run_mdf_analysis.py             command-line driver
compare_to_tian2017.py          overlay the prediction on measured metabolome data
mdf_vs_titer.py                 MDF versus ethanol titer, for any set of pathways
cofactor_titration.py           MDF and the redox ratios versus titer, one figure per mode
combined_figure.py              all pathways in one figure, one row each, sorted by MDF
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

### Reactions

The reaction IDs below are what appear on the horizontal axis of every driving-force
plot. The workbook is authoritative for the formulas; this table adds the enzyme names.

| ID | Enzyme | Gene | Formula |
|---|---|---|---|
| `glk` | glucokinase | | glc + atp ⇌ g6p + adp |
| `pgi` | phosphoglucose isomerase | | g6p ⇌ f6p |
| `pfk` | phosphofructokinase | | atp + f6p ⇌ adp + fbp |
| `fba` | fructose-bisphosphate aldolase | | fbp ⇌ dhap + g3p |
| `tpi` | triose-phosphate isomerase | | dhap ⇌ g3p |
| `gap` | glyceraldehyde-3-phosphate dehydrogenase | | pi + nad + g3p ⇌ nadh + bpg |
| `pgk` | phosphoglycerate kinase | | bpg + adp ⇌ 3pg + atp |
| `gpm` | phosphoglycerate mutase | | 3pg ⇌ 2pg |
| `eno` | enolase | | 2pg ⇌ pep + h2o |
| `pyk` | pyruvate kinase | | adp + pep ⇌ atp + pyr |
| `pepck` | PEP carboxykinase | | pep + adp + co2 ⇌ oaa + atp |
| `mdh` | malate dehydrogenase | | oaa + nadh ⇌ nad + mal |
| `mae` | malic enzyme (NADP-linked) | | mal + nadp ⇌ nadph + pyr + co2 |
| `pfor` | pyruvate:ferredoxin oxidoreductase | `pforA` | coa + pyr + 2 fdox ⇌ accoa + co2 + 2 fdred |
| `pdc` | pyruvate decarboxylase | | pyr ⇌ acald + co2 |
| `pdh` | pyruvate dehydrogenase | | pyr + coa + nad ⇌ accoa + nadh + co2 |
| `rnf` | Rnf ferredoxin:NAD⁺ oxidoreductase (ion-translocating) | `rnf` | 2 fdred + nad + 0.5 adp + 0.5 pi ⇌ 0.5 atp + 0.5 h2o + nadh + 2 fdox |
| `bif-hyd` | electron-bifurcating hydrogenase | `hydA` | 2 fdred + nadh ⇌ 2 fdox + 2 h2 + nad |
| `hyd` | ferredoxin-only hydrogenase | `hfsD` | 2 fdred ⇌ 2 fdox + h2 |
| `nfn` | NfnAB electron-bifurcating transhydrogenase | `nfnAB` | 2 fdred + nadh + 2 nadp ⇌ 2 nadph + nad + 2 fdox |
| `fnor` | ferredoxin:NAD⁺ oxidoreductase | | 2 fdred + nad ⇌ nadh + 2 fdox |
| `fnorp` | ferredoxin:NADP⁺ oxidoreductase | `nfnB` alone, `cac_0764` | 2 fdred + nadp ⇌ nadph + 2 fdox |
| `xhyd` | pyridine nucleotide transhydrogenase | | nadh + nadp ⇌ nadph + nad |
| `aldh` | acetaldehyde dehydrogenase, NADH-linked | `adhE` WT | accoa + nadh ⇌ coa + acald + nad |
| `aldhp` | acetaldehyde dehydrogenase, NADPH-linked | `adhE` PROSS | accoa + nadph ⇌ coa + acald + nadp |
| `adh` | alcohol dehydrogenase, NADH-linked | `adhE` WT | acald + nadh ⇌ nad + etoh |
| `adhp` | alcohol dehydrogenase, NADPH-linked | `adhE` mutants, `adhA` | acald + nadph ⇌ nadp + etoh |
| `adhfd` | alcohol dehydrogenase, ferredoxin-linked | | acald + 2 fdred ⇌ 2 fdox + etoh |
| `aldh_adh` | channeled ALDH/ADH, NADH-linked | `adhE` with substrate channeling | accoa + 2 nadh ⇌ coa + 2 nad + etoh |
| `aldhp_adhp` | channeled ALDH/ADH, NADPH-linked | `adhE` with substrate channeling | accoa + 2 nadph ⇌ coa + 2 nadp + etoh |

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

> **On the channeled modes.** Modes 10 and 11 assume acetaldehyde is transferred directly
> between the two domains of AdhE. Fabri et al. (2025) call that substrate channeling into
> question, so those two modes should be read as an upper bound on what channeling could
> buy rather than as a description of the enzyme as it is. The same paper notes that the
> route from reduced ferredoxin to NAD(P)⁺, and which cofactor ALDH and ADH use in vivo,
> both remain open questions — the *T. saccharolyticum* modes here are hypotheses about
> that, not settled fact.

### Conditions

All pathways are evaluated under one set of conditions (`data/thermodynamic_config.tsv`):

| Option | Value |
|---|---|
| pH | 7.0 |
| ionic strength | 0.25 M |
| pMg | 3.0 |
| temperature | 328 K (55 °C) |
| dg_confidence | 0.0 |
| ln_conc_confidence | 0.95 |

**ΔrG′° is used as a point estimate** (`dg_confidence = 0`). This matters. With a
non-zero value, equilibrator-pathway treats the ΔrG′° vector as a *decision variable*
constrained to that confidence ellipsoid of the component-contribution covariance matrix
and **maximizes** over it — so the result is a best-case bound, not a worst-case one, and
not the expected value. At the library default of 0.95 it inflates every MDF here by
0.7–5.7 kJ/mol and makes both infeasible malate-shunt pathways appear feasible.

`ln_conc_confidence` is left at its default of 0.95; with explicit min/max bounds read
from a model file it reproduces those bounds exactly and has no other effect.

The `stdev_factor = 1.96` row inherited from the original SBtab files is retained for
provenance but is **inert** — equilibrator-pathway 0.8.1 does not read it in the MDF code
path (it survives only as a TODO in the ECM code). Do not mistake it for uncertainty
handling.

Concentrations are free to vary between 0.001 and 10 mM except where fixed to represent a
boundary condition or a buffered pool: glucose 10 mM, **ethanol 2000 mM**, phosphate
10 mM, CoA 1 mM, NAD⁺ / NADP⁺ / ADP / oxidized ferredoxin 0.1 mM, H₂ 0.01 mM, water
55 M. Holding the oxidized partner of each redox pair fixed while letting
the reduced partner range over 0.001–10 mM constrains each cofactor ratio to between
1:100 and 100:1; ATP is bounded at 0.1–10 mM, limiting ATP/ADP to 1:1–100:1.

**CO₂ is free over 0.001–10 mM**, the same default range as most other metabolites,
rather than fixed. This matters for the malate shunt, whose first step (PEP
carboxykinase) consumes CO₂. Dash et al. (2019) used a low fixed value at which the
shunt is thermodynamically infeasible; but Olson et al. (2017) showed that a *ppdk*
deletion strain of *C. thermocellum* still produces ethanol, which requires flux through
the malate shunt, so the shunt must be feasible in vivo. Leaving CO₂ free lets each
pathway sit where its own thermodynamics allow: the optimizer pushes CO₂ down to ~0.002–
0.5 mM in the pathways that produce it through PFOR, PDC or PDH, and up to ~0.6–1.4 mM
in the malate-shunt pathways that consume it at PEP carboxykinase.

> **Caveat.** CO₂ is exchanged with the gas phase, so in a real fermentation every
> pathway would see the same dissolved concentration. Letting it float per pathway gives
> each the benefit of the doubt, which is the right choice when the question is whether a
> pathway is feasible at all, but it is more permissive than fixing a common value.

H₂ is fixed at 0.01 mM, low enough that hydrogen cycling costs driving force; above
roughly 2 mM it becomes thermodynamically equivalent to a direct ferredoxin:NAD
oxidoreductase.

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

Re-running in a fresh clone reproduces the committed `mdf_summary.csv`, SBtab models,
detail tables and PNG figures byte-for-byte. Only the PDFs differ, because matplotlib
embeds a creation timestamp; set `SOURCE_DATE_EPOCH` (which matplotlib honors) to make
those reproducible too:

```bash
SOURCE_DATE_EPOCH=1727395200 python run_mdf_analysis.py
```

Pinned versions are in `requirements.txt`. The results were produced with
equilibrator-api 0.8.1 and Python 3.12.14 on Windows 11. `equilibrator-cache-data` ships
the component-contribution training data as a Python package, so no download is needed on
first use. The MDF linear program is solved through CVXPY; `equilibrator-pathway`
defaults to the CLARABEL solver, which can be overridden by adding a `solver` row to
`thermodynamic_config.tsv`.

## Results

MDF values use the component-contribution point estimates (`dg_confidence = 0`).

| Mode | Pathway | ATP/glc | MDF (kJ/mol) | Bottleneck | NADH/NAD⁺ | NADPH/NADP⁺ |
|-----:|---------|--------:|-------------:|------------|----------:|------------:|
|  4 | NFN-only pathway with engineered AdhE | 2 | 6.57 | fba, tpi, gap, pgk, gpm, eno | 0.010 | 31.0 |
| 11 | Tsac NADPH substrate-channeling | 2 | 6.57 | fba, tpi, gap, pgk, gpm, eno | 0.010 | 11.9 |
|  9 | Tsac fd adh | 2 | 6.54 | fba, tpi, gap, pgk, gpm, eno, aldh | 0.011 | — |
| 15 | Cth malate shunt FNORp adhp | 2 | 4.40 | upper glycolysis, pepck, mdh, mae, pfor, fnorp, aldhp, adhp | 0.010 | 20.2 |
|  3 | Tsac NADH+NADPH ethanologen fnorp | 2 | 3.11 | fba, tpi, gap, aldh, adhp | 0.296 | 100 |
|  6 | Cth malate shunt adhp | 3 | 2.78 | upper glycolysis, pepck, mdh, mae, pfor, rnf, aldh, adhp | 0.227 | 100 |
|  2 | Tsac NADH+NADPH ethanologen | 2 | 2.64 | fba, tpi, gap, bif-hyd, nfn, aldh, adhp | 0.430 | 47.2 |
|  7 | PDC ethanol pathway | 2 | 2.15 | fba, tpi, gap, adh | 0.640 | — |
|  8 | PDC Xhyd ethanol pathway | 2 | 1.61 | fba, tpi, gap, xhyd, adhp | 0.987 | 0.61 |
| 10 | Tsac NADH substrate-channeling | 2 | 0.92 | fba, tpi, gap, aldh_adh | 1.721 | — |
|  1 | Tsac NADH ethanologen | 2 | 0.77 | fba, tpi, gap, aldh, adh | 1.949 | — |
|  5 | Cth ethanologen WT adhE | 3 | 0.77 | fba, tpi, gap, aldh, adh | 1.949 | — |
| 12 | PDH ethanol pathway | 2 | 0.77 | fba, tpi, gap, aldh, adh | 1.949 | — |
| 13 | PFOR FNOR ethanol pathway | 2 | 0.77 | fba, tpi, gap, aldh, adh | 1.949 | — |
| 14 | PFOR FNOR Xhyd pathway | 2 | 0.66 | fba, tpi, gap, xhyd, aldh, adhp | 2.130 | 1.93 |

All fifteen pathways are thermodynamically feasible, including both malate-shunt modes.

### The bottleneck is shared between the two ends of the pathway

"Upper glycolysis" is the fba/tpi/gap segment, binding in all fifteen pathways. In most
of them the terminal aldehyde and alcohol dehydrogenase steps are binding as well: the
bottleneck is **shared** between the two ends, not localized to glycolysis.

The two ends limit together because they are coupled through a single cofactor ratio.
GAPDH reduces NAD⁺, so its driving force rises as NADH/NAD⁺ falls. ALDH and ADH oxidize
NADH, so theirs rises as NADH/NAD⁺ rises. When all three draw on the same pool no value
of the ratio satisfies both, and the optimum is the compromise at which they become
limiting together — NADH/NAD⁺ = 1.95 in the native *C. thermocellum* pathway.

### Excess driving force is trapped at PFOR

The native pathway is not short of driving force overall. At its optimum (MDF 0.77):

```
glk  28.81   fba   0.77 *  pgk   4.23   pyk  11.36   aldh  0.77 *
pgi   5.63   tpi   0.77 *  gpm   2.37   pfor  7.87   adh   0.77 *
pfk   9.37   gap   0.77 *  eno   4.45   rnf   3.05         (* = at the MDF)
```

PFOR runs at 7.87 kJ/mol, ten times the MDF. That surplus is **trapped**: RNF can only
discharge it into NADH/NAD⁺, the very ratio GAPDH needs kept low, so moving more
electrons into that pool relieves the terminal steps only by constraining GAPDH to the
same degree. The excess stays stranded upstream, contributing nothing to the achievable
titer.

### NfnAB unlocks it by separating the two redox ratios

NfnAB reduces NADP⁺ at the expense of reduced ferredoxin *and* NADH, which lets the two
pyridine nucleotide ratios be driven apart instead of held in compromise. With an
NADPH-linked AdhE, NADH/NAD⁺ falls to 0.010 — the 1:100 bound — while NADPH/NADP⁺ rises
to 31, a ~3,000-fold separation. GAPDH and the terminal reductions can then both be given
high driving force, using the ferredoxin surplus previously stranded at PFOR, and the MDF
rises from 0.77 to 6.57 kJ/mol.

The MDF column tracks the NADH/NAD⁺ column across the whole table. Every pathway that
decouples the terminal reductions from NADH — by NfnAB plus an NADPH-linked AdhE (mode 4),
a ferredoxin-linked ADH (mode 9), or channeled NADPH-linked activities (mode 11) — drives
NADH/NAD⁺ to ~0.01 and reaches ≈6.5 kJ/mol. Every pathway retaining NADH-linked terminal
steps settles near NADH/NAD⁺ ≈ 1.9 and ≈ 0.77 kJ/mol.

### Sensitivity to ethanol titer

| ethanol (mM) | M05 Cth WT adhE | M04 NFN + engineered AdhE |
|---:|---:|---:|
|   10 | 2.96 | 6.57 |
|  100 | 2.01 | 6.57 |
|  500 | 1.34 | 6.57 |
| 1000 | 1.06 | 6.57 |
| 1500 | 0.89 | 6.57 |
| 2000 | 0.77 | 6.57 |
| 2500 | 0.68 | 6.57 |

The MDF of the NFN pathway is independent of ethanol concentration, because its terminal
steps draw on a separate, highly reduced NADPH pool and have driving force to spare. In
the wild-type pathway the ALDH and ADH steps are part of the bottleneck, so every
increase in ethanol is paid for directly out of the pathway's driving force.

### Cofactor specificity of the ADH sets the titer response

`python mdf_vs_titer.py --modes 1 2 --name tsac_adh_cofactor` compares the two
*T. saccharolyticum* routes that differ only at the terminal reduction: mode 1, with
NADH-linked ALDH and ADH, against mode 2, which adds NfnAB to supply NADPH to an
NADPH-linked ADH (`results/figures/tsac_adh_cofactor.*`).

| ethanol (g/L) | MDF, NADH-ADH | MDF, NADPH-ADH | NADH/NAD⁺ (M01) | NADH/NAD⁺ (M02) |
|---:|---:|---:|---:|---:|
| 1 | 2.64 | 4.88 | 0.43 | 0.07 |
| 20 | 1.40 | 3.40 | 1.17 | 0.23 |
| 55 | 0.98 | 2.90 | 1.64 | 0.35 |
| 100 | 0.74 | 2.60 | 2.00 | 0.45 |

The NADPH-linked route holds ~2 kJ/mol more driving force at every titer. Put another
way, at 100 g/L it still has more driving force than the NADH-linked route has at 1 g/L.
An NADH-linked ADH has to share the NAD(H) pool with GAPDH, which wants it poised the
other way, so sustaining the terminal reduction against accumulating product drives
NADH/NAD⁺ from 0.43 to 2.00 and charges every increment to GAPDH. An NADPH-linked ADH
does not compete for that pool.

The same script handles any pathway pair, e.g. `--modes 5 4` for the native
*C. thermocellum* route against NfnAB with an engineered AdhE. This offers a thermodynamic rationale for
the finding of Fabri et al. (2025) that AdhE mutations which disrupt the ADH domain, and
so reduce NADH-linked ADH activity, dramatically increase ethanol tolerance. Those
mutations do not change AdhE's cofactor specificity — none of the variants showed
substantial NADPH-linked activity — but by removing the NADH-linked route they leave
ethanol production to an NADPH-linked ADH, moving the organism from mode 1 toward mode 2.
Note that the measured effect is on *tolerance*; titer data are not yet available, and
the two are not interchangeable.

### Every pathway in one figure

`python combined_figure.py` writes `combined_pathway_overview_part1.*` and `_part2.*` to
`results/figures/`: one row per pathway, sorted by MDF at the reference titer, split
seven rows and eight. Three panels across — the cumulative driving-force plot at 2 M
ethanol, the MDF against ethanol titer, and the three redox ratios overlaid on one log
axis against titer. The MDF panels share one vertical scale across every row and both
pages, so heights are directly comparable, and each ratio panel carries its legend above
the axes rather than inside, where a curve would eventually reach it.

Pages are 6 x 11 in — narrower than the paper they sit on, which keeps the panels close
to square and enlarges every label when the figure is scaled to a document's text width.

`--modes` selects a subset, `--split 0` puts it all on one page, and `--name` sets the
file stem. The six-pathway version used in the manuscript is

```bash
python combined_figure.py --modes 14 5 1 4 3 2 --split 0 --name selected_pathway_overview
```

which contrasts NADPH from NfnAB (M04, 6.57 kJ/mol) against NADPH from a
transhydrogenase (M14, 0.66) and the NADH-linked routes in between.

The top three rows (modes 11, 4 and 9) show the pattern the analysis turns on: the middle
panel is a flat line, and in the right panel NADH/NAD⁺ sits pinned on the 1:100 bound
while NADPH/NADP⁺ and Fd(red)/Fd(ox) run two to four orders of magnitude above it. The
bottom rows (modes 1, 13, 12, 5) show the opposite — a decaying MDF, and all available
ratios converging near 1.

Pass `--reuse-titration` to read the curves from `results/tables/cofactor_titration.csv`
rather than re-solving them; it is the same calculation, and it turns a ten-minute
rebuild into a few seconds.

### How each pathway re-poises its cofactors as product accumulates

`python cofactor_titration.py` writes one four-panel figure per pathway to
`results/figures/titration/`, showing the MDF and the NADH/NAD⁺, NADPH/NADP⁺ and
Fd(red)/Fd(ox) ratios against ethanol titer, with the full curves in
`results/tables/cofactor_titration.csv`. Ratio panels are on a log scale with dotted
lines at the 1:100 and 100:1 bounds; a pathway that lacks a pool leaves that panel empty
rather than dropping it, so the figures stay comparable.

Read together they separate the pathways into two kinds. Those whose terminal reductions
draw on NADH must walk NADH/NAD⁺ up as ethanol accumulates and lose driving force doing
it — modes 1, 5, 12 and 13 all run 0.20 → 2.00 while the MDF falls 3.59 → 0.74 over
0.1–100 g/L. Those that decouple the terminal steps hold NADH/NAD⁺ pinned at the 1:100
bound and their MDF does not move at all: modes 4, 9 and 11 are flat across the entire
range.

### Comparison with measured metabolome data

`python compare_to_tian2017.py` plots the measured and predicted NADH/NAD⁺ ratios side by
side on a common scale (`results/figures/tian2017_comparison.*`). The measurements are
from Tian et al. (2017), who added ethanol to a growing *C. thermocellum* culture at
~5 g/L/h and quantified intracellular metabolites by LC-MS against purified standards.

| ethanol (g/L) | measured NADH/NAD⁺ | predicted |
|---:|---:|---:|
| 0.5 | 0.30 | 0.33 |
| 9.2 | 0.45 | 0.90 |
| 21.7 | 0.54 | 1.20 |
| 33.6 | 0.94 | 1.39 |
| 48.4 | 1.39 | 1.57 |
| 55.3 | 1.73 | 1.64 |

Both rise about five-fold across the same range of titers, and converge at high ethanol.
At intermediate titers the measurement lags the prediction, which is the expected
direction of disagreement: MDF describes a cell operating exactly at its thermodynamic
optimum, which a real culture approaches only once thermodynamics, and not kinetics, is
the binding constraint.

The five-fold span is structural rather than fitted. NADH/NAD⁺ ∝ [ethanol]^⅓ follows from
the binding constraints, and the experiment spans 120-fold in titer, so the predicted
span is 120^⅓ = 4.93 regardless of parameters — enabling the dG′° confidence ellipsoid,
fixing CO₂ anywhere in 0.001–10 mM, fixing NAD⁺ anywhere in 0.01–1 mM, or removing the
ATP/ADP constraint all leave it unchanged. Only the dG′° treatment shifts the absolute
level, and by the same factor at both ends.

The close agreement at *low* titer is coincidence, and should not be claimed as support:
that sample was taken before ethanol was added, from a cell that is not thermodynamically
limited and has no reason to sit at its MDF optimum.

The same study independently identified GAPDH as the site of the bottleneck, by two routes
this analysis does not use: metabolites accumulated upstream of GAPDH and were depleted
downstream of it, and the purified *C. thermocellum* enzyme lost more than half its
activity at NADH/NAD⁺ = 0.2 and essentially all of it at 1.0, where the
*T. saccharolyticum* enzyme retained ~30%. Expressing the *T. saccharolyticum* `gapdh` in
*C. thermocellum* improved both ethanol tolerance and production. The measured ratio
crosses 1.0 at roughly 35 g/L ethanol — the point at which, in vitro, the native enzyme
stops working.

### Net ΔrG′° is identical within a stoichiometry class

Every 2-ATP pathway has a net Δ*r*G′° of −151.55 kJ/mol and both 3-ATP (RNF) pathways
−119.93 kJ/mol, as they must, since they catalyse the same overall conversion. The
cumulative Δ*r*G′ at the *optimum* still differs slightly between pathways, but only
through the ATP/ADP ratio — the one species in the net reaction that is not pinned by the
concentration bounds. Pathways that settle at ATP/ADP = 1.0 all end at exactly
−182.35 kJ/mol.

(With `dg_confidence > 0` this is no longer true: the optimizer picks a different Δ*r*G′°
assignment within the confidence ellipsoid for each pathway, shifting the net Δ*r*G′° by
as much as +28 kJ/mol. That is one reason the point estimates are used here.)

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
- Olson DG, Hörl M, Fuhrer T, Cui J, Zhou J, Maloney MI, Amador-Noguez D, Tian L, Sauer U,
  Lynd LR (2017). Glycolysis without pyruvate kinase in *Clostridium thermocellum*.
  *Metabolic Engineering* 39:169–180. doi:10.1016/j.ymben.2016.11.011
- Tian L, Perot SJ, Stevenson D, Jacobson T, Lanahan AA, Amador-Noguez D, Olson DG,
  Lynd LR (2017). Metabolome analysis reveals a role for glyceraldehyde 3-phosphate
  dehydrogenase in the inhibition of *C. thermocellum* by ethanol.
  *Biotechnol Biofuels* 10:276. doi:10.1186/s13068-017-0961-3
- Fabri JHTM, Pech-Canul A, Ziegler SJ, Burgin TE, Richardson ID, Maloney MI, Bomble YJ,
  Lynd LR, Olson DG (2025). The role of AdhE mutations in *Thermoanaerobacterium
  saccharolyticum*. *Journal of Bacteriology*. doi:10.1128/jb.00015-25
- Dash S, Olson DG, Joshua Chan SH, Amador-Noguez D, Lynd LR, Maranas CD (2019).
  Thermodynamic analysis of the pathway for ethanol production from cellobiose in
  *Clostridium thermocellum*. *Metabolic Engineering* 55:161–169.
  doi:10.1016/j.ymben.2019.06.006

## Citing this work

Please cite the paper:

> [CITATION TO BE ADDED]

Each GitHub release is archived at Zenodo. Cite the archive instead only if you need to
refer to the code or data specifically: <https://doi.org/10.5281/zenodo.23005004>. That
DOI covers all versions and always resolves to the most recent one; the Zenodo record for
an individual release carries its own DOI if you need to pin a particular version.

## License

MIT — see `LICENSE`.
