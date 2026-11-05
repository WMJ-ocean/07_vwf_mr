# VWF and venous thromboembolism — Mendelian randomization analysis code

Analysis code for the manuscript:

> **Genetically predicted von Willebrand factor and shear-adhesion genes in venous thromboembolism: a Mendelian randomization study**

This repository contains the full analysis pipeline, derived instrumental-variable (IV) lists, and instrument-level results for a two-sample Mendelian randomization (MR) study that estimates the causal effects of (i) five blood-cell rheological traits and (ii) 24 shear-response / platelet-adhesion / coagulation genes on venous thrombotic outcomes (VTE, DVT, PE).

## Repository layout

```
.
├── code/          # analysis scripts (Python + R + bash)
├── instruments/   # derived instrumental variables (clumped IV lists, cis-eQTL instruments)
├── results/       # MR output tables (IVW / weighted median / MR-Egger / heterogeneity / pleiotropy / Steiger)
├── figures/       # study-design diagram and forest plots
└── README.md
```

Raw summary statistics and the LD reference panel are **not** stored here (see [Data sources](#data-sources)).

## Dependencies

- **R 4.5.2** with packages `TwoSampleMR`, `MendelianRandomization`, `dplyr`
- **Python 3.13** (standard library only; `pandas` optional for plotting)
- **PLINK 1.9** (used only for LD clumping)

## Reproducing the analysis

> **Note on paths.** The scripts were run on a specific Windows machine and contain a hard-coded base path `C:/Users/Hero/WorkBuddy/07_rheo_MR` (or `/d/R-4.5.2/...` for Rscript). Before re-running on your own machine, change the `BASE` / `RSCRIPT` variables at the top of each script to your local working directory.

1. **Obtain raw summary statistics** (links in [Data sources](#data-sources)) and place them under a working directory with the following structure:
   - `01_raw_exposure/` — Astle et al. blood-cell trait summary statistics (`hct`, `plt`, `mpv`, `rbc`, `pct`)
   - `02_raw_outcome/` — FinnGen R12 outcome files (`I9_VTE`, `I9_PHLETHROMBDVTLOW` = DVT, `I9_PULMEMB` = PE)
   - `tools/` — PLINK 1.9 and the 1000G Phase 3 EUR reference panel

2. **Extract instruments** for the five blood-cell traits:
   ```bash
   python code/01_extract_iv.py
   ```

3. **Strict LD clumping** (PLINK + 1000G Phase 3 EUR, r² < 0.001, 10 Mb window):
   ```bash
   bash code/07_ld_clump.sh
   ```

4. **Extract gene-level cis-eQTL instruments** (eQTLGen):
   ```bash
   python code/10_extract_gene_eqtl.py
   ```

5. **Run MR** (blood-cell traits → outcomes, clumped instruments):
   ```bash
   bash code/08_run_mr_clumped.sh
   ```

6. **Run gene-level MR** (single cis-eQTL per gene, Wald ratio + Steiger):
   ```bash
   python code/11_gene_mr.py
   ```

7. **Summarise and plot**:
   ```bash
   python code/05_summarize.py
   python code/06_plot_forest.py
   python code/13_plot_fig3.py
   ```

## Data sources

| Data | Consortium / reference | Access |
|---|---|---|
| Blood-cell traits (HCT, PLT, MPV, RBC, PCT) | Astle et al., *Cell* 2016;167:1415–1429 | `https://www.ebi.ac.uk/gwas/` (see original publication) |
| cis-eQTL instruments | eQTLGen consortium | `https://www.eqtlgen.org/` |
| VTE / DVT / PE outcomes | FinnGen release 12 | `https://www.finngen.fi/en/access_results` |
| 1000 Genomes Phase 3 EUR reference (LD clumping) | Zenodo record 7768714 | `https://doi.org/10.5281/zenodo.7768714` |

## License

MIT — see [LICENSE](LICENSE).
