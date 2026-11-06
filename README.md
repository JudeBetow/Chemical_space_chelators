# Fragment-Based Analysis of Natural Product and Synthetic Chelators

A reproducible, RDKit-based cheminformatics workflow for comparing the BRICS fragment composition, enrichment, and fingerprint-based diversity of natural product chelators (NC) and synthetic chelators (SC).

## Project overview

This repository contains a Jupyter notebook that re-analyses the supplementary chelator dataset associated with Feeney *et al.* (2026), *Mapping the Chemical Space of Synthetic and Natural Product Chelators*. The workflow classifies compounds using the `NP info` annotation, standardises their SMILES, fragments the cleaned structures with BRICS, and compares the resulting fragment populations.

**Important:** This repository contains a secondary fragment-level analysis of the published supplementary dataset. It is not the original authors' descriptor-based chemical-space analysis and should not be interpreted as a reproduction of every result in the paper.

## Research objectives

1. Prepare a consistent set of organic ligand structures from the supplied chelator dataset.
2. Compare BRICS fragments observed in natural product and synthetic chelators.
3. Identify fragments statistically enriched in either class.
4. Compare fragment fingerprint diversity using Tanimoto distance.
5. Visualise fragment overlap, enrichment, frequency, and fingerprint-based PCA.

## Repository contents

```text
Chemical_space_chelators/
├── README.md
├── requirements.txt
├── .gitignore
├── fragment_analysis_NC_vs_SC.ipynb
├── ao6c06383_si_001.csv
├── plots/
│   ├── 01_raw_class_counts.png
│   ├── 02_preprocessing_funnel.png
│   ├── 03_venn_fragments.png
│   ├── 04a_top_nc_fragments.png
│   ├── 04b_top_sc_fragments.png
│   ├── 04c_discussed_nc_fragments.png
│   ├── 04d_discussed_sc_fragments.png
│   ├── 05_volcano.png
│   ├── 06_diversity_violin.png
│   ├── 07_frags_per_mol.png
│   ├── 08_pca_fragments.png
│   └── 09_tanimoto_heatmap.png
└── reports/
    ├── diversity_summary.csv
    ├── fragment_enrichment.csv
    ├── pca_coordinates.csv
    └── processed_compounds.csv
```

The `brics_cache.pkl` file from the working archive is intentionally not included. It is a regenerable intermediate cache, is not required to interpret the supplied results, and Python pickle files should not be loaded from untrusted sources. On a fresh run, the notebook can regenerate the cache.

## Dataset and classification

The input file, `ao6c06383_si_001.csv`, is the supplementary CSV associated with the Feeney *et al.* paper cited below. The notebook applies this operational classification:

- **NC (natural product chelator):** `NP info` is present and non-empty.
- **SC (synthetic chelator):** `NP info` is missing or empty.

This classification follows the annotation supplied in the input dataset. It does not independently establish the biosynthetic origin of every compound.

## Workflow

### 1. Data loading and classification

The supplementary CSV is loaded with pandas. The `NP info` field is used to assign each record to NC or SC.

### 2. SMILES preprocessing

The notebook applies the following steps:

1. Parse the input SMILES using RDKit.
2. Disconnect organometallic bonds using `MetalDisconnector`.
3. Remove atoms outside the notebook's predefined organic-atom set.
4. Retain the largest remaining molecular fragment.
5. Sanitise and canonicalise the resulting SMILES.
6. Remove duplicate canonical SMILES within each class.
7. Remove structures that occur in both classes.

The preprocessing decisions affect the structures entering the fragment analysis. In particular, retaining only the largest organic fragment is a deliberate simplification and may remove smaller organic components from some records.

### 3. BRICS fragmentation

Cleaned structures are fragmented using RDKit's BRICS implementation with a minimum fragment size of three heavy atoms. To make the calculation more robust, the notebook:

- skips molecules with more than 40 BRICS-cleavable bonds;
- applies a 60-second per-molecule timeout where `SIGALRM` is available;
- records processing status; and
- checkpoints completed results in `brics_cache.pkl`.

Molecules skipped or failing at this stage are excluded from downstream fragment analysis. Timeout behaviour is platform-dependent; the notebook notes that the signal-based timeout is unavailable on Windows.

### 4. Morgan fingerprints

Fragments are represented using Morgan fingerprints generated with:

- Radius: **2**
- Fingerprint size: **2,048 bits**

These fingerprints are used for Tanimoto comparisons and PCA.

### 5. Fragment enrichment

For each fragment, the notebook compares its occurrence counts in NC and SC using a two-sided Fisher's exact test. P-values are adjusted using the Bonferroni method. The notebook defines class-enriched fragments using both statistical significance (`adjusted p < 0.05`) and an absolute log2 enrichment threshold greater than 1 in the relevant direction.

### 6. Fragment diversity

Pairwise Tanimoto distances (`1 − similarity`) are calculated from Morgan fingerprints. To limit memory requirements, the diversity calculation samples up to 2,000 fragments per class using the fixed random seed.

### 7. Visualisation and export

The notebook generates figures for class counts, preprocessing, fragment overlap, enriched fragments, enrichment statistics, diversity, fragments per molecule, and PCA. It exports processed compounds, fragment-enrichment statistics, diversity summaries, and PCA coordinates as CSV files.

## Main results in the supplied outputs

The following values are taken from the notebook outputs and the CSV reports included in this repository.

| Measure | Result |
|---|---:|
| Raw input records | 35,548 |
| Records after SMILES parsing | 35,532 |
| Records after within-class deduplication | 35,423 |
| Records after cross-class duplicate removal | 35,409 |
| Compounds retained after BRICS processing | 35,403 |
| Natural product chelators (NC) | 2,703 |
| Synthetic chelators (SC) | 32,700 |
| Unique NC fragments | 2,980 |
| Unique SC fragments | 29,569 |
| Shared fragments | 776 |
| NC-enriched fragments | 161 |
| SC-enriched fragments | 24 |
| Mean NC pairwise Tanimoto distance | 0.9149 |
| Mean SC pairwise Tanimoto distance | 0.9256 |

The reported mean Tanimoto distances are based on the notebook's sampled pairwise-distance calculation. A larger distance indicates lower average fingerprint similarity within the sampled fragment set; it should not be interpreted as a direct measure of whole-molecule chemical-space area.

## Installation

Python 3.10 or later is recommended. RDKit is most straightforward to install through conda-forge.

### Option A: Conda (recommended)

```bash
conda create -n chelator-fragments -c conda-forge python=3.11 rdkit jupyterlab ipykernel numpy pandas matplotlib seaborn matplotlib-venn scipy statsmodels scikit-learn
conda activate chelator-fragments
jupyter lab
```

### Option B: pip

Create and activate a virtual environment, then install the packages:

```bash
python -m venv .venv
```

On Windows:

```powershell
.venv\Scripts\activate
```

On macOS/Linux:

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

If RDKit installation fails with pip on your platform, use the conda installation above.

## Running the notebook

1. Clone or download this repository.
2. Install the dependencies.
3. Start JupyterLab from the repository's root directory:

   ```bash
   jupyter lab
   ```

4. Open `fragment_analysis_NC_vs_SC.ipynb`.
5. Run the notebook cells from top to bottom.

The notebook expects `ao6c06383_si_001.csv` in the repository root. It creates or updates the `plots/` and `reports/` directories and writes the BRICS checkpoint to the root directory. The first BRICS run may take some time; subsequent runs can reuse the checkpoint if it exists and its parameters match.

For a clean re-analysis, remove the generated `brics_cache.pkl` before running the notebook. The notebook will then recompute the BRICS fragments.

## Output files

- `reports/processed_compounds.csv` — processed compound records, canonical SMILES, BRICS fragments, and BRICS status.
- `reports/fragment_enrichment.csv` — fragment counts, frequencies, enrichment ratios, Fisher-test p-values, adjusted p-values, and significance flags.
- `reports/diversity_summary.csv` — summary statistics for NC and SC fragment Tanimoto distances.
- `reports/pca_coordinates.csv` — PCA coordinates and class labels for the fragments included in the PCA.
- `plots/` — figures generated by the notebook.

## Interpretation and limitations

- NC/SC labels are inherited from the source dataset's `NP info` field.
- The results depend on the preprocessing choices, including metal disconnection, atom filtering, largest-fragment selection, and duplicate removal.
- BRICS fragmentation is rule-based and represents one particular way to decompose molecular structures; it does not uniquely define chemically meaningful pharmacophores or metal-binding motifs.
- The two classes are strongly imbalanced in size. The enrichment analysis uses class-normalised frequencies and Bonferroni-adjusted Fisher tests, but results should still be interpreted in light of class size and fragment occurrence.
- Morgan fingerprints encode local substructural environments. Their Tanimoto distances and PCA projection are fragment-level representations, not a direct comparison of the complete molecular structures.
- The diversity calculation samples up to 2,000 fragments per class, so the reported means are sample-based.
- The notebook's 60-second timeout relies on `SIGALRM` and may not operate on Windows.
- The supplementary CSV contains fields and annotations inherited from the source publication; consult the original article and supporting information for the dataset's full provenance and limitations.

## Citation

If you use the source dataset, please cite the original publication:

Feeney SV, Matthews ST, Codd R. Mapping the Chemical Space of Synthetic and Natural Product Chelators. *ACS Omega*. 2026;11(35):53318–53329. https://doi.org/10.1021/acsomega.6c06383

If you use or adapt this notebook, please also cite this repository or the version of the notebook you used. Add a repository DOI (for example, through Zenodo) if you want a persistent, citable release.

## Reproducibility

The notebook sets Python and NumPy random seeds to `42`. Package versions are not pinned in the notebook, so exact numerical reproducibility across different RDKit, NumPy, scikit-learn, and plotting-library versions is not guaranteed. For a publication or formal release, record the tested package versions and the Git commit hash used to generate the final outputs.

## Author(s)

**Bakary N’tji Diallo, **Jude Y. Betow**, Boris D. Bekono, Solange A. Tanyi, José L. Medina-Franco5*, Fidele Ntie-Kang**. Further exploration of the chemical space of synthetic and natural chelators: scaffold-, fragment- and functional groups-based analysis (2026)

Research interests: cheminformatics, computational chemistry, natural product chemistry, and AI-assisted drug discovery.

GitHub: [@JudeBetow](https://github.com/JudeBetow)
