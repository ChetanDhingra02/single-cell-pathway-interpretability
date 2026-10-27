# Data dictionary

## Processed Kang files

### `data/kang/kang_raw.h5ad`

24,562 cells x 15,706 genes after the main QC filter of 200 detected genes.

- `.X`: raw UMI counts
- `.obs`: donor, condition, cell_type, n_genes, n_counts

### `data/kang/kang_log1p.h5ad`

Same cells and genes as `kang_raw.h5ad`.

- `.X`: log1p(counts per 10,000)
- `layers["counts"]`: raw counts

### `data/kang/kang_log1p_noqc.h5ad`

24,673 cells before the 200-gene cell cutoff, used only for the QC sensitivity analysis.

## Kang observation columns

| Column | Meaning |
|---|---|
| `donor` | donor identifier, renamed from upstream `replicate` |
| `condition` | `ctrl` or `stim` |
| `cell_type` | upstream cell-type annotation |
| `n_genes` | number of genes with count > 0 |
| `n_counts` | total UMI count |

## Gene-set files

### `data/gene_sets/hallmark_membership.csv`

Columns:

- `pathway`
- `gene`

MSigDB Hallmark 2026.1.Hs membership, intersected with the measured Kang genes.

### `data/gene_sets/overlap_aware_membership.csv`

Same format, except that the original IFN-alpha and IFN-gamma pathways are replaced by:

- `IFN_SHARED`
- `IFNA_UNIQUE`
- `IFNG_UNIQUE`

## Donor-CV fold files

### `data/folds/fold_<i>.csv`

Columns:

- `cell_id`
- `donor`
- `split`

`split` is one of `train`, `val`, `test`.

### `data/folds/genes_fold_<i>.csv`

One `gene` column, in the exact model input order for fold `i`.

### `data/folds/pathways_fold_<i>.csv`

One `pathway` column, in latent-unit order for D/E.

### `data/folds/pathways_fold_<i>_overlap_aware.csv`

One `pathway` column, in latent-unit order for D2.

## Model files

Donor-CV model files are stored under `models/fold_<i>/` and are kept out of Git.

| File | Meaning |
|---|---|
| `scaler.npz` | training-donor mean, SD and gene order |
| `<M>_seed<s>.pt` | PyTorch state dictionary |
| `<M>_seed<s>_latent.npy` | latent scores for all cells |
| `<M>_seed<s>_units.csv` | latent-unit names |

Model codes:

- A = gene-set average
- B = PCA
- C = unrestricted autoencoder
- D = pathway-constrained autoencoder
- E = gene-identity-shuffled pathway control
- D2 = overlap-aware pathway autoencoder
- A2 = overlap-aware gene-set average
- mean = training-mean reconstruction baseline

## Final HIRISA transfer model files

Stored under `models/final_hirisa/`:

| File | Meaning |
|---|---|
| `D.pt` | final original pathway autoencoder |
| `D2.pt` | final overlap-aware autoencoder |
| `D_metadata.json` | gene order, pathway order, epoch count and training metadata |
| `D2_metadata.json` | same for D2 |
| `scaler.npz` | final Kang mean, SD and 2,016-gene order |
| `kang_reference_sd.csv` | Kang SD used to scale D/D2/A/A2 external responses |

## HIRISA metadata

### `data/hirisa/download_stage1.csv`

The 20 primary CD14 files.

### `data/hirisa/stage1_file_verification.csv`

Frozen file identifiers and verification information for stage 1.

### `data/hirisa/download_stage2.csv`

The 65 stage-2 files used for the B, NK and T analyses.

### `data/hirisa/download_stage2_verification.csv`

Frozen verification information for the stage-2 files.

### `data/hirisa/kang_hirisa_shared_genes.csv`

The 2,018 genes shared by the Kang analysis space and the HIRISA probe panel, before the final Kang detection and pathway filters.

## Key result columns

| Column | Definition |
|---|---|
| `difference` | held-out donor mean(STIM score) - mean(CTRL score) |
| `std_difference` | `difference` divided by training-donor score SD |
| `spearman` | Spearman correlation used for faithfulness or latent association |
| `attribution_purity` | fraction of total absolute Integrated Gradients attribution assigned to member genes |
| `attribution_enrichment` | attribution purity divided by member-gene input fraction |
| `ablation_ratio` | pathway-gene ablation effect divided by mean equal-size random non-member ablation effect |
| `jaccard` | pathway gene-set intersection / union |
| `alpha_specificity` | type-I alpha-unit response minus IFN-gamma alpha-unit response |
| `gamma_specificity` | IFN-gamma gamma-unit response minus type-I gamma-unit response |
| `SEP` | alpha-specificity + gamma-specificity |

## Main output files

- `results/results_summary.md`
- `results/low_faithfulness_units.md`
- `results/external/external_results.md`
- `results/tables/donor_cv_performance.csv`
- `results/tables/attribution_purity.csv`
- `results/tables/overlap_latent_association.csv`
- `results/tables/hirisa_specificity.csv`
- `results/tables/hirisa_stage2_specificity.csv`
