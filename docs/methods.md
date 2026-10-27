# Methods overview

## Kang preprocessing

The Kang single-cell dataset is loaded from the supplied `.h5ad` file. Only the upstream donor, condition and cell-type labels are kept; upstream embeddings and clustering outputs are not used.

Cells with fewer than 200 detected genes are removed. Counts are normalized per cell to 10,000 total counts and `log1p` transformed, and the raw counts are kept in a separate layer.

## Hallmark pathways

MSigDB Hallmark 2026.1.Hs gene sets are intersected with the genes measured in Kang. A pathway needs at least 10 available genes to enter a model. For donor CV, gene detection filtering is fit on training donors only.

## Cross-validation

Eight leave-one-donor-out folds are used. Each fold has six training donors, one validation donor and one held-out test donor.

The validation donor is used for neural early stopping. The held-out donor is used only for the final fold evaluation.

## Models

### A: gene-set average

For each pathway, average the standardized expression of its member genes.

### B: PCA

Fit PCA on training-donor cells with as many components as there are retained pathway units.

### C: unrestricted autoencoder

Encoder:

`genes -> 256 ReLU hidden units -> K latent units`

Decoder:

`K latent units -> genes`

No pathway mask is applied.

### D: pathway-constrained autoencoder

The encoder is unrestricted. The decoder is masked so that latent unit `k` can reconstruct only genes in pathway `k`, and the allowed decoder weights are constrained to be non-negative.

### E: shuffled-gene control

Same architecture as D. One column permutation is applied to the pathway mask, which randomizes gene identities while preserving:

- pathway sizes
- all pairwise intersections
- the full Jaccard-overlap topology

### D2: overlap-aware autoencoder

Same architecture as D, with the original IFN-alpha and IFN-gamma pathways replaced by shared and unique modules.

## Training

Models are trained with Adam on mean squared reconstruction error. Donor-CV neural models use three random seeds, and early stopping uses validation MSE.

## Held-out performance

Reconstruction MSE is evaluated on the held-out donor.

Cell-type and condition information are measured with logistic-regression probes trained on the training-donor representations and tested on the held-out donor. The probes measure information content and play no part in training the autoencoders.

## Pathway response

For each held-out donor and eligible cell type:

`response = mean(STIM score) - mean(CTRL score)`

To make units comparable:

`standardized response = response / training-donor SD`

## Faithfulness

### Gene-set correlation

For each pathway unit, calculate the Spearman correlation between the latent score and the simple average of the unit's own member genes on held-out cells.

### Integrated Gradients

Integrated Gradients is calculated for the encoder output of each pathway unit. Attribution purity is the fraction of absolute attribution assigned to member genes.

### Ablation

Member genes are set to zero in standardized input space, and the resulting change in the latent unit is compared with the change from 20 equal-size random non-member gene sets.

## Pathway overlap

For every retained pathway pair, calculate the Jaccard overlap and the absolute Spearman correlation of the two scores on held-out cells. The analysis is repeated on CTRL-only cells.

## QC sensitivity

The main model-D seed-0 analysis is repeated on the unfiltered Kang cells, before the 200-gene cutoff. Predefined stability checks compare cell-type eligibility, IFN-alpha signs, H1 verdicts and median faithfulness.

## HIRISA transfer

The final D/D2 Kang models are trained on the 2,016-gene input space shared with HIRISA after the final Kang filters.

HIRISA cells are normalized with CP10k + `log1p`, put in the same gene order and standardized with the frozen Kang scaler. No model fitting uses HIRISA data.

Specificity is evaluated within donor, using treatment-versus-control contrasts and Kang reference SDs.
