# External validation plan: HIRISA / GSE306664

The design and decision rules below were fixed before the HIRISA scores were calculated in the local reproduction.

## 1. Purpose

Kang contains IFN-beta stimulation only. That is enough to show interferon responsiveness and latent entanglement, but it cannot establish IFN-alpha versus IFN-gamma specificity.

HIRISA is an independent multi-interferon dataset. It is used to test whether splitting the overlapping IFN-alpha and IFN-gamma gene sets improves specificity.

## 2. Frozen representations

No model is retrained on HIRISA. The following Kang-derived representations are scored:

| Representation | Type-I unit | Type-II unit |
|---|---|---|
| D | Hallmark IFN-alpha | Hallmark IFN-gamma |
| D2 | IFNA_UNIQUE | IFNG_UNIQUE |
| A | mean of original IFN-alpha genes | mean of original IFN-gamma genes |
| A2 | mean of IFNA_UNIQUE genes | mean of IFNG_UNIQUE genes |

`IFN_SHARED` is reported descriptively for D2/A2.

## 3. Gene space and scaling

The final transfer models use genes shared between Kang and HIRISA:

- frozen shared list: 2,018 genes
- final model input after Kang filtering: 2,016 genes

HIRISA cells are normalized independently per cell with CP10k + `log1p` and then restricted to the exact final Kang gene order. The Kang scaler is applied to HIRISA; HIRISA mean and SD values are never used to rescale the model inputs.

External response magnitudes are divided by the corresponding Kang reference SD.

## 4. Stage 1: primary analysis

The primary cell population is CD14 monocytes.

Design:

- 5 donors
- unstimulated control
- IFN-alpha
- IFN-beta
- IFN-gamma
- 20 files total

The 20 files were chosen by experimental design, not by model performance, and all five donors are included. Each treatment is contrasted with the matched unstimulated sample from the same donor.

## 5. Specificity score

For donor `d`, unit `u` and treatment `X`:

`response(u | X) = [mean score under X - mean score under matched control] / SD_Kang(u)`

IFN-alpha is the primary type-I stimulus. IFN-beta is a secondary robustness stimulus.

For each donor:

`alpha_specificity = response(alpha | IFN-alpha) - response(alpha | IFN-gamma)`

`gamma_specificity = response(gamma | IFN-gamma) - response(gamma | IFN-alpha)`

`SEP = alpha_specificity + gamma_specificity`

A larger SEP means type-I and type-II interferon responses are better separated.

## 6. Primary decision rule

H5 compares D2 with D in CD14 monocytes using IFN-alpha versus IFN-gamma. It is supported when:

1. D2 has a larger SEP than D in at least 4 of 5 donors, and
2. median(D2 - D) is positive.

The same rule is evaluated for:

- D2 vs A
- D2 vs A2
- A2 vs A

### Credit rule

Let:

- `delta_network = median(SEP_D2 - SEP_D)`
- `delta_split = median(SEP_A2 - SEP_A)`

If `delta_split >= delta_network`, the improvement is credited mainly to the gene-set decomposition rather than to the neural network.

## 7. Stage 2: secondary replication

The same frozen models, scaler, response definition and comparison rule are applied separately to B cells, NK cells and T cells. Stage 2 is secondary and does not change the primary stage-1 verdict.

For B-cell donor 3491BW, the IFN-gamma sample has no same-pool control. The prespecified analysis uses the available donor-matched control, and a sensitivity analysis that excludes this donor is also reported.

## 8. Minimum-cell rule

A donor/treatment comparison is included only when the treated and control groups each have at least 20 QC-passing cells.

## 9. Prohibited adaptations

Once external scoring begins, none of the following is allowed:

- retraining on HIRISA
- HIRISA-based early stopping or hyperparameter tuning
- changing pathway definitions
- changing the IFN decomposition
- fitting a HIRISA scaler
- dropping donors because of their results
- choosing IFN-alpha versus IFN-beta after seeing the scores
- changing the donor-level decision threshold

## 10. External result

Primary CD14 monocytes:

- D2 vs D: 5/5 donors, median SEP difference +0.666
- D2 vs A: 5/5 donors, +1.519
- D2 vs A2: 0/5 donors, -0.996
- A2 vs A: 5/5 donors, +2.498

The primary D2 > D result is supported in CD14 monocytes. The neural model does not beat A2, though, so the improvement is credited to the shared/unique gene-set decomposition.

Stage 2:

- D2 vs D: 0/5 donors in B, NK and T cells
- D2 vs A2: 0/5 donors in B, NK and T cells
- A2 vs A: 5/5 donors in B, NK and T cells

In this dataset the D2 > D advantage appears only in monocytes, while the simpler split-set representation improves specificity in every cell type tested.
