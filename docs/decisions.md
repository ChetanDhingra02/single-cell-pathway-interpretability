# Analysis decisions

## D1: Gene-set source

Use MSigDB Hallmark 2026.1.Hs from `raw/msigdb/h.all.v2026.1.Hs.symbols.gmt`. All 50 Hallmark pathways are considered at the start.

## D2: Cell QC

Keep cells with at least 200 detected genes, counted from raw counts. The main Kang dataset has 24,562 cells after this filter. A version without the cutoff is kept only for the QC sensitivity analysis.

## D3: Normalization

Kang and HIRISA are both normalized per cell to counts per 10,000 and then `log1p` transformed. Raw Kang counts are preserved in the `counts` layer.

## D4: Donor-held-out evaluation

Use 8 leave-one-donor-out folds. For fold `i`:

- test donor = donor `i`
- validation donor = the next donor in the sorted donor list
- training donors = the remaining six

No test-donor information is used for gene filtering, scaler fitting, early stopping or model selection.

## D5: Eligible cell types

A cell type is eligible for pathway-response analysis when at least 6 of 8 donors have at least 20 CTRL and at least 20 STIM cells.

The eligible Kang cell types are:

- B cells
- CD4 T cells
- CD14+ Monocytes
- Dendritic cells
- FCGR3A+ Monocytes
- NK cells

## D6: Model input genes

Within each donor-CV fold:

1. start from genes in the retained Hallmark pathways
2. require at least 1% detection in training-donor cells
3. retain pathways with at least 10 available genes
4. keep genes covered by the retained pathways

The scaling mean and SD are fitted on training donors only.

## D7: Neural models

The shared encoder is unrestricted. Models D and D2 use a non-negative masked decoder, with decoder weights outside the pathway mask fixed at zero.

Model E has the same architecture as D, but one random permutation of gene identities is shared across all pathway masks. This keeps pathway sizes and the pairwise overlap topology intact.

## D8: Network seeds

Donor-CV neural models use seeds 0, 1 and 2. Attribution is expensive, so it uses all three seeds for D and D2 and only seed 0 for E.

The final Kang-to-HIRISA transfer models use seed 0 only. Their training length is the median of the locally reproduced best epochs across the donor-CV runs.

## D9: Response scaling

Internal pathway responses are divided by the SD of the score over training-donor cells in the same fold.

External HIRISA responses are divided by a Kang reference SD calculated from the frozen final Kang representation. HIRISA statistics are never used to set the scale.

## D10: Attribution

Integrated Gradients uses a zero baseline in standardized expression space.

Attribution purity is:

`sum absolute attribution on member genes / sum absolute attribution on all input genes`

## D11: Ablation

For each pathway unit:

- zero the pathway's member genes
- measure the mean absolute change in that latent unit
- compare it with 20 equal-size random non-member gene sets

The reported ablation ratio is the pathway ablation effect divided by the mean random-ablation effect.

## D12: Interferon overlap decomposition

The measured Hallmark IFN sets contain:

- 72 shared genes
- 23 IFN-alpha-unique genes
- 123 IFN-gamma-unique genes

D2 replaces the original IFN-alpha and IFN-gamma pathways with shared and unique modules and leaves every other retained Hallmark pathway unchanged.

## D13: QC sensitivity

Repeat the model-D seed-0 analysis without the 200-gene cell cutoff. The result counts as stable when:

- eligible cell types are unchanged
- IFN-alpha response signs agree
- H1 verdicts for D and A are unchanged
- the median pathway-faithfulness shift is below 0.05

All checks passed in the local reproduction.

## D14: External validation hierarchy

Primary external analysis:

- HIRISA / GSE306664
- CD14 monocytes
- 5 donors
- none, IFN-alpha, IFN-beta and IFN-gamma
- IFN-alpha vs IFN-gamma as the primary specificity contrast

Secondary replication:

- B cells
- NK cells
- T cells

The primary CD14 result is not redefined after stage 2.

## D15: External model transfer

The final Kang models are trained only on genes available in both Kang and the HIRISA probe panel. The frozen shared list has 2,018 genes, and 2,016 of them pass the Kang detection and pathway rules and enter the final D/D2 models.

No HIRISA expression values are used for retraining, early stopping, pathway selection or scaler fitting.
