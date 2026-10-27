# Analysis hypotheses

This file lists the project's hypotheses and decision criteria, and separates the confirmatory checks from descriptive follow-up analyses.

## H1: Donor robustness

The pathway-constrained autoencoder (model D) should preserve the main interferon response in donors it has never seen.

For each eligible cell type and held-out donor:

- response = mean(STIM score) − mean(CTRL score)
- the response is divided by the SD of that score in the fold's training-donor cells
- neural-network seeds are averaged within donor before donor-level summaries

Criterion: for every eligible cell type, the `HALLMARK_INTERFERON_ALPHA_RESPONSE` unit must be positive in all eligible donors except at most one, and its median donor response must rank in the top three retained pathways.

The same check is applied to the simple gene-set average (A). If A also meets the criterion, H1 is treated as a sanity check rather than as evidence that the neural model is better.

## H2: Pathway-label faithfulness is imperfect

A pathway name attached to a latent unit does not guarantee that the encoder computes that unit mainly from genes in the named pathway.

Faithfulness is evaluated with:

1. Spearman correlation between each latent unit and the average of its own pathway genes
2. Integrated Gradients attribution purity
3. pathway-gene ablation compared with equal-size random non-member ablation

Criterion: H2 is supported if at least 25% of retained model-D pathways have a median attribution purity below 0.5.

Pathways with a median ablation ratio <= 1 are reported separately as faithfulness warnings.

## H3: Pathway overlap is associated with latent correlation

Pathway pairs with more gene-set overlap are expected to have higher absolute latent-score correlation.

For every pathway pair:

- overlap = Jaccard similarity of the pathway gene sets
- latent association = absolute Spearman correlation of the corresponding scores on held-out donor cells

The association is also calculated on CTRL-only cells.

Criterion: H3 is supported if the Spearman correlation between pathway Jaccard overlap and absolute latent correlation is positive for model D.

Pathway pairs are not independent observations, so this is treated as a descriptive association rather than a confirmatory p-value test.

Model E is a gene-identity-randomized control. One shared column permutation is applied to the pathway mask, which randomizes gene identities while keeping pathway sizes and the full pairwise overlap topology.

## H4: IFN-alpha / IFN-gamma entanglement

The original Hallmark IFN-alpha and IFN-gamma latent units are expected to be substantially correlated, because the two pathway definitions share many genes and interferon stimulation activates both.

Criterion: H4 is supported if the median held-out-cell Spearman correlation between the two model-D units is at least 0.5.

Kang contains IFN-beta stimulation only, so it can show entanglement but cannot establish IFN-alpha versus IFN-gamma specificity.

## H5: External overlap-decomposition test

The overlap-aware representation splits the original interferon pathways into:

- `IFN_SHARED`
- `IFNA_UNIQUE`
- `IFNG_UNIQUE`

The primary external test uses HIRISA CD14 monocytes and compares IFN-alpha with IFN-gamma.

For each donor, the separation index is:

`SEP = alpha_specificity + gamma_specificity`

where:

- `alpha_specificity = response(alpha | IFN-alpha) - response(alpha | IFN-gamma)`
- `gamma_specificity = response(gamma | IFN-gamma) - response(gamma | IFN-alpha)`

Primary criterion: D2 must have a greater SEP than D in at least 4 of 5 donors, and the median donor-level difference must be positive.

The same donor-level rule is also applied to:

- D2 vs A
- D2 vs A2
- A2 vs A

If the split gene-set average A2 improves over A at least as much as D2 improves over D, the gain is attributed mainly to the gene-set decomposition rather than to the neural network.

B, NK and T cells are secondary replication analyses and do not change the primary CD14 verdict.
