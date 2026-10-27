# Limitations

## Architecture and pathway labels

The decoder is pathway constrained, but the encoder is not, so a latent unit can depend on genes outside the pathway named in its label. Pathway names are therefore treated as hypotheses about what a latent variable means, not as guaranteed mechanistic identities.

## Pathway overlap

Hallmark pathways overlap substantially. The overlap makes pathway scores statistically dependent and complicates interpretation, so the pathway-pair analyses are descriptive and do not treat the 1,081 pairs as independent observations in any inferential test.

## Gene-set averages are strong baselines

Model A also reproduces the IFN-alpha response sanity check, and in HIRISA the overlap-aware average A2 generalizes more consistently than D2. The project does not claim that the neural model is uniformly better than simpler pathway scoring.

## Kang perturbation design

Kang contains IFN-beta stimulation only. It can test donor robustness and interferon entanglement, but IFN-alpha versus IFN-gamma specificity needs the independent HIRISA analysis.

## External validation size

The primary HIRISA test has five donors, so the external decision rule relies on donor-level consistency and median differences instead of asymptotic p-values. The result is external replication evidence, not a definitive population-level estimate.

## Cross-platform transfer

Kang and HIRISA differ in assay platform, donor population and stimulation design, so absolute score levels are not assumed to carry over. External comparisons use within-donor treatment-versus-control contrasts and Kang-derived reference SDs.

## Final network seed

The final external-transfer D and D2 models use one fixed neural-network seed, with the training length taken from donor CV. The external neural result therefore says nothing about across-seed uncertainty in the final transferred model.

## Cell-type differences

The stage-2 B, NK and T results are descriptive replications. The project does not infer a biological mechanism from D2 > D appearing in monocytes but not in the other cell types tested.

## Attribution methods

Integrated Gradients depends on the choice of baseline, and this project uses zero in standardized expression space. Attribution is one diagnostic of model dependence, not a uniquely correct biological explanation.

## Ablation

Setting genes to zero in standardized expression space perturbs the model; it is not a biological knockout. Ablation ratios show how much the model depends on those inputs and should not be read as causal gene effects in cells.

## Hallmark-only input universe

The pathway models use only genes covered by the retained Hallmark pathways after fold-specific filtering. The findings describe pathway-constrained representations within that input universe, not the whole transcriptome.

## Numerical reproducibility

Neural early stopping and downstream probe classifiers can differ slightly across PyTorch, BLAS and CPU environments. The project's conclusions are donor-level and qualitative, and they do not require bit-identical neural weights across machines.
