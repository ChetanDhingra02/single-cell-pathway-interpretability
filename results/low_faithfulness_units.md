# Low-faithfulness units: what are they computed from? (model D)

Held-out donors; medians over folds with network seeds averaged within each fold.

Column guide: `attribution_purity` is the share of absolute Integrated Gradients attribution on the unit's own genes; `pathway_input_fraction` is the fraction of input genes belonging to that pathway; `other_ifn_attribution_share` is attribution assigned to IFN-alpha/gamma genes outside the unit's own set; `median_ablation_ratio` compares ablation of the pathway genes with equal-size random non-member gene sets; `median_spearman` is the correlation between the latent unit and its own gene-set average.

| unit                      |   attribution_purity |   pathway_input_fraction |   attribution_enrichment |   other_ifn_attribution_share |   other_ifn_input_fraction |   median_ablation_ratio |   folds_ratio_le_1 |   median_spearman |   B cells |   CD14+ Monocytes |   CD4 T cells |   Dendritic cells |   FCGR3A+ Monocytes |   NK cells |   own_genes_in_model |   own_genes_also_ifn_posthoc |
|:--------------------------|---------------------:|-------------------------:|-------------------------:|------------------------------:|---------------------------:|------------------------:|-------------------:|------------------:|----------:|------------------:|--------------:|------------------:|--------------------:|-----------:|---------------------:|-----------------------------:|
| IL6_JAK_STAT3_SIGNALING   |                0.244 |                    0.028 |                    8.588 |                         0.146 |                      0.08  |                   8.006 |                  0 |            -0.178 |    -0.125 |            -1.867 |        -0.035 |            -1.331 |              -1.108 |     -0.106 |                   61 |                           19 |
| APICAL_SURFACE            |                0.126 |                    0.006 |                   22.552 |                         0.072 |                      0.089 |                  29.657 |                  0 |             0.142 |     0.171 |             1.732 |         0.02  |             1.061 |               1.038 |      0.038 |                   12 |                            1 |
| ESTROGEN_RESPONSE_EARLY   |                0.192 |                    0.027 |                    7.224 |                         0.072 |                      0.089 |                   7.909 |                  0 |             0.099 |     0.158 |            -0.379 |         0.035 |            -0.221 |              -0.38  |      0.025 |                   56 |                            2 |
| UNFOLDED_PROTEIN_RESPONSE |                0.245 |                    0.041 |                    5.971 |                         0.068 |                      0.088 |                   6.623 |                  0 |             0.158 |    -0.067 |            -0.422 |        -0.002 |            -0.194 |              -0.432 |      0.144 |                   89 |                            3 |
| BILE_ACID_METABOLISM      |                0.197 |                    0.017 |                   11.677 |                         0.064 |                      0.089 |                  17.179 |                  0 |             0.212 |     0.098 |             0.588 |         0.075 |             0.452 |               1.006 |      0.076 |                   37 |                            1 |
| KRAS_SIGNALING_DN         |                0.261 |                    0.01  |                   26.467 |                         0.095 |                      0.088 |                  61.3   |                  0 |             0.298 |     0.038 |             1.484 |         0.091 |             0.953 |               0.796 |      0.007 |                   20 |                            3 |
| INTERFERON_ALPHA_RESPONSE |                0.596 |                    0.042 |                   14.327 |                         0.1   |                      0.048 |                  52.712 |                  0 |             0.906 |     1.942 |             1.847 |         1.682 |             2.266 |               1.761 |      1.746 |                   89 |                           89 |
| MEDIAN OVER ALL UNITS (D) |                0.372 |                    0.038 |                   10.498 |                         0.062 |                      0.088 |                  30.058 |                  0 |             0.673 |     0.005 |            -0.174 |         0.007 |            -0.143 |              -0.049 |      0.006 |                  nan |                          nan |

## Top input genes per unit

Genes most often appearing in the top-20 absolute attribution values across donor folds and network seeds.

### IL6_JAK_STAT3_SIGNALING

Category counts among these 20: {'own member': 12, 'other': 7, 'IFN-alpha/gamma member': 1}

| gene   |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category               |   runs_total |
|:-------|------------------:|------------:|-----------------------:|:-----------------------|-------------:|
| CD14   |                24 |      4.4167 |                 0.4814 | own member             |           24 |
| CXCL3  |                24 |      5.9583 |                 0.5005 | own member             |           24 |
| IL1B   |                23 |      4.0435 |                 0.5245 | own member             |           24 |
| CD36   |                21 |      7.2857 |                 0.3863 | own member             |           24 |
| NFKBIA |                21 |     13.2857 |                 0.2963 | IFN-alpha/gamma member |           24 |
| CD9    |                20 |     10.75   |                 0.3305 | own member             |           24 |
| CFLAR  |                19 |     12      |                 0.3269 | other                  |           24 |
| STAT1  |                17 |      3.4706 |                 0.5293 | own member             |           24 |
| PLEK   |                17 |      7.2353 |                 0.416  | other                  |           24 |
| CCL4   |                17 |     13.4118 |                 0.3334 | other                  |           24 |
| LTB    |                16 |      7.0625 |                 0.411  | own member             |           24 |
| NINJ1  |                16 |      8.5625 |                 0.3915 | other                  |           24 |
| SOCS1  |                15 |      1      |                 0.8301 | own member             |           24 |
| IL2RG  |                15 |      9.0667 |                 0.3916 | own member             |           24 |
| CXCL1  |                15 |      9.1333 |                 0.3341 | own member             |           24 |
| DRAM1  |                15 |     13.4    |                 0.3411 | other                  |           24 |
| IL1R2  |                14 |     14.7143 |                 0.2297 | own member             |           24 |
| SOCS3  |                14 |     14.9286 |                 0.2379 | own member             |           24 |
| KLF6   |                14 |     15.6429 |                 0.3066 | other                  |           24 |
| SAT1   |                12 |     16.1667 |                 0.3064 | other                  |           24 |

### APICAL_SURFACE

Category counts among these 20: {'other': 16, 'own member': 4}

| gene    |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category   |   runs_total |
|:--------|------------------:|------------:|-----------------------:|:-----------|-------------:|
| HSPB1   |                24 |      1.8333 |                 1.3395 | own member |           24 |
| LYN     |                24 |      4.4583 |                 0.8885 | own member |           24 |
| IL2RG   |                20 |      1.2    |                 1.8115 | own member |           24 |
| B4GALT1 |                20 |      6.8    |                 0.7065 | own member |           24 |
| PKM     |                17 |      3.8824 |                 1.0146 | other      |           24 |
| SLC16A3 |                17 |      5.8235 |                 0.8137 | other      |           24 |
| IER3    |                17 |      6.8235 |                 0.7514 | other      |           24 |
| TALDO1  |                16 |      8.625  |                 0.6157 | other      |           24 |
| IL13RA1 |                16 |     12.625  |                 0.467  | other      |           24 |
| ALDOA   |                15 |      9.2667 |                 0.6024 | other      |           24 |
| TPI1    |                15 |     11.1333 |                 0.54   | other      |           24 |
| VCAN    |                15 |     12.4667 |                 0.4852 | other      |           24 |
| ENO1    |                14 |     14.7143 |                 0.4233 | other      |           24 |
| QSOX1   |                13 |     14.9231 |                 0.4018 | other      |           24 |
| CD44    |                12 |     11.3333 |                 0.5531 | other      |           24 |
| GPX1    |                11 |     11.9091 |                 0.4089 | other      |           24 |
| HSP90B1 |                 9 |     13      |                 0.3616 | other      |           24 |
| CALR    |                 9 |     15      |                 0.332  | other      |           24 |
| LAPTM5  |                 7 |      8.7143 |                 0.5134 | other      |           24 |
| ARPC3   |                 7 |     15.8571 |                 0.2781 | other      |           24 |

### ESTROGEN_RESPONSE_EARLY

Category counts among these 20: {'own member': 11, 'other': 9}

| gene    |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category   |   runs_total |
|:--------|------------------:|------------:|-----------------------:|:-----------|-------------:|
| SNX10   |                24 |      6.25   |                 0.4096 | other      |           24 |
| RAB31   |                21 |     11.5238 |                 0.3149 | own member |           24 |
| PMAIP1  |                18 |      3.7778 |                 0.6323 | own member |           24 |
| ELOVL5  |                18 |     12      |                 0.2743 | own member |           24 |
| BHLHE40 |                17 |      6.6471 |                 0.4229 | own member |           24 |
| FKBP4   |                15 |      6.2    |                 0.4047 | own member |           24 |
| FOS     |                15 |     10.1333 |                 0.3054 | own member |           24 |
| ANXA2   |                14 |      4.9286 |                 0.5015 | other      |           24 |
| S100A11 |                13 |      4.7692 |                 0.5318 | other      |           24 |
| CD44    |                13 |     10.1538 |                 0.3951 | own member |           24 |
| S100A4  |                12 |      8.4167 |                 0.4311 | other      |           24 |
| FDFT1   |                12 |      9.9167 |                 0.2742 | own member |           24 |
| CXCL16  |                11 |      5.3636 |                 0.5177 | other      |           24 |
| ANXA5   |                11 |      8.5455 |                 0.4376 | other      |           24 |
| ATF5    |                11 |     10.5455 |                 0.4046 | other      |           24 |
| GAPDH   |                11 |     13.4545 |                 0.3549 | other      |           24 |
| XBP1    |                10 |      9.2    |                 0.2985 | own member |           24 |
| GSTO1   |                10 |     13.1    |                 0.3258 | other      |           24 |
| IGFBP4  |                 9 |     10.1111 |                 0.2495 | own member |           24 |
| TIPARP  |                 9 |     13      |                 0.2294 | own member |           24 |

### UNFOLDED_PROTEIN_RESPONSE

Category counts among these 20: {'own member': 13, 'other': 7}

| gene    |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category   |   runs_total |
|:--------|------------------:|------------:|-----------------------:|:-----------|-------------:|
| RPS14   |                24 |      4.5417 |                 0.4917 | own member |           24 |
| EIF4A2  |                24 |      5.7917 |                 0.4369 | own member |           24 |
| HSPA5   |                23 |      1      |                 1.1084 | own member |           24 |
| XBP1    |                23 |      2.913  |                 0.7041 | own member |           24 |
| IFIT1   |                22 |      9.6364 |                 0.302  | own member |           24 |
| CEBPG   |                22 |     11.7273 |                 0.2733 | own member |           24 |
| DDIT4   |                20 |      4.2    |                 0.5473 | own member |           24 |
| EEF2    |                20 |      8.8    |                 0.3219 | own member |           24 |
| MAFB    |                17 |      6.8235 |                 0.4595 | other      |           24 |
| DNAJB9  |                16 |     11.25   |                 0.2534 | own member |           24 |
| FCER1G  |                16 |     12.5    |                 0.3061 | other      |           24 |
| HSP90B1 |                15 |     11.4667 |                 0.2613 | own member |           24 |
| DNAJC3  |                15 |     11.6    |                 0.2709 | own member |           24 |
| CTSB    |                15 |     11.8    |                 0.3166 | other      |           24 |
| GOSR2   |                15 |     12.6667 |                 0.2359 | own member |           24 |
| CD63    |                13 |      6.7692 |                 0.4631 | other      |           24 |
| FCN1    |                13 |     14.7692 |                 0.2736 | other      |           24 |
| ANXA5   |                12 |     14.75   |                 0.3093 | other      |           24 |
| SPCS1   |                11 |     10.3636 |                 0.2284 | own member |           24 |
| CTSL    |                11 |     13      |                 0.3144 | other      |           24 |

### BILE_ACID_METABOLISM

Category counts among these 20: {'own member': 12, 'other': 8}

| gene    |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category   |   runs_total |
|:--------|------------------:|------------:|-----------------------:|:-----------|-------------:|
| SOD1    |                21 |      3.8571 |                 0.7382 | own member |           24 |
| OPTN    |                21 |      4.1905 |                 0.6509 | own member |           24 |
| TSPO    |                16 |     16.125  |                 0.2554 | other      |           24 |
| PRDX5   |                15 |      8.8    |                 0.3195 | own member |           24 |
| AQP9    |                14 |      6.1429 |                 0.4149 | own member |           24 |
| PRDX1   |                14 |     13.9286 |                 0.3038 | other      |           24 |
| GPX4    |                13 |     11.3846 |                 0.3682 | other      |           24 |
| GSTK1   |                12 |      3.9167 |                 0.6082 | own member |           24 |
| GPX1    |                12 |      5.8333 |                 0.5394 | other      |           24 |
| IDI1    |                12 |     10.3333 |                 0.2839 | own member |           24 |
| TXN     |                12 |     10.5833 |                 0.3838 | other      |           24 |
| LCK     |                11 |      3.2727 |                 0.5814 | own member |           24 |
| IDH2    |                11 |      5.6364 |                 0.5438 | own member |           24 |
| ISOC1   |                11 |      6.7273 |                 0.2899 | own member |           24 |
| GCLM    |                11 |     10.8182 |                 0.3438 | own member |           24 |
| CH25H   |                10 |      5.1    |                 0.4363 | own member |           24 |
| PKM     |                10 |      9      |                 0.4304 | other      |           24 |
| ECH1    |                10 |      9.2    |                 0.2083 | other      |           24 |
| ALDH9A1 |                10 |     10.5    |                 0.2803 | own member |           24 |
| TALDO1  |                10 |     12.3    |                 0.3454 | other      |           24 |

### KRAS_SIGNALING_DN

Category counts among these 20: {'own member': 10, 'other': 8, 'IFN-alpha/gamma member': 2}

| gene   |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category               |   runs_total |
|:-------|------------------:|------------:|-----------------------:|:-----------------------|-------------:|
| TGM1   |                24 |      1.5417 |                 0.7383 | own member             |           24 |
| SNN    |                24 |      2.2083 |                 0.6579 | own member             |           24 |
| CCNA1  |                24 |      2.75   |                 0.5968 | own member             |           24 |
| CD80   |                24 |      3.875  |                 0.4564 | own member             |           24 |
| RSAD2  |                24 |      6.5    |                 0.1394 | own member             |           24 |
| BARD1  |                24 |      7.125  |                 0.1287 | own member             |           24 |
| IFI44L |                22 |      9.5455 |                 0.0917 | own member             |           24 |
| SKIL   |                21 |      5.3333 |                 0.2291 | own member             |           24 |
| JUNB   |                21 |     10.1905 |                 0.0753 | other                  |           24 |
| RHOA   |                18 |     12.5    |                 0.0484 | other                  |           24 |
| IL2RG  |                16 |     14.25   |                 0.0357 | other                  |           24 |
| MX1    |                14 |      9.9286 |                 0.0621 | own member             |           24 |
| HSPA5  |                13 |     14.4615 |                 0.0329 | other                  |           24 |
| FKBP1A |                12 |     13.4167 |                 0.0494 | other                  |           24 |
| PLEK   |                12 |     16.1667 |                 0.0319 | other                  |           24 |
| SOD1   |                11 |     13.3636 |                 0.039  | other                  |           24 |
| NINJ1  |                11 |     16.9091 |                 0.0294 | other                  |           24 |
| GTF3C5 |                 9 |     11.8889 |                 0.0476 | own member             |           24 |
| SOCS1  |                 8 |     16      |                 0.0306 | IFN-alpha/gamma member |           24 |
| SOD2   |                 7 |     15.4286 |                 0.0359 | IFN-alpha/gamma member |           24 |

### INTERFERON_ALPHA_RESPONSE

Category counts among these 20: {'own member': 20}

| gene    |   n_runs_in_top20 |   mean_rank |   mean_abs_attribution | category   |   runs_total |
|:--------|------------------:|------------:|-----------------------:|:-----------|-------------:|
| OAS1    |                24 |      1.25   |                 0.5015 | own member |           24 |
| MX1     |                24 |      1.875  |                 0.4759 | own member |           24 |
| IFI44L  |                24 |      3.7917 |                 0.417  | own member |           24 |
| IFIT3   |                24 |      4.8333 |                 0.3841 | own member |           24 |
| ISG20   |                24 |      5.2917 |                 0.3918 | own member |           24 |
| SAMD9   |                24 |      6.2083 |                 0.3609 | own member |           24 |
| LY6E    |                24 |      6.4167 |                 0.3613 | own member |           24 |
| EIF2AK2 |                24 |      7.75   |                 0.3408 | own member |           24 |
| UBE2L6  |                24 |      9.875  |                 0.3167 | own member |           24 |
| PSMB9   |                24 |     10.7917 |                 0.304  | own member |           24 |
| B2M     |                24 |     11.6667 |                 0.2962 | own member |           24 |
| ISG15   |                24 |     11.75   |                 0.2962 | own member |           24 |
| PLSCR1  |                24 |     12.2083 |                 0.2916 | own member |           24 |
| SP110   |                24 |     13.625  |                 0.2769 | own member |           24 |
| IRF7    |                24 |     14.1667 |                 0.2749 | own member |           24 |
| IFI44   |                24 |     15.5417 |                 0.2583 | own member |           24 |
| PSMB8   |                24 |     15.9583 |                 0.261  | own member |           24 |
| TAP1    |                21 |     18.381  |                 0.2233 | own member |           24 |
| EPSTI1  |                16 |     19      |                 0.2255 | own member |           24 |
| NMI     |                16 |     19      |                 0.2204 | own member |           24 |

