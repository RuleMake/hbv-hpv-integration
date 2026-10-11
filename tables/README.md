# Supplemental numerical tables

| Table | Contents | Rows |
|---|---|---:|
| S1 | Class enrichment at ±100 kb: 24 paired-call tests and 24 donor-signature sensitivity tests | 48 |
| S2 | Local composition displacements and HBV–HPV contrasts under genomic and breakpoint-density matching | 24 |
| S3 | Composition under consensus, HeLa-S3 and HepG2 replication-timing references | 18 |
| S4 | Replication-timing class atlas: 48 paired-call tests and 48 donor-signature sensitivity tests | 96 |
| S5 | Donor-average duplication proportions by detected integration status, including histology-standardized HPV analysis | 3 |
| S6 | Duplication enrichment near non-hotspot HPV sites under gene-status and gene-length matching | 24 |
| S7 | Descriptive summary of the supporting detection-threshold registry | 11 |

`paired_sample` retains sample-level SV calls with both breakends. `donor_signature` groups identical chromosome, position and orientation pairs within each donor. Paired calls are the primary representation in S1–S4; donor signatures are primary in S6. Both representations use the same underlying data.

Call-count estimators sum distinct calls within each cluster; a call can contribute to more than one cluster. The calibrated cluster-sign statistic counts clusters whose observed counts exceed their own spatial-null means. Its randomized values use leave-one-out null means. The S4/S6 field labels `event_incidence` and `cluster_event_incidence` refer to these call counts.

Composition scores use (TRA + INV) − (DEL + DUP), where each term is a local breakend fraction. The HBV-minus-HPV contrast compares the two observed-minus-null scores. S2 and S3 retain unscaled scores; multiplying by 100 gives the score scale used in the figures. The corresponding difference in the TRA + INV fraction is half the score difference. In S2, `stratum` denotes genomic-covariate matching and `density` adds local breakpoint-density matching. `full` uses the complete eligible set; `strict` uses the subset requiring no density-pool fallback.

Spatial reference ranges in S1–S4 and S6 are derived from null distributions. S5 reports 95% donor-bootstrap confidence intervals. Detectable-effect fields describe approximate sensitivity under the stated power assumptions. All numerical strings retain the precision of the saved analysis outputs.

Use the two-sided adjusted P columns for atlas inference: `emp_p_two_sided_bonferroni_24` in S1 and `p_bonferroni_48` in S4. S3 uses `p_bonferroni_3` for the three primary paired-call HBV-minus-HPV contrasts. S5 uses `p_holm_two_pooled_tests` for the two pooled donor comparisons. S6 uses `p_two_sided_bonferroni2` for the two prespecified protein-coding, gene-length-matched donor-signature estimators. Blank adjusted P fields identify rows outside the corresponding test family. Spatial empirical P values use inclusive upper and lower tails with add-one correction; the two-sided value is twice the smaller tail, capped at 1. S5 uses the absolute permutation statistic.

## S1. Structural-variant class enrichment near viral integration sites

The table reports call-count and calibrated cluster-sign enrichment within ±100 kb of HBV, HPV and EBV integration clusters. It contains 24 tests for paired sample-level calls and 24 sensitivity tests for within-donor breakpoint-pair signatures. A call is counted once within a cluster and may contribute to several clusters. The two-sided P values are adjusted across 24 tests within each representation; fold ranges describe spatial randomization.

| Column | Definition |
|---|---|
| `variant` | SV representation: paired sample-level calls or within-donor breakpoint-pair signatures. |
| `svclass` | Structural-variant class: deletion (DEL), duplication (DUP), inversion (INV) or translocation (TRA). |
| `n_events` | Number of reference SV calls in this class and representation. |
| `virus` | Viral catalogue, or the stated between-virus contrast. |
| `estimator` | bc_uniqSV_blockNull: summed within-cluster call counts; f_clusterSignTest: calibrated number of clusters with an excess. |
| `n_coords` | Number of retained integration coordinates; blank for a between-virus contrast. |
| `n_clusters` | Number of retained integration clusters; blank for a between-virus contrast. |
| `seed` | Random-number seed. |
| `window_bp` | Query-window half-width in base pairs. |
| `observed` | Summed within-cluster call counts, or number of clusters whose counts exceed their own null means. |
| `expected` | Mean of the randomized statistic. |
| `null_sd` | Standard deviation of the randomized statistic. |
| `fold` | Observed statistic divided by its spatial-null mean. |
| `fold_null_interval_low` | Observed statistic divided by the 97.5th null percentile. |
| `fold_null_interval_high` | Observed statistic divided by the 2.5th null percentile. |
| `null_quantile_025` | 2.5th percentile of the randomized statistic. |
| `null_quantile_975` | 97.5th percentile of the randomized statistic. |
| `z` | Observed-minus-null-mean difference divided by the null standard deviation. |
| `emp_p` | Upper-tail empirical P value with add-one correction. |
| `emp_p_lower` | Lower-tail empirical P value with add-one correction. |
| `emp_p_two_sided_equal_tail` | Twice the smaller empirical tail P value, capped at 1. |
| `n_perm` | Number of spatial randomizations. |
| `emp_p_bonferroni_24` | Upper-tail empirical P value multiplied by 24 and capped at 1. |
| `emp_p_two_sided_bonferroni_24` | Two-sided empirical P value multiplied by 24 and capped at 1; use for atlas conclusions. |

## S2. Local rearrangement composition under genomic and breakpoint-density matching

The table compares local composition around integration clusters with matched spatial references, using SVs from donors without detected integration. The score is (TRA + INV) − (DEL + DUP), averaged equally across clusters over ±1 Mb. Rows report individual viral catalogues and the HBV-minus-HPV difference in observed-minus-null scores. Genomic-stratum matching is compared with added breakpoint-density matching, including the strict subset requiring no density-pool fallback. Scores and contrasts are unscaled.

| Column | Definition |
|---|---|
| `variant` | SV representation: paired sample-level calls or within-donor breakpoint-pair signatures. |
| `virus` | Viral catalogue, or the stated between-virus contrast. |
| `arm` | stratum: genomic-covariate matching; density: additional local breakpoint-density matching. |
| `window_bp` | Query-window half-width in base pairs. |
| `n_coords` | Number of retained integration coordinates; blank for a between-virus contrast. |
| `n_clusters` | Number of retained integration clusters; blank for a between-virus contrast. |
| `n_perm` | Number of spatial randomizations. |
| `null_fallback_pct` | Percentage of candidate-placement draws that used a fallback pool; blank for between-virus contrasts. |
| `rung` | full: complete eligible set; strict: subset requiring no density-pool fallback. |
| `statistic` | Local composition score or HBV-minus-HPV contrast of observed-minus-null scores. |
| `observed` | Observed score; contrast rows give the HBV-minus-HPV difference in observed-minus-null scores. |
| `null_mean` | Mean of the randomized statistic. |
| `null_reference_q025` | 2.5th percentile of the randomized statistic. |
| `null_reference_q975` | 97.5th percentile of the randomized statistic. |
| `interval_type` | Label identifying a spatial-randomization reference range. |
| `excess` | Observed statistic minus its spatial-null mean. |
| `fold` | Unused for these signed statistics; values are blank. |
| `null_sd` | Standard deviation of the randomized statistic. |
| `z` | Observed-minus-null-mean difference divided by the null standard deviation. |
| `n_null_ge_observed` | Number of randomized values greater than or equal to the observed statistic. |
| `n_null_le_observed` | Number of randomized values less than or equal to the observed statistic. |
| `emp_p_greater` | Upper-tail empirical P value with add-one correction. |
| `emp_p_less` | Lower-tail empirical P value with add-one correction. |
| `emp_p_two_sided` | Twice the smaller empirical tail P value, capped at 1. |
| `mde_excess_approx_2p80sd` | Approximate detectable score difference: 2.80 times the null standard deviation. |
| `mde_fold_approx` | Unused for these signed statistics; values are blank. |
| `mde_interpretation` | Label identifying approximate sensitivity rather than an effect bound. |
| `null_stream_corr` | Pearson correlation of the HBV and HPV centred null streams; populated for contrast rows. |

## S3. Replication-timing sensitivity of local rearrangement composition

The table reports composition analyses using the consensus, HeLa-S3 and HepG2 replication-timing references on a common set of eligible clusters. Each reference is applied to both viruses. The 18 rows comprise two SV representations, three references and three contrasts: HBV, HPV and HBV minus HPV. Bonferroni adjustment covers the three paired-call HBV-minus-HPV contrasts. Excess-reference limits are the observed statistic minus the opposite spatial-null quantile.

| Column | Definition |
|---|---|
| `variant` | SV representation: paired sample-level calls or within-donor breakpoint-pair signatures. |
| `track` | Replication-timing reference: 15-cell-line consensus, HeLa-S3 or HepG2. |
| `contrast` | HBV or HPV composition score, or HBV-minus-HPV contrast of observed-minus-null scores. |
| `n_clusters` | Number of retained integration clusters; blank for a between-virus contrast. |
| `n_coordinates` | Number of retained integration coordinates; blank for a between-virus contrast. |
| `observed` | Observed score; HBV-minus-HPV rows give the difference in observed-minus-null scores. |
| `null_mean` | Mean of the randomized statistic. |
| `excess` | Observed statistic minus its spatial-null mean. |
| `null_sd` | Standard deviation of the randomized statistic. |
| `null_reference_q025` | 2.5th percentile of the randomized statistic. |
| `null_reference_q975` | 97.5th percentile of the randomized statistic. |
| `excess_reference_low` | Observed statistic minus the 97.5th null percentile. |
| `excess_reference_high` | Observed statistic minus the 2.5th null percentile. |
| `p_upper` | Upper-tail empirical P value with add-one correction. |
| `p_lower` | Lower-tail empirical P value with add-one correction. |
| `p_two_sided` | Twice the smaller empirical tail P value, capped at 1. |
| `n_permutations` | Number of spatial randomizations. |
| `interval_type` | Label identifying a spatial-randomization reference range. |
| `p_bonferroni_3` | Two-sided P value adjusted across the three paired-call HBV-minus-HPV contrasts; blank for other rows. |

## S4. Replication-timing sensitivity of the structural-variant class atlas

The table reports call-count and calibrated cluster-sign enrichment within ±100 kb under three replication-timing references, using clusters eligible under every reference. Its 96 rows contain 48 paired-call tests and 48 donor-signature sensitivity tests across two viruses, four SV classes and two estimators. The paired-call P values are adjusted across 48 tests. Fold ranges are derived from spatial-null quantiles.

| Column | Definition |
|---|---|
| `variant` | SV representation: paired sample-level calls or within-donor breakpoint-pair signatures. |
| `virus` | Viral catalogue, or the stated between-virus contrast. |
| `track` | Replication-timing reference: 15-cell-line consensus, HeLa-S3 or HepG2. |
| `svclass` | Structural-variant class: deletion (DEL), duplication (DUP), inversion (INV) or translocation (TRA). |
| `estimator` | event_incidence: summed within-cluster call counts; calibrated_sign: calibrated number of clusters with an excess. |
| `n_clusters` | Number of retained integration clusters; blank for a between-virus contrast. |
| `n_coordinates` | Number of retained integration coordinates; blank for a between-virus contrast. |
| `observed` | Summed within-cluster call counts, or number of clusters whose counts exceed their own null means. |
| `expected` | Mean of the randomized statistic. |
| `fold` | Observed statistic divided by its spatial-null mean. |
| `null_quantile_025` | 2.5th percentile of the randomized statistic. |
| `null_quantile_975` | 97.5th percentile of the randomized statistic. |
| `fold_null_low` | Observed statistic divided by the 97.5th null percentile. |
| `fold_null_high` | Observed statistic divided by the 2.5th null percentile. |
| `p_upper` | Upper-tail empirical P value with add-one correction. |
| `p_lower` | Lower-tail empirical P value with add-one correction. |
| `p_two_sided` | Twice the smaller empirical tail P value, capped at 1. |
| `n_permutations` | Number of spatial randomizations. |
| `interval_type` | Label identifying a spatial-randomization reference range. |
| `p_bonferroni_48` | Two-sided P value adjusted across 48 paired-call tests; blank for donor-signature rows. |

## S5. Donor-level duplication proportions by detected integration status

The table compares donor means of the fraction of large duplications with a breakend within ±100 kb of the shared viral integration catalogue. Donors are grouped by detected integration status within the specified histologies and contribute equally to each pooled mean. The HPV common-support analysis standardizes stratum-specific donor means to fixed histology weights. Results include donor-bootstrap confidence intervals and studentized label-permutation P values; Holm adjustment covers the two pooled comparisons.

| Column | Definition |
|---|---|
| `analysis` | Named pooled or histology-standardized comparison. |
| `virus` | Viral catalogue, or the stated between-virus contrast. |
| `coordinate_set` | Shared, jointly eligible integration coordinates after amplification-peak masking. |
| `n_coordinates` | Number of coordinates in the shared viral catalogue. |
| `window_bp` | Query-window half-width in base pairs. |
| `histologies` | Included histologies, in the order used by stratum_weights. |
| `n_detected` | Number of donors with detected integration in the comparison. |
| `n_comparison` | Number of comparison donors without detected integration. |
| `mean_detected` | Mean donor duplication proportion in the detected group, using the specified stratum weights. |
| `mean_comparison` | Mean donor duplication proportion in the comparison group, using the specified stratum weights. |
| `mean_difference` | Detected-group mean minus comparison-group mean. |
| `difference_percentage_points` | Mean difference multiplied by 100, in percentage points. |
| `diff_ci_lo` | Lower limit of the 95% donor-bootstrap percentile interval for the mean difference. |
| `diff_ci_hi` | Upper limit of the 95% donor-bootstrap percentile interval for the mean difference. |
| `mean_ratio` | Detected-group mean divided by comparison-group mean. |
| `ratio_ci_lo` | Lower limit of the 95% donor-bootstrap percentile interval for the mean ratio. |
| `ratio_ci_hi` | Upper limit of the 95% donor-bootstrap percentile interval for the mean ratio. |
| `welch_statistic` | Mean difference divided by its estimated standard error, with the specified stratum weights. |
| `p_two_sided_studentized_permutation` | Two-sided label-permutation P value for the studentized mean difference. |
| `p_two_sided_raw_difference_permutation` | Two-sided label-permutation P value for the unstudentized mean difference. |
| `p_holm_two_pooled_tests` | Holm-adjusted studentized P value for the two pooled comparisons; blank for the standardized analysis. |
| `n_bootstrap` | Number of donor-bootstrap replicates. |
| `n_permutations` | Number of label permutations. |
| `bootstrap_zero_comparison_means` | Number of bootstrap replicates with a zero comparison-group mean. |
| `seed_base` | Base random-number seed used to derive analysis-specific streams. |
| `stratum_weights` | Fixed histology weights in histologies order; 1 denotes the single pooled group. |

## S6. Gene-matching sensitivity of duplication enrichment near non-hotspot HPV sites

The table compares duplication enrichment near non-hotspot HPV coordinates under genomic matching, gene-overlap-status matching and additional gene-length matching. Within each annotation scope, all three arms use the same eligible coordinates and clusters. Donor signatures are the primary SV representation; paired sample-level calls provide sensitivity results. Bonferroni adjustment covers the two estimators in the protein-coding, gene-length-matched donor-signature analysis. Fold ranges describe spatial randomization.

| Column | Definition |
|---|---|
| `scope` | Gene annotation scope: protein-coding genes or all annotated genes. |
| `variant` | SV representation: paired sample-level calls or within-donor breakpoint-pair signatures. |
| `arm` | joint_geometry: genomic matching; gene_status: added overlap-status matching; gene_length: added gene-length matching. |
| `estimator` | cluster_event_incidence: summed within-cluster call counts; calibrated_cluster_sign: calibrated number of clusters with an excess. |
| `n_coordinates` | Number of retained integration coordinates; blank for a between-virus contrast. |
| `n_clusters` | Number of retained integration clusters; blank for a between-virus contrast. |
| `observed` | Summed within-cluster call counts, or number of clusters whose counts exceed their own null means. |
| `expected` | Mean of the randomized statistic. |
| `fold` | Observed statistic divided by its spatial-null mean. |
| `null_quantile_025` | 2.5th percentile of the randomized statistic. |
| `null_quantile_975` | 97.5th percentile of the randomized statistic. |
| `fold_null_low` | Observed statistic divided by the 97.5th null percentile. |
| `fold_null_high` | Observed statistic divided by the 2.5th null percentile. |
| `p_upper` | Upper-tail empirical P value with add-one correction. |
| `p_lower` | Lower-tail empirical P value with add-one correction. |
| `p_two_sided` | Twice the smaller empirical tail P value, capped at 1. |
| `n_permutations` | Number of spatial randomizations. |
| `interval_type` | Label identifying a spatial-randomization reference range. |
| `p_two_sided_bonferroni2` | Two-sided P value adjusted across the two primary protein-coding, gene-length-matched donor-signature estimators; otherwise blank. |

## S7. Detection-threshold registry

S7 summarizes saved registry entries rather than independent observations. Group quantiles combine several windows and estimators. The power-constant 10th–90th percentile range is 2.788162–2.809058; 2.640000–2.874791 is the minimum–maximum range. Field-coverage percentages overlap.
