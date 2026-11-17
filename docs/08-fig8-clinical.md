# 08 · Fig. 8: clinical associations

[← Fig. 7](07-fig7-spatial.md) · [README](../README.md) · [Next: Evidence and reproducibility →](09-evidence-and-reproducibility.md)

**Short answer.** In bulk fat samples from 864 people, sub-cluster scores fall into three groups by how they correlate with metabolic measures, and most of them move back towards normal after weight loss. Everything here is association.

The code for Fig. 8 is not in the repository. The description follows Methods, the figure legend, the source data, and the script that uses the same function for Fig. 2c (`scripts/Bulk_deconvolution/Deconvolution_Depots.R`).

## Contents

- [Why this figure exists](#why-this-figure-exists)
- [How is the bulk deconvolution done?](#how-is-the-bulk-deconvolution-done)
- [What the score is and is not](#what-the-score-is-and-is-not)
- [From scores to the figure](#from-scores-to-the-figure)
- [How are correlations pooled across cohorts?](#how-are-correlations-pooled-across-cohorts)
- [What is concluded](#what-is-concluded)
- [Where the evidence stops](#where-the-evidence-stops)

## Why this figure exists

Single-cell studies have few donors, so the clinical relevance of a newly named cell type is hard to judge. Bulk transcriptomes of fat exist for hundreds of people with clinical measurements. The idea is to carry the sub-cluster definitions over to those bulk data through their marker genes.

## How is the bulk deconvolution done?

It uses the **marker mode of BisqueRNA**. This mode does not estimate cell proportions.

| Input | Content |
|---|---|
| Marker table | For every integrated sub-cluster, the genes with `avg_log2FC > 0` and adjusted p < 0.05 from the marker Excel files of Fig. 2–4, as three columns: sub-cluster, gene, log2FC |
| Bulk expression | Eight cohorts of subcutaneous fat, 864 people in total. Each cohort is processed on its own |

```r
BisqueRNA::MarkerBasedDecomposition(eset, markers, unique_markers = FALSE,
                                    weighted = TRUE, min_gene = 6, w_col = "avg_log2FC")
```

| Argument | Meaning |
|---|---|
| `eset` | The bulk expression set of one cohort |
| `markers` | The three-column marker table |
| `unique_markers = FALSE` | A gene may be a marker of several sub-clusters |
| `weighted = TRUE`, `w_col = "avg_log2FC"` | Genes are weighted by their log2 fold change |
| `min_gene` | Sub-clusters with fewer markers than this are skipped. Methods says 6. The depot script uses 12 and 8 |

For each sub-cluster, the function does three things:

1. It takes the expression of that sub-cluster's marker genes across all bulk samples of the cohort.
2. It runs a PCA on that sub-matrix, with the genes weighted.
3. It returns the first principal component (PC1) as one score per sample.

**Illustrative example.** Three LAM markers, with expression standardised across samples. The real calculation is a PCA. Thinking of it as an average is a simplification that works when the genes move together.

| | TREM2 | CD9 | LPL | Score (approx.) |
|---|---|---|---|---|
| Sample 1 | +1.0 | +0.8 | +1.2 | +1.0 |
| Sample 2 | −1.0 | −0.8 | −1.2 | −1.0 |

The three genes rise and fall together, so PC1 is close to their average. Sample 1 has more of this marker set and is inferred to have relatively more LAMs.

## What the score is and is not

- **Signed and without unit.** It is a relative value within a cohort, not a percentage. The source data of Fig. 2 show this: for one person, the adipocyte score is −3.43 in the omental sample and +3.41 in the subcutaneous sample.
- **Comparable only within one cohort and one sub-cluster.** It can say that sample 1 has more LAM signal than sample 2. It cannot say that LAMs are more abundant than M2 macrophages, and scores from different cohorts cannot be put side by side.
- **More cells or more expression?** A high score means the marker set is expressed more in the tissue. That can come from more cells of that type, or from each cell expressing the markers more strongly. The score cannot separate the two.
- **Related sub-clusters move together.** With `unique_markers = FALSE`, one gene can serve several sub-clusters. The M2-like macrophage clusters share many markers, so their scores are correlated by construction.

## From scores to the figure

**Fig. 8a: the cohorts.** Eight cohorts with 864 people: Imbert 332, Armenise 189, Arner P#2 80, Arner P#1 78, Kerr 77, Arner E 56, Petrus 29, Krieg 23 (from `Figure_8a.txt`). The panel shows the spread of age, BMI and HOMA-IR.

**Fig. 8b: the heat map.**

1. Within each cohort, compute a score for every sub-cluster.
2. Within each cohort, compute the Spearman correlation between every sub-cluster score and every clinical measure (`Hmisc`).
3. Pool the correlations across cohorts by meta-analysis (`meta`), with both a common-effect and a random-effects model.
4. Arrange the pooled correlations as a matrix. `Figure_8b.txt` has 58 rows and 11 columns.
5. Cluster the matrix. This gives groups A, B and C.

| Dimension | Content |
|---|---|
| 58 sub-clusters | 2 B cell, 12 vascular, 16 myeloid, 17 subcutaneous FAP, 11 T and NK |
| 11 measures | Age, BMI, waist-to-hip ratio, HOMA-IR, HDL cholesterol, triglycerides, leptin, fat cell volume, basal lipolysis, isoprenaline-stimulated lipolysis, stimulated-to-basal ratio |

**Fig. 8c: forest plots** (`meta::forest`). For one sub-cluster and one measure, the correlation and 95% confidence interval of each cohort, plus the two pooled estimates. All examples shown have p < 0.0001.

**Fig. 8d–e: weight loss.** Two cohorts have bulk data before and after bariatric surgery: Petrus (15 people, before and 2 years after) and Kerr (37 people, before, 2 years and 5 years after). In `Figure_8d_8e.xlsx` the sample suffixes B, C and D stand for baseline, 2 years and 5 years. Fig. 8d shows whether groups A and B shift as a whole. Fig. 8e shows single sub-clusters as paired box plots with a two-sided paired t test.

## How are correlations pooled across cohorts?

`meta::metacor` uses Fisher's z transformation by default.

$$z_i=\mathrm{atanh}(r_i),\qquad w_i=n_i-3,\qquad \bar z=\frac{\sum_i w_i z_i}{\sum_i w_i},\qquad \bar r=\tanh(\bar z)$$

$r_i$ and $n_i$ are the correlation and sample size of cohort $i$.

**Worked example, common-effect model.**

| | r | n | z = atanh(r) | Weight = n − 3 |
|---|---|---|---|---|
| Cohort 1 | 0.5 | 103 | 0.549 | 100 |
| Cohort 2 | 0.3 | 53 | 0.310 | 50 |

Weighted mean: (0.549 × 100 + 0.310 × 50) / 150 = 0.469. Back-transformed: tanh(0.469) ≈ 0.44. The larger cohort has more weight, so the pooled value lies closer to 0.5. A random-effects model adds the variation between cohorts to the weights, which makes them more equal.

**Why correlate within each cohort first.** The scores are relative within a cohort, and the cohorts used different platforms (arrays and RNA-seq). Raw scores cannot be stacked across cohorts. Pooling at the level of correlation coefficients is the approach that fits this kind of score.

## What is concluded

| Group | Sub-clusters | Relation to metabolic state | Main members |
|---|---|---|---|
| A | 6 | Higher score with a metabolically healthier profile | Capillary endothelial cells (vC0) and intermediate or late FAP states (sfC04, sfC14, sfC16) |
| B | 15 | Higher score with higher BMI and insulin resistance | Mostly immune cells: LAMs (myC02), MMe (myC10), class 2 dendritic cells (myC03), six of the nine M2-like macrophage clusters (myC0, 01, 07, 08, 11, 13), early CD8 T cells (lyC02), CD16+ NK cells (lyC06). Among FAPs, only the CD55/PI16 progenitors (sfC02) |
| C | 37 | Weak associations | The rest |

How the authors read the groups:

- **Group B** fits earlier findings that certain immune cells are enriched in obesity and insulin resistance.
- **Capillary endothelial cells in group A** may reflect capillary rarefaction, the reduced capillary density reported in fat tissue of people with obesity.
- **FAPs on both sides.** Progenitors (sfC02) are in group B and some intermediate or late states are in group A. The authors suggest that "adipogenesis may be impacted by metabolic health at several different levels".
- **Age** is not associated with any sub-cluster.

Further results:

- **Adipocyte subtypes** (Supplementary Fig. 8a), scored with spatial markers because of [Fig. 5](05-fig5-adipocytes.md). The PLIN subtype falls with BMI, insulin resistance and leptin. The LEP subtype rises with all of them. The SAA subtype shows weak links.
- **Weight loss.** In both surgery cohorts (52 people), the scores of groups A and B "normalise" after weight loss. Capillary endothelial cells (vC0) are the exception and stay unchanged.
- **Weight regain.** When people are split into thirds by weight regain between years 2 and 5, several sub-clusters follow body weight, again except vC0 (Supplementary Fig. 8c).

The paper's summary: the cellular landscape of fat tissue is dynamic and mirrors changes in fat mass.

## Where the evidence stops

- **Association only.** The data are cross-sectional correlations plus before-and-after changes in two cohorts. The paper uses the word "association" throughout.
- **A score is not a cell count.** See above. The paper does not calibrate the scores against true proportions, for example with samples that also have flow cytometry or histology.
- **Fifteen is not fifteen.** The members of group B are not 15 independent findings. Six M2-like clusters plus LAM, MMe and dendritic cells share markers, and if Fig. 8 used `unique_markers = FALSE` as the depot script does, their scores move together.
- **Mechanisms are suggestions.** Capillary rarefaction and effects on adipogenesis are offered as explanations. Neither was tested.
- **vC0 after weight loss.** Why capillary endothelial cells do not recover is left unexplained.
- **Settings cannot be checked.** Methods states a minimum of six marker genes, and the depot script uses 12 and 8. Whether the heat map shows the common-effect or random-effects estimate, and which distance and linkage the clustering used, are not stated.
- **Incomplete columns.** Not every cohort measured all 11 variables. Lipolysis and fat cell volume exist in some cohorts only, so some cells of the heat map rest on fewer cohorts.
- **One kind of intervention.** Only bariatric surgery was studied. The authors list this as a limitation.

[← Fig. 7](07-fig7-spatial.md) · [README](../README.md) · [Next: Evidence and reproducibility →](09-evidence-and-reproducibility.md)
