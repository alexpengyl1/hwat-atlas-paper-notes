# 09 · Evidence, reproducibility, takeaways

[← Fig. 8](08-fig8-clinical.md) · [README](../README.md) · [Next: Glossary →](10-glossary.md)

**Short answer.** Nearly every downstream result passes through the sub-cluster marker tables, and the repository cannot regenerate several of the figures. The method itself is a short, reusable skeleton.

## Contents

- [What the paper is betting on](#what-the-paper-is-betting-on)
- [What the authors gave up](#what-the-authors-gave-up)
- [Which design choices carry the paper?](#which-design-choices-carry-the-paper)
- [How far does each claim reach?](#how-far-does-each-claim-reach)
- [What happens if you try to rerun the repository?](#what-happens-if-you-try-to-rerun-the-repository)
- [Where the scripts and Methods disagree](#where-the-scripts-and-methods-disagree)
- [Takeaways](#takeaways)
- [Still open](#still-open)

## What the paper is betting on

**The problem it reacts to.** Earlier single-cell studies of human fat each had fewer than 15 donors, used different platforms and pipelines, and named their cell types differently. The names did not line up between papers.

**The bet.** Clusters defined independently in each study share enough marker genes that the following order works: cluster each dataset, link clusters by marker overlap, split into broad classes, and only then integrate. Fig. 1a (the Jaccard column) and Fig. 1b (the network) carry this bet, which is why [chapter 01](01-fig1ab-jaccard-network.md) goes into such detail.

## What the authors gave up

Most of these are stated in the text. A few appear only as comments in the code.

| Given up | Source |
|---|---|
| The whole Hildreth dataset (83,546 cells) | Low marker overlap and failed scANVI integration. Main text and a comment in the network script |
| Choosing network edges by the number of shared genes | Script comment: `### NUMBER not good selection` |
| Integrating all cells together | The metrics are better after splitting into four classes. Main text |
| Analysing subcutaneous and omental FAPs together | Script comment in `FAP_MSL.R`: integration not good between sc and om, due to mesothelial-like cells |
| Defining adipocyte subtypes from snSeq | The whole of Fig. 5 argues that this cannot be done. Adipocytes are pooled into one class |

## Which design choices carry the paper?

| Choice | What breaks without it | Evidence for it |
|---|---|---|
| Cluster each dataset first and define four classes from marker overlap (Fig. 1b) | No basis for per-class integration, and no reason to drop Hildreth | Fig. 1a–b |
| Integrate each class separately | The metrics get worse | Fig. 1e, all-cells columns against per-class columns in the source data |
| scVI as the integration method | Unknown. rPCA scores better on ARI and LISI, and the sub-clusters were not compared between methods | Fig. 1e |
| Split FAPs by depot | Subcutaneous and omental FAPs do not mix | A script comment and one sentence in the text. No quantified control |
| Pool adipocytes and use spatial markers for their subtypes | Platform bias would be read as subtypes | Fig. 5 |
| BisqueRNA marker scores as the bridge from single cells to bulk | Fig. 2c, 2g, 3g and 8 all rest on it | No comparison with another approach, such as reference-based deconvolution |
| cell2location | Less critical. Six tools agree at the level of the five classes | Supplementary Fig. 7a |

Parts that are closer to decoration: the circle plots in Fig. 6b (the sub-clusters are hand-picked), the overlay in Fig. 7d (one section), and the colour schemes of the UMAPs.

## How far does each claim reach?

| Claim | Evidence | Limit |
|---|---|---|
| Depot enrichment (Fig. 2b, 2f, 3f) | Pooled cell proportions, supported by marker scores in two paired bulk cohorts | No per-person test in the single-cell data. Perivascular snSeq comes from one cohort of 3 men, which the authors point out |
| New sub-clusters (myC08, myC12, vC05, vC08, vC09) | Markers, presence across cohorts and methods, flow cytometry for some | No functional test. The authors write "may represent" and "it is possible that" |
| Two FAP routes | Pseudotime plus marker enrichment during in-vitro differentiation | Route 1 has matching in-vitro data. Route 2 is not detected in vitro, and the reading that culture cannot reproduce it is an interpretation |
| Cell–cell communication | Inference from expression | Secreted signalling only. No spatial or protein constraint |
| Spatial niches | Within-spot correlation and two immunofluorescence stains | Subcutaneous fat only, 10 people |
| Clinical links | Cross-sectional correlations and before-and-after data from two surgery cohorts | Association only. No other kind of intervention. A score is not a cell count |

**Not reported in the paper, within what was read.**

- The rule for choosing each clustering resolution. The text says only that several were tried.
- The distribution of cell numbers per person.
- A calibration of BisqueRNA scores against true cell proportions.
- The effect of age, BMI and similar variables on the single-cell atlas itself. The authors say explicitly that they did not analyse this.

## What happens if you try to rerun the repository?

The repository documents the analysis. It is not a pipeline that runs end to end.

**Missing inputs.** `Seurat_Merged.rds`, the integrated objects, the marker Excel files, `WAT_scRNA_MTX.h5ad` and the bulk expression sets are not included. They have to be rebuilt from the Mendeley deposit (doi:10.17632/y3pxvr4xbf.2) and GEO.

**Missing steps.**

- Per-dataset preprocessing
- ARI, LISI and kBET
- Capping of cell numbers per dataset
- CellTypist training (only the trained models are provided)
- Fig. 4g, Fig. 5 and Fig. 8
- Statistics and plots of Fig. 1c–d

**Script errors and slips.**

| File | Problem |
|---|---|
| `Deconvolution_Depots.R` | Two `write.table(...)` calls lack a closing bracket (L89, L238) |
| `rPCA.R` | Uses `red_use` and `dims_use`, which are never defined |
| `FAP_MSL.R` | The argument `n.cpvatponents` (L133) is a casualty of replacing "om" with "pvat" in `n.components`. monocle3 and SeuratWrappers are not loaded |
| `Immune_Cells.R`, `Vascular_Cells.R` | Totals are built with `aggregate(..., by = list(depot = ...))` and then read as `$tissue`. Run as written this gives missing values, yet the source data hold normal percentages |
| `Spatial_Deconvolution/work.sh` | Section 51 appears twice and section 48 is missing, but the downstream scripts read results for section 48 |
| `FAP_comparisons.R` | The subcutaneous-against-perivascular loop (L116–129) sits outside the outer loop and compares only the last cluster. Fig. 4j is not affected |
| `FAP_MSL.R` | A feature named `"LY6"` is requested. It is not a valid human gene symbol and Seurat skips it with a warning |
| `Jaccard_and_Network_Fig1.R` | The `geom_jitter` layer is not attached to the plot. Three datasets skip the marker filter loop |

The fourth row is the most telling. It shows that the published scripts differ from the versions that produced the figures.

## Where the scripts and Methods disagree

| Setting | Methods | Script |
|---|---|---|
| Minimum marker genes for BisqueRNA | 6 | 12 and 8 (depot comparison) |
| Ranking of network edges | Five highest Jaccard values | Five highest average overlap percentages |
| Marker threshold for Fig. 1 | log2FC > 0.5 | log2FC > 0 |
| CellChat database | CellChatDB | Secreted signalling subset only |
| Clusters in the co-localisation map | Three | Tree cut into four |

## Takeaways

1. **The skeleton transfers to other tissues.** Cluster each dataset on its own. Link clusters by marker overlap to define broad classes. Integrate each class with scVI. Re-cluster on the latent space. Treat the resulting marker tables as a common interface to bulk data (BisqueRNA), spatial data (cell2location) and communication analysis (CellChat).
2. **One interface carries almost everything.** Nearly every downstream conclusion passes through the sub-cluster marker tables. Their thresholds, and whether genes shared between sub-clusters are removed, decide how far Fig. 2c to Fig. 8 can be trusted.
3. **Fig. 5 is the most useful figure for method.** It shows how to use bulk and spatial data as outside references to judge whether heterogeneity seen in single-nucleus data is real.
4. **Read each "deconvolution" for what it computes.** The paper uses the word for three different things: a PCA score on marker genes (BisqueRNA), a probabilistic mixture model (cell2location), and a marker-set fold enrichment (Fig. 4g).
5. **A benchmark answers the question its metrics ask.** Fig. 1e measures batch mixing. It says nothing about whether cell types stay separate.

## Still open

- How many of the 60-plus sub-clusters survive a different integration method, rPCA in particular.
- The quantitative relation between a marker score and the true proportion of a cell type.
- What route 2 FAPs become in the tissue.
- Whether vC05 is a real transitional state between endothelial cells and preadipocytes.
- How much of the negative spatial correlation reflects real exclusion, and how much comes from proportions summing to 1.

[← Fig. 8](08-fig8-clinical.md) · [README](../README.md) · [Next: Glossary →](10-glossary.md)
