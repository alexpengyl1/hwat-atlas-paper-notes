# 10 · Glossary

[← Evidence and reproducibility](09-evidence-and-reproducibility.md) · [README](../README.md)

Abbreviations and names used in the paper, the scripts and these notes.

## Tissue and depots

| Term | Meaning |
|---|---|
| WAT | White adipose tissue |
| `sc` | Subcutaneous fat, under the skin. In "scSeq" the same two letters mean single-cell |
| `om` | Omental fat, a visceral depot inside the abdomen |
| `pv`, `pvat` | Perivascular fat, around blood vessels |
| Depot | An anatomical location of fat tissue |
| SVF | Stromal vascular fraction: the non-adipocyte cells obtained after digesting fat tissue and removing the floating fat cells |

## Sequencing methods and data types

| Term | Meaning |
|---|---|
| scSeq | Single-cell RNA sequencing of whole cells. It does not capture mature adipocytes |
| snSeq | Single-nucleus RNA sequencing. It captures adipocytes, but measures nuclear RNA only |
| STx | Spatial transcriptomics. Here, 10x Visium |
| Visium spot | One 55 µm capture area on a Visium slide. It holds the mixed RNA of the cells lying on it |
| Bulk RNA-seq | Sequencing of a whole tissue sample, with no single-cell resolution |
| Object | The paper's word for one cell, one nucleus or one spot |
| Pseudobulk | Single-cell data averaged over many cells to imitate a bulk sample |

## Cell classes and cell types

| Term | Meaning |
|---|---|
| FAP | Fibroblast and adipogenic progenitor. The stromal class that includes fibroblasts and cells that can become adipocytes |
| CPA | Committed preadipocyte |
| APC | Adipose precursor cell |
| MSL | Mesothelial-like cell |
| LAM | Lipid-associated macrophage. Markers include TREM2 and CD9 |
| MMe | Metabolic-regulated (metabolically activated) macrophage |
| Mox | Redox-regulatory metabolic macrophage |
| M1-like, M2-like | Broad labels for pro-inflammatory and tissue-remodelling macrophage states |
| DC | Dendritic cell |
| NK, NKT | Natural killer cell, natural killer T cell |
| Th1 | T helper 1 cell |
| TREG | Regulatory T cell |
| TRM | Tissue-resident memory T cell |
| EC | Endothelial cell |
| VSMC | Vascular smooth muscle cell |
| Areg | Adipogenesis-regulatory cell, an anti-adipogenic FAP population described in mice |

## Cluster names

| Prefix | Class | Range |
|---|---|---|
| `lyC` | Lymphoid (T, NK, NKT) | lyC0–10 |
| `myC` | Myeloid (monocytes, macrophages, dendritic cells) | myC0–15 |
| `bC` | B cells | bC0–01 |
| `vC` | Vascular | vC0–11 |
| `sfC` | Subcutaneous FAP | sfC0–16 |
| `ofC` | Omental FAP | ofC0–14 |
| `pfC` | Perivascular FAP | pfC0–07 |
| `A_3` and similar | Cluster 3 of dataset A, before integration (network script only) |

Numbers follow cluster size. Cluster 0 is the largest.

## Statistics and metrics

| Term | Meaning | Chapter |
|---|---|---|
| Marker gene | A gene expressed more in one cluster than in the rest | 01 |
| log2FC | Log2 fold change of expression, cluster against the rest | 01 |
| Adjusted p | P value corrected for multiple testing | 01 |
| Jaccard index | Shared items divided by all items of two sets. 0 means no overlap, 1 means identical | 01 |
| Batch | A source of technical variation. Here a dataset, split further by depot and method (`orig.ident`) | 02 |
| Integration | Removing batch differences so that cells from different sources can be analysed together | 02 |
| Latent space | A low-dimensional representation of each cell. scVI gives 10 dimensions here | 02 |
| ARI | Adjusted Rand index. Agreement between two groupings of the same cells | 02 |
| LISI | Local inverse Simpson index. The effective number of batches around a cell | 02 |
| kBET | k-nearest-neighbour batch effect test. Whether the local batch mix matches the overall mix | 02 |
| Resolution | The clustering parameter that sets how many clusters come out. Higher gives more | 03 |
| UMAP | A two-dimensional map of cells for display | 03 |
| Pseudotime | A number per cell: its distance along a learned trajectory from a chosen starting point | 04 |
| z-score | A value expressed in standard deviations from the mean | 05 |
| Communication probability | CellChat's score for one ligand–receptor pair between a sender and a receiver group | 06 |
| Deconvolution | Estimating cell-type content of a mixed sample. The paper uses the word for three different calculations | 07, 08 |
| q05 abundance | The 5% quantile of cell2location's estimated abundance. A cautious lower bound | 07 |
| Pearson correlation | Linear correlation. Used across spots in Fig. 7c | 07 |
| Spearman correlation | Rank correlation. Used between scores and clinical measures in Fig. 8 | 08 |
| Fisher's z | A transformation that lets correlation coefficients be averaged | 08 |
| Common-effect and random-effects models | Two ways of pooling results across cohorts. The second allows true effects to differ between cohorts | 08 |

## Tools

| Tool | Role in the paper |
|---|---|
| Seurat | Per-dataset processing, clustering, markers, plotting |
| SCTransform | Normalisation in Seurat |
| Harmony | Batch correction. Used per dataset for subject effects, and as one benchmark method |
| rPCA | Seurat's reciprocal PCA integration. Benchmark method |
| BBKNN | Batch-balanced k-nearest neighbours. Benchmark method |
| scVI | Variational autoencoder for single-cell data. The integration method chosen |
| scANVI | Semi-supervised version of scVI that uses cell-type labels. Used for the first trial integration |
| sceasy, reticulate | Conversion between Seurat and AnnData, and calling Python from R |
| CellTypist | Automated cell-type annotation. Used to check immune labels. The authors publish trained models |
| igraph, RCy3, Cytoscape | Building and drawing the marker-overlap network |
| monocle3 | Trajectory and pseudotime analysis |
| BisqueRNA | Marker-based scores of cell types in bulk data |
| CellChat | Inference of cell–cell communication from ligand and receptor expression |
| cell2location | Spatial deconvolution, the main tool for Fig. 7 |
| RCTD, Tangram, DestVI, stereoscope, SPOTlight | Other spatial deconvolution tools, used for comparison |
| Hmisc, meta | Spearman correlations, and meta-analysis across cohorts |
| pdfCluster, lisi, kBET | The three integration metrics |

## Clinical measures

| Term | Meaning |
|---|---|
| BMI | Body mass index |
| WHR | Waist-to-hip ratio |
| HOMA-IR | Homeostatic model assessment of insulin resistance, calculated from fasting glucose and insulin |
| HDL | High-density lipoprotein cholesterol |
| TG | Triglycerides |
| Leptin | A hormone secreted by fat cells, higher with more fat mass |
| Fat cell volume | Average size of adipocytes |
| Lipolysis | Release of fatty acids and glycerol from fat cells. Measured at rest (basal) and after isoprenaline, a drug that stimulates it |
| Bariatric surgery | Weight-loss surgery |

## File types

| File | Content |
|---|---|
| `.rds` | A saved R object, here usually a Seurat object |
| `.h5ad` | An AnnData object, the Python format for single-cell data |
| `.cys` | A Cytoscape session. It is a zip archive with node tables, edge tables and visual styles |
| `.pkl` | A trained CellTypist model |
| `Source Data/` | The numbers behind each figure panel, as provided by the authors |

[← Evidence and reproducibility](09-evidence-and-reproducibility.md) · [README](../README.md)
