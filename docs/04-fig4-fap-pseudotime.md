# 04 · Fig. 4: FAPs and pseudotime

[← Fig. 2–3](03-fig2-3-immune-vascular.md) · [README](../README.md) · [Next: Fig. 5 →](05-fig5-adipocytes.md)

**Short answer.** Subcutaneous FAPs branch from one CD55/PI16 progenitor cluster into two routes. The experiments support route 1 as the adipogenic route. Route 2 is not confirmed.

FAP stands for fibroblast and adipogenic progenitor. Scripts: `scripts/Cell_Types/FAP_MSL.R`, `FAP_compare_with_Dong_etal.R`, `FAP_comparisons.R`.

## Contents

- [How were the FAP sub-clusters defined?](#how-were-the-fap-sub-clusters-defined)
- [How was pseudotime computed?](#how-was-pseudotime-computed)
- [What does the pseudotime analysis conclude?](#what-does-the-pseudotime-analysis-conclude)
- [Did the experiments confirm it?](#did-the-experiments-confirm-it)
- [Fig. 4g: how are the in-vitro time courses scored?](#fig-4g-how-are-the-in-vitro-time-courses-scored)
- [Fig. 4h: do human FAPs match mouse FAPs?](#fig-4h-do-human-faps-match-mouse-faps)
- [Fig. 4i–j: are subcutaneous and omental FAPs the same cells?](#fig-4ij-are-subcutaneous-and-omental-faps-the-same-cells)
- [Fig. 4k: do the sub-clusters exist at the protein level?](#fig-4k-do-the-sub-clusters-exist-at-the-protein-level)
- [Where the evidence stops](#where-the-evidence-stops)

## How were the FAP sub-clusters defined?

(Fig. 4a, 4i, `FAP_MSL.R`)

All FAPs were clustered together first, at resolution 0.4. Subcutaneous and omental cells did not mix well. A script comment (L26) attributes this to mesothelial-like cells, which are abundant in omental fat. The FAPs were therefore split by the `tissue` column, and each depot was re-clustered on the same scVI latent space.

| Depot | Resolution | Clusters |
|---|---|---|
| Subcutaneous | 0.6 | 17 (sfC0–16) |
| Omental | 0.5 | 15 (ofC0–14) |
| Perivascular | 0.5 | 8 (pfC0–07) |

The subcutaneous clusters form one very large group (86.5% of all subcutaneous FAPs) and four small ones:

| Cluster | Identity |
|---|---|
| sfC02 | CD55/PI16-expressing progenitors, inside the large group |
| sfC0 | APOD/CFD-expressing cells, read as committed preadipocytes, inside the large group |
| sfC03 | Enriched for cell adhesion molecules. The only FAP cluster without CD34, PDGFRA or PDGFRB |
| sfC11 | Mesothelial-like cells |
| sfC13 | CD74+ stromal cells |
| sfC16 | Late committed preadipocytes |

The colours in Fig. 4a are not arbitrary. A comment at L232 says they follow four categories: stem-like, transient, committed, and not regulated in vitro. That classification comes from the time-course data in Fig. 4g.

Fig. 4b–c are marker UMAPs (`FeaturePlot`) for PDGFRA, PDGFRB, CD34, CD74, EZR, CD55 and PI16.

## How was pseudotime computed?

(Fig. 4d, `FAP_MSL.R` L254–275)

```r
subc.cds <- as.cell_data_set(subc)                  # SeuratWrappers
subc.cds <- cluster_cells(subc.cds)
subc.cds <- learn_graph(subc.cds, use_partition = TRUE)
subc.cds <- order_cells(subc.cds, reduction_method = "UMAP")
```

- **Where the trajectory is learned.** monocle3 learns a graph through the cells on the UMAP coordinates, which themselves come from the scVI latent space.
- **Where it starts.** `order_cells` is called without a root. In that form monocle3 opens an interactive window and the user clicks the starting node. The authors picked sfC02, the CD55/PI16 cluster, because the literature describes such cells as universal fibroblast progenitors that can become more specialised cells, adipocytes included.
- **What pseudotime is.** Each cell gets a number: its distance along the graph from the root. The value is written back to the Seurat object and coloured with `FeaturePlot(subc, "Pseudotime")`.
- **Where the routes come from.** Route 1 and route 2 are read off the branches of the graph. No code extracts them into variables.

## What does the pseudotime analysis conclude?

From sfC02 the graph splits at the first step into two branches that run in opposite directions.

| Route | Clusters along it | End point |
|---|---|---|
| Route 1 | sfC04, sfC05, sfC07, sfC08 (transient states) | sfC0, high in APOD and CFD, "most likely" committed preadipocytes |
| Route 2 | sfC01, sfC06, sfC09, sfC10, sfC12, sfC14, sfC15 | No identity given |

In one sentence: the CD55/PI16 progenitors have two fates, one leading to preadipocytes and one whose destination is unknown.

1. **sfC02 is the starting point.** This is assumed from the literature. Pseudotime did not derive it.
2. **Route 1 is the adipogenic route.** In-vitro adipogenesis shows the same marker sets changing in the same order.
3. **Route 2 is a separate branch.** It is not seen during in-vitro adipogenesis. The authors conclude that standard culture reproduces route 1 only.

Point 3 is the message of this part of the paper: FAPs in human subcutaneous fat are more diverse than culture models and mouse studies show, and a whole branch is outside what those models cover.

## Did the experiments confirm it?

| Panel | What was done | Result | Supports |
|---|---|---|---|
| 4e–g | CD55+ human adipose-derived stem cells were differentiated in vitro. Bulk RNA-seq at each time point was scored with the markers of each FAP cluster | sfC02 markers are high before induction and fall afterwards. Markers of the route 1 transient clusters rise briefly. sfC0 markers rise late | Route 1 |
| Supp. 4b | The same in three other cell models | Same trend | Route 1 |
| 4h | Marker overlap with FAP clusters from mouse inguinal fat | Mouse early progenitors match sfC02, committed preadipocytes match sfC0, the Areg population matches sfC04 | Route 1 |
| 4k | Flow cytometry for CD55, APOD, CD74 and EZR | The cell types exist at the protein level | Existence only |

**Route 2.**

- None of its clusters is detected during in-vitro adipogenesis.
- None overlaps with the mouse FAP clusters.
- The authors suggest that standard culture conditions cannot reproduce the full diversity. That reading is not tested.
- One hint points the other way. DPP4, an accepted marker of pro-adipogenic FAPs in mouse inguinal fat, falls on route 2 (sfC06 and sfC10). Whether route 2 cells can form fat is left open.

## Fig. 4g: how are the in-vitro time courses scored?

There is no code for this panel. Methods says that for each FAP cluster and each time point, the expression of its top 30 marker genes is computed relative to the background of that bulk dataset. Four systems were used: adipose-derived stem cells and primary stromal cells (both from FANTOM5), SGBS cells and hMADS cells.

My reading, since the paper gives no formula:

$$\text{score}_{c,t}=\frac{\text{mean expression of the top 30 markers of cluster } c \text{ at time } t}{\text{background expression of the dataset}}$$

**Illustrative example.** The 30 markers of sfC02 average 120 on day 0 against a background of 40, a score of 3.0. On day 8 they average 30, a score of 0.75. Read as: the CD55/PI16 gene set is switched off after induction.

The legend calls this "deconvolution" too. It is a simple marker-set enrichment and not the BisqueRNA score used elsewhere. Whether the background is all genes or all time points is not stated.

## Fig. 4h: do human FAPs match mouse FAPs?

(`FAP_compare_with_Dong_etal.R`)

1. A mouse marker table with seven populations (P1.1 to P4) is read in. Gene symbols are converted to human with `nichenetr::convert_mouse_to_human_symbols`. Genes that cannot be converted are dropped.
2. Human subcutaneous FAP markers are filtered at `avg_log2FC > 0.5` and adjusted p < 0.05. Mouse markers are filtered at `> 0`.
3. A Jaccard index is computed for each pair, giving a 17 × 7 matrix drawn with `pheatmap`.

| Mouse population | Human cluster | Meaning |
|---|---|---|
| P1 (early progenitors) | sfC02 (CD55/PI16) | The starting point matches |
| P2 (committed preadipocytes) | sfC0 (APOD/CFD) | The end point matches |
| P3 (Areg, anti-adipogenic) | sfC04 | One of the route 1 transient clusters |

All matches lie on route 1. The authors' explanation for the missing route 2 matches is that the mouse studies pre-selected cells by surface markers (Lin−/SCA+) and may not have captured such cells.

This panel does two jobs. It helps name the human clusters, and it gives a second line of support for route 1.

## Fig. 4i–j: are subcutaneous and omental FAPs the same cells?

Fig. 4i is the same display as 4a for omental fat. It sets up the comparison in 4j.

Fig. 4j (`FAP_comparisons.R`, L159–170) repeats the Jaccard procedure for 17 subcutaneous against 15 omental clusters, with log2FC > 0.5 on both sides. Rows and columns are ordered by the default clustering of `pheatmap`.

| Subcutaneous | Omental | Cell type |
|---|---|---|
| sfC02 | ofC07 | CD55/PI16 progenitors |
| sfC0 | ofC02 | Committed preadipocytes |
| Route 2 clusters | ofC05, ofC09, ofC14 | Route 2-like cells exist in omental fat too |
| sfC11 | ofC03, ofC06, ofC11, ofC12 | Mesothelial-like cells |
| sfC13 | ofC10 | Both high in CD74, modest overall overlap |

**Conclusion.** The two depots hold largely the same FAP types. They differ in proportions and in single genes. Mesothelial-like cells are one cluster in subcutaneous fat and four in omental fat. DPP4 marks route 2 clusters in subcutaneous fat and mesothelial-like clusters in omental fat, so the same gene marks different cells in the two depots.

The perivascular depot was compared as well. pfC0 resembles the progenitors and pfC03 the committed preadipocytes. The other perivascular clusters have no clear counterpart.

## Fig. 4k: do the sub-clusters exist at the protein level?

Flow cytometry on stromal vascular cells from both depots, gated on CD45− CD31− CD34+ with fluorescence-minus-one controls.

| Marker | Cluster |
|---|---|
| CD55 | sfC02 / ofC07, progenitors |
| APOD | sfC0 / ofC02, committed preadipocytes |
| CD74 | sfC13 / ofC10 |
| EZR | Mesothelial-like cells |

All four are detected in both depots. The raw files are in the `FACS/FAP panel/` folder of the authors' repository.

## Where the evidence stops

- **This is not lineage tracing.** Fig. 4g uses bulk data. It shows that marker sets rise in the order "sfC02, then transient clusters, then sfC0", the order pseudotime predicts. It cannot show that single cells travel that path.
- **No sorting experiment.** sfC02 and sfC0 were not sorted and compared for adipogenic capacity.
- **The direction is assumed.** The root is picked by hand. With another root the direction changes. The experiments agree with the chosen direction, but the direction is not derived from the data.
- **Pseudotime orders cells by similarity.** It is a continuous ranking of expression profiles and not evidence of descent.
- **Route 2 "not detected in vitro" is absence of evidence.** The interpretation that culture cannot reproduce these cells is one explanation among several.
- **Fig. 4g cannot be reproduced.** The code is missing and the background is undefined.
- **Unequal thresholds in 4h.** Human markers use log2FC > 0.5 and mouse markers > 0. Mouse genes without a human symbol are dropped. Both affect the Jaccard values.
- **Fig. 4k shows one representative sample** and covers four cell types. No route 2 cluster was checked by flow cytometry.
- **Script details.** `FAP_MSL.R` has a mistyped argument (`n.cpvatponents`, an accident of replacing "om" with "pvat") and does not load monocle3 or SeuratWrappers. In `FAP_comparisons.R` the subcutaneous-against-perivascular loop sits outside the outer loop and compares only the last cluster. Fig. 4j is not affected.

[← Fig. 2–3](03-fig2-3-immune-vascular.md) · [README](../README.md) · [Next: Fig. 5 →](05-fig5-adipocytes.md)
