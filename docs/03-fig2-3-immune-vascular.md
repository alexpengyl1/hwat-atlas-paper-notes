# 03 · Fig. 2–3: immune and vascular cells

[← Fig. 1c–e](02-fig1cde-integration-benchmark.md) · [README](../README.md) · [Next: Fig. 4 →](04-fig4-fap-pseudotime.md)

**Short answer.** Both classes integrated well and were re-clustered into 11 lymphoid, 16 myeloid and 12 vascular sub-clusters. Depot differences seen in pooled cell proportions were then checked in paired bulk samples.

Scripts: `scripts/Cell_Types/Immune_Cells.R`, `scripts/Cell_Types/Vascular_Cells.R`, `scripts/Bulk_deconvolution/Deconvolution_Depots.R`.

## Contents

- [The recipe shared by Fig. 2, 3 and 4](#the-recipe-shared-by-fig-2-3-and-4)
- [Fig. 2: immune cells, panel by panel](#fig-2-immune-cells-panel-by-panel)
- [Fig. 3: vascular cells, panel by panel](#fig-3-vascular-cells-panel-by-panel)
- [How are depot differences shown?](#how-are-depot-differences-shown)
- [How were the depot differences checked in bulk data?](#how-were-the-depot-differences-checked-in-bulk-data)
- [What is concluded](#what-is-concluded)
- [Where the evidence stops](#where-the-evidence-stops)

## The recipe shared by Fig. 2, 3 and 4

Every sub-clustering in the paper follows the same three lines. Once this is clear, the cell-type scripts read quickly.

```r
obj <- RunUMAP(obj, dims = 1:10, reduction = "scvi")
obj <- FindNeighbors(obj, reduction = "scvi", dims = 1:10)
obj <- FindClusters(obj, algorithm = 1, resolution = res)   # Louvain
```

- **Everything runs on the scVI latent space.** UMAP, neighbours and clusters all use the 10 scVI dimensions.
- **Subsets are not integrated again.** When a subset is taken, for example T cells, only the UMAP, neighbours and clustering are rerun, on the same latent space.
- **Cluster numbers are sizes.** Seurat numbers clusters from largest to smallest. Names such as lyC0 to lyC10 are therefore ordered by abundance automatically.
- **Markers are written to Excel.** `FindAllMarkers(only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.15)` gives the top 20 genes per cluster for display. A separate `FindMarkers(logfc.threshold = 0.1, min.pct = 0.1)` per cluster is written to Excel, one sheet per cluster, for example `MARKER/Tcells_res065.xlsx`.
- **Colours are assigned by hand.** Four colour ramps (blue, orange, green, purple) are distributed over the clusters by biological group.

Those Excel files matter more than they look. They are the input for bulk deconvolution ([chapter 08](08-fig8-clinical.md)) and the later marker comparisons, and their cluster names become the labels of the spatial reference ([chapter 07](07-fig7-spatial.md)).

## Fig. 2: immune cells, panel by panel

The whole immune class is clustered first at resolution 0.6, which gives 19 clusters (0–18). Groups of those clusters are then taken out and re-clustered.

| Panel | What is done | Where |
|---|---|---|
| Supp. 3a–b (overview) | The 19 clusters are merged at resolution 0.2 for an overview plot. Violin genes: MS4A1, CPA3, MRC1, ITGAX, FABP4, CD3D | L20–48 |
| 2a, top: T, NK and NKT cells | Clusters 1, 2, 4 and 16 are taken out. Resolutions 0 to 1.2 are scanned and 0.65 is chosen. A comment says this value separates NK from NKT cells after validation. Result: 11 clusters, lyC0–10 | L141–170 |
| 2a, bottom: marker UMAPs | `FeaturePlot(..., cols = viridis(15), max.cutoff = 3, order = T)` | L172–174 |
| 2b: proportions by depot | Counts per cluster and depot, turned into percentages within each depot, stacked bars | L258–290 |
| 2d: myeloid cells | Clusters 0, 3, 5, 6, 7, 9, 10, 11, 13 and 14 are taken out. Resolution 0.6, with a comment that LAM, MMe and dendritic cells only separate at this value. Result: 16 clusters, myC0–15 | L293–323 |
| 2d–e: marker UMAPs | TIMP1, CALD1, FOXO1, ACACB, BCOR, CRISPLD2, MRC1, ITGAX, PPARG, CD1C, CD9, FCGR3A | L335–346 |
| 2f: myeloid proportions by depot | Same code as 2b | L351–381 |
| Supp. 3e–f | Mast cells are cluster 12 and show no clear sub-groups. B cells are cluster 15 and split in two at resolution 0.5 | L458–579 |

Composition of the immune class: monocytes, macrophages and dendritic cells 67.3%, T, NK and NKT cells 28.8%, mast cells 2.94%, B cells 0.87%, plasma cells 0.16%.

**Lymphoid sub-clusters (lyC0–10).** Th1-polarised CD4 T cells (lyC0), tissue-resident memory CD4 T cells (lyC01), early (lyC02) and late (lyC03) differentiated CD8 T cells, naive or early differentiated CD4 T cells (lyC04), NKT cells (lyC05), CD16+ NK cells (lyC06), regulatory T cells (lyC08) and CD16− NK cells (lyC09).

**Myeloid sub-clusters (myC0–15).** Several M2-like macrophage populations (myC0, 01, 04, 07, 08, 09, 11, 12, 13), of which myC08 and myC12 had not been reported in WAT before. Also lipid-associated macrophages (LAM, myC02), class 2 dendritic cells (myC03), non-classical monocytes (myC05), mixed M1/M2-like macrophages (myC06), metabolic-regulated macrophages (MMe, myC10), classical monocytes (myC14) and redox-regulatory macrophages (Mox, myC15).

## Fig. 3: vascular cells, panel by panel

| Panel | What is done |
|---|---|
| 3a: UMAP | Vascular class on the scVI latent space, resolution 0.6, 12 clusters (vC0–11). L20–42 |
| 3b: stacked violins | One representative gene per cluster: PCDH10, CD48, TYROBP, TIMP1, ACTA2, LYVE1, APOD, NCKAP5, LDLR, EFNB2, ACKR1, CA4. L45–46 |
| 3c–e: marker UMAPs | vC08 (TIMP1, ITLN1), vC09 (TYROBP, FCER1G), vC05 (CXCL14, APOD, CFD, PDGFRA). These lines are not in the script |
| 3f: proportions by depot | Same code as Fig. 2b. L61–85 |
| 3g: bulk check | Clusters 6 and 8 in `Deconvolution_Depots.R` |
| Supp. 2e: pie charts | For each cluster, the relative share of snSeq and scSeq cells. It shows that no cluster belongs to one method only. L88–111 |

**Vascular sub-clusters (vC0–11).** Blood endothelial cells: capillary (vC0), venous (vC01, vC11), arterial (vC02), mixed capillary and venous (vC03). Lymphatic endothelial cells (vC06). Vascular smooth muscle cells (vC07). Pericytes (vC04). Three unusual endothelial populations sit apart:

| Cluster | Markers | The authors' reading |
|---|---|---|
| vC08 | TIMP1, ITLN1. Shares markers with the macrophage cluster myC08 (Jaccard 0.33) | May inhibit new vessel formation |
| vC09 | TYROBP, FCER1G | Resembles blood-derived "early endothelial progenitor cells" that promote vessel formation |
| vC05 | CXCL14, APOD, CFD, and the only vascular cluster with PDGFRA | Possibly a state between endothelial cells and adipocyte precursors. "Needs further functional studies" |

## How are depot differences shown?

Panels 2b, 2f and 3f are stacked bars with one bar per depot. The code, from `Vascular_Cells.R` (the immune script does the same):

```r
cluster_prop <- vascular[[]] %>% group_by(seurat_clusters, tissue) %>% dplyr::count()
# N_sc, N_om, N_pvat = total cells of this class in each depot
cluster_prop$percent <- ifelse(cluster_prop$tissue == "sc", cluster_prop$n / N_sc * 100,
                        ifelse(cluster_prop$tissue == "om", cluster_prop$n / N_om * 100,
                                                            cluster_prop$n / N_pvat * 100))
ggplot(cluster_prop1, aes(y = percent, x = tissue, fill = seurat_clusters)) +
  geom_bar(position = "stack", stat = "identity", color = "black")
```

Each bar adds up to 100%. A segment is one sub-cluster as a share of all cells of that class in that depot.

**Example from the source data (`Figure_2b.txt`).**

| Depot | lyC0 cells | Share of T and NK cells in that depot |
|---|---|---|
| Omental | 773 | 28.5% |
| Subcutaneous | 910 | 12.3% |
| Perivascular | 139 | 52.5% |

"lyC0 is enriched in omental fat" compares shares within a depot. The absolute count is actually higher in subcutaneous fat.

The denominators pool the cells of all cohorts. No test is run per person or per cohort. That is why the bulk check below is needed.

## How were the depot differences checked in bulk data?

Two published cohorts have subcutaneous and omental bulk samples from the same people: Arner et al. with 80 pairs and Krieg et al. with 23 pairs.

1. The marker Excel files are read. The summary sheet is dropped, and per cluster only genes with `avg_log2FC > 0` and adjusted p < 0.05 are kept. They are combined into one long table with three columns: cluster, gene, log2FC (L30–70).
2. BisqueRNA's marker mode is run on each bulk expression set (L86, L235):

```r
BisqueRNA::MarkerBasedDecomposition(eset, markers, unique_markers = FALSE,
                                    weighted = TRUE, min_gene = 12, w_col = "avg_log2FC")
```

   `min_gene` is 12 for the larger cohort and 8 for the smaller one.
3. The scores of the chosen clusters are drawn as violins, subcutaneous against omental.
4. The test is a two-sided Wilcoxon signed-rank test on the pairs. It is not in the script.

What the score means, and what it cannot say, is explained in [chapter 08](08-fig8-clinical.md#what-the-score-is-and-is-not). In short: it is a relative score per sample, it can be negative, and it is not a cell proportion. The source data show this directly. For one person, the adipocyte score is −3.43 in the omental sample and +3.41 in the subcutaneous sample.

## What is concluded

| Finding | Single-cell proportions | Bulk check |
|---|---|---|
| lyC0 (Th1 CD4 T cells) is enriched in omental fat | Fig. 2b | Fig. 2c, both cohorts |
| lyC03 (late CD8 T cells) is enriched in subcutaneous fat | Fig. 2b | Fig. 2c, both cohorts |
| M2-like myC07, myC08 and myC12 are enriched in omental fat | Fig. 2f | Fig. 2g |
| LAMs (myC02) are in subcutaneous and omental fat, nearly absent in perivascular fat | Fig. 2f | — |
| Mox (myC15) are found only in subcutaneous fat | Fig. 2f | — |
| Lymphatic (vC06) and TIMP1+ (vC08) endothelial cells are enriched in omental and perivascular fat | Fig. 3f | Fig. 3g |

All lymphoid and vascular sub-clusters are present in all three depots, in different proportions. Flow cytometry confirms several lymphoid and myeloid sub-clusters at the protein level in subcutaneous fat (Supplementary Fig. 3c–d).

## Where the evidence stops

- **Mixed states could be technical.** vC05, vC08 and vC09 carry markers of preadipocytes, macrophages and blood-derived cells. The paper reads them as intermediate or mixed states. Doublets or ambient RNA could produce the same picture. The newly generated datasets were cleaned with CellBender and scDblFinder, but Methods does not say whether the public datasets were treated the same way.
- **Resolutions were picked by hand** so that known cell types separate. The script comments say so. There is no stated rule.
- **Perivascular proportions come from one dataset** (3 men).
- **The 0.33 Jaccard between vC08 and myC08** is a comparison across classes. There is no code for it.
- **A bug in the published script.** The depot code builds the totals with `aggregate(..., by = list(depot = ...))` and then reads `tissue_spot_count$tissue`, a column that does not exist under that name. Run as written, the percentages would be missing values. The source data hold normal percentages, so the version that was run must have differed from the one published.
- **Flow cytometry confirms that cell types exist.** It does not address the depot differences.

[← Fig. 1c–e](02-fig1cde-integration-benchmark.md) · [README](../README.md) · [Next: Fig. 4 →](04-fig4-fap-pseudotime.md)
