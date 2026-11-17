# 07 · Fig. 7: spatial deconvolution and co-localisation

[← Fig. 6](06-fig6-cell-communication.md) · [README](../README.md) · [Next: Fig. 8 →](08-fig8-clinical.md)

**Short answer.** Mapping the sub-clusters onto 10 tissue sections shows that myeloid cells gather in patches and that two FAP sub-clusters keep different company: sfC12 sits with endothelial cells and sfC08 with lipid-associated macrophages.

Scripts: everything under `scripts/Spatial_Deconvolution/`.

## Contents

- [What goes in?](#what-goes-in)
- [Is the reference the integrated single-cell data?](#is-the-reference-the-integrated-single-cell-data)
- [How is the deconvolution run?](#how-is-the-deconvolution-run)
- [What does "run at two levels of detail" mean?](#what-does-run-at-two-levels-of-detail-mean)
- [Fig. 7a: myeloid patches](#fig-7a-myeloid-patches)
- [Fig. 7c: how is co-localisation computed?](#fig-7c-how-is-co-localisation-computed)
- [Fig. 7d: two cell types on one section](#fig-7d-two-cell-types-on-one-section)
- [What is concluded](#what-is-concluded)
- [Where the evidence stops](#where-the-evidence-stops)

## What goes in?

| Input | Content |
|---|---|
| Spatial data | 10 Visium sections of subcutaneous fat from 10 people, from an earlier study (Bäckdahl et al. 2021). Section numbers 42, 44, 46, 48–52, 54, 55 |
| Reference | `WAT_scRNA_MTX.h5ad`, with cell-type labels in the column `celltype_new` |

A Visium spot is 55 µm across and captures the RNA of whatever cells lie on it. Deconvolution estimates how much of each cell type a spot contains, given a reference expression profile per cell type.

## Is the reference the integrated single-cell data?

Its labels are. Its values are not.

| Aspect | What the code does |
|---|---|
| Labels | The final sub-cluster names from Fig. 2–4, such as `myC02_LAM`, `sfC08_FAP` and `vC01_Blood_EC_vein`, plus one pooled `adipose` class. These names exist only after integration, re-clustering and annotation |
| Cells | Not all of them. At most 1,500 randomly chosen cells per sub-cluster. Smaller sub-clusters are used whole |
| Values | The expression matrix stored in the h5ad file. The 10-dimensional scVI latent space is not used |
| Batch | `Cell2location_pipeline.py` passes only `labels_key` to the model (line 48). No `batch_key` is given |

The last row matters. Integration is used only to **name** the cells. When the reference profile of each sub-cluster is estimated, cells from different studies and from snSeq and scSeq are pooled without batch correction. cell2location accepts a batch variable for exactly this purpose, and it is not used here.

**Illustrative example.** Suppose `myC02_LAM` has 5,000 cells from 8 datasets in the integrated data.

1. 1,500 of them are drawn at random.
2. Their expression is used to estimate one reference profile for LAMs, one value per gene.
3. That profile, together with the profiles of all other sub-clusters, is used to explain the expression of each spot. The output is the abundance of LAMs in every spot.

**What cannot be confirmed from the repository.**

- The reference file and the code that built it, the sub-sampling included, are missing.
- Whether the reference holds only subcutaneous cells. The sections are subcutaneous, the co-localisation script uses subcutaneous FAP names (`sfC`), and the file name resembles the subcutaneous object used for CellChat. So it probably does, but that is an inference.
- Whether the matrix holds raw counts. cell2location expects raw counts, and the script reads the matrix without converting it.

## How is the deconvolution run?

`Spatial_Deconvolution.py` takes a section number and calls `Benchmarking/DeconvolutionSpot.py`, which starts one script per tool. The structure and comments of this wrapper suggest it was adapted from a published benchmarking framework for spatial deconvolution.

```python
RNA_h5ad     = '~/Spatial_Deconvolution/WAT_scRNA_MTX.h5ad'
Spatial_h5ad = '~/Spatial_Deconvolution/st_s' + sys.argv[1] + '_MTX.h5ad'
celltype_key = 'celltype_new'

test = DeconvolutionSpot.Deconvolutions(RNA_h5ad=RNA_h5ad, Spatial_h5ad=Spatial_h5ad,
                                        celltype_key=celltype_key, output_path=output_path)
Result = test.Dencon(['Cell2location', 'Stereoscope', 'Tangram', 'DestVI'])
# RCTD and SPOTlight run from h5Seurat files in a second call
```

**cell2location**, the tool used for the main figure (`Codes/Deconvolution/Cell2location_pipeline.py`, condensed):

```python
# 1. reference: one expression signature per sub-cluster
selected = filter_genes(adata_ref, cell_count_cutoff=5,
                        cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
RegressionModel.setup_anndata(adata=adata_ref, labels_key=celltype_key)
mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# 2. spatial: explain every spot as a mixture of those signatures
Cell2location.setup_anndata(adata=adata_vis)
mod = Cell2location(adata_vis, cell_state_df=inf_aver,
                    N_cells_per_location=30, detection_alpha=200)
mod.train(max_epochs=30000, ...)

# 3. output
adata_vis.obsm['q05_cell_abundance_w_sf'].to_csv(output_file_path + '/cell2location_result.csv')
```

| Setting | Meaning |
|---|---|
| `N_cells_per_location = 30` | A prior for how many cells a spot holds. It depends on the tissue |
| `detection_alpha = 200` | How much the RNA detection efficiency is allowed to vary between spots |
| `q05_cell_abundance_w_sf` | The 5% quantile of the estimated abundance of each cell type in each spot. A deliberately cautious value: the model is 95% confident that the true abundance is at least this high |

**The other five tools** are used only to check agreement (Supplementary Fig. 7a):

| Tool | Key settings |
|---|---|
| RCTD | `doublet_mode = 'full'`. Reference types with fewer than 25 cells are dropped |
| Tangram | `mode = 'clusters'`. Training genes are the union of the top 200 genes per type |
| DestVI | `CondSCVI`, then `DestVI.from_rna_model`, 2,500 epochs on the spatial data |
| stereoscope | 7,000 variable genes. 100 epochs on the reference, 10,000 on the spatial data |
| SPOTlight | 100 cells per type, 3,000 variable genes, non-negative matrix factorisation |

The comparison (`Compare_Spatial_Deconvolution_Tools.R`): each tool's result is rescaled to proportions per spot, values below 0.03 are set to 0, and within each section all tool-by-cell-type columns are correlated with each other. The 10 correlation matrices are averaged. If the same cell type correlates highly across tools, the tools agree.

## What does "run at two levels of detail" mean?

Deconvolution was run twice on the same sections. The two runs differ only in how finely the reference cells are labelled.

| | Coarse run | Fine run |
|---|---|---|
| Labels | 5 classes: adipose, FAP, lymphocyte, myeloid, vascular | Several dozen sub-clusters, such as `adipose`, `lyC00_CD4_Th1`, `myC02_LAM`, `sfC08_FAP` |
| Each spot gets | 5 numbers | Several dozen numbers |
| Tools | All six | cell2location only |
| Source data file | `Figure_7a_and_ED_Figure_7a.xlsx`, one sheet per section and tool | `Figure_7cd_and_ED_Figure_7bc.xlsx`, one sheet per section |
| Used for | The overview in Fig. 7a and the six-tool comparison | Fig. 7c, 7d and Supplementary Fig. 7b–c |
| Answers | Are there many myeloid cells in this spot? | Are they LAMs or M2-like? Which FAP sits next to them? |

**Why two runs.** Classes differ strongly in expression, so they are easy to tell apart and suited to comparing tools. Sub-clusters within a class are similar to each other, which makes them harder to separate, but the questions of Fig. 7 can only be answered at that level.

**The catch.** The agreement between six tools was shown for the five classes only. The sub-cluster results come from cell2location alone. In addition, the repository script names a single label column, `celltype_new`, which corresponds to the fine run. Which column or file was used for the coarse run cannot be seen in the code.

## Fig. 7a: myeloid patches

(`Compare_Spatial_Deconvolution_Tools.R`, second half, and the "Zoom in" block of `Co_localization.R`)

1. **Whole section.** `SpatialFeaturePlot(st_s52, features = "cell2location_Myeloid", image.alpha = 0, min.cutoff = "q1", max.cutoff = "q99")` for two people, sections 52 and 55.
2. **Zoomed region.** The section is cropped by pixel coordinates, for example rows 4200–4800 and columns 7000–7600 for section 55.
3. **Box plot.** For all spots inside the crop, the sub-cluster abundances of the myeloid types (`cell2location_myC*`) are drawn with `ggpubr::ggboxplot`, sorted by column total.

The patches are rich in classical monocytes (myC14), LAMs, MMe and Mox, and poor in M2-like macrophages. The two regions are boxes picked by eye. Fig. 7b, immunofluorescence for the LAM markers TREM2 and CD9, is an independent check that LAMs form clusters.

## Fig. 7c: how is co-localisation computed?

(`Co_localization.R`)

1. In each section, the sub-cluster abundances of every spot are rescaled to proportions that sum to 1. Proportions below 0.03 are set to 0. Cell types that are 0 everywhere are removed.
2. **Within each section**, the Pearson correlation between every pair of cell types is computed across spots.
3. The 10 correlation matrices are averaged. The script starts from a matrix of zeros, adds each section's matrix and divides by 10.
4. `pheatmap(clustering_method = "ward.D2", cutree_rows = 4, cutree_cols = 4)`. Fig. 7c is a zoom on the first 16 types of the tree, and the full map is Supplementary Fig. 7b.

$$r_{AB}=\frac{\sum_s (a_s-\bar a)(b_s-\bar b)}{\sqrt{\sum_s (a_s-\bar a)^2}\,\sqrt{\sum_s (b_s-\bar b)^2}}$$

$a_s$ and $b_s$ are the proportions of cell types A and B in spot $s$.

**Illustrative example with four spots.**

| Spot | A (raw) | A (after threshold) | B |
|---|---|---|---|
| 1 | 0.40 | 0.40 | 0.30 |
| 2 | 0.30 | 0.30 | 0.20 |
| 3 | 0.10 | 0.10 | 0.10 |
| 4 | 0.02 | 0 | 0 |

The means are 0.20 for A and 0.15 for B. The products of the deviations add up to 0.03 + 0.005 + 0.005 + 0.03 = 0.07. The denominator is √(0.10 × 0.05) ≈ 0.0707. So r ≈ 0.99. The two cell types rise and fall in the same spots, which is read as co-localisation.

**What "co-localised" means here.** Present in the same 55 µm spot. It does not mean present in neighbouring spots.

Supplementary Fig. 7c draws, for cell-type pairs with |r| > 0.1, the correlation of each person separately. It shows that the associations do not come from one or two sections.

## Fig. 7d: two cell types on one section

The last lines of `Co_localization.R` use a trick:

```r
st_s49@reductions$umap@cell.embeddings[, 1] <-  st_s49$col
st_s49@reductions$umap@cell.embeddings[, 2] <- -st_s49$row
FeaturePlot(st_s49, features = c("cell2location_sfC08_FAP", "cell2location_myC02_LAM"),
            blend = TRUE, cols = c("darkblue", "green", "magenta"),
            min.cutoff = "q1", max.cutoff = "q99", blend.threshold = 0) & DarkTheme()
```

The UMAP coordinates are overwritten with the row and column positions of the Visium array. `FeaturePlot(blend = TRUE)` can then overlay two abundances in space, one in green and one in magenta, with mixed colour where both are present. Three pairs are drawn: sfC08 with sfC12, sfC08 with LAMs, and sfC12 with venous endothelial cells (vC01). Only one section (49) is shown.

## What is concluded

- **Uneven distribution (7a–b).** FAPs, lymphocytes and adipocytes are spread through the tissue. Vascular and myeloid cells are concentrated in specific areas. The myeloid patches hold classical monocytes and M1-like macrophages (LAM, MMe, Mox) and lack M2-like macrophages.
- **FAP neighbours (7c–e).** The correlation map has groups made of endothelial cells and of two sets of FAPs. sfC12 is found close to endothelial cells, confirmed by staining for its marker SLIT2 together with the endothelial protein CD31. The other FAPs relate negatively to vascular cells. sfC08 is found close to LAMs.
- **Interpretation.** The authors read this as FAPs having specific positions in the tissue, "possibly to form different types of adipogenic niches".

## Where the evidence stops

- **One tool at sub-cluster level.** The sub-cluster results rest on cell2location alone, plus immunofluorescence for two findings: LAM clustering, and sfC12 next to endothelial cells.
- **Proportions force negative correlations.** The proportions in a spot sum to 1, so a high share of one type pushes the others down. Part of the negative correlation between most FAPs and vascular cells may come from that arithmetic and not from biology. This caution is mine, not the paper's. The positive correlations are less affected.
- **The prior may not fit fat.** A spot is 55 µm across, and mature adipocytes are often around 100 µm. Whether 30 cells per spot is a suitable prior for this tissue is not discussed. Methods says default settings were used.
- **No batch term in the reference.** See above.
- **Hand-picked regions.** The patches in 7a are two boxes chosen by eye, and 7d shows a single section.
- **Three or four clusters.** The text reports three clusters in the correlation map. The script cuts the tree into four.
- **A slip in the job script.** `work.sh` lists section 51 twice and omits section 48, yet the downstream scripts read results for section 48.
- **Scope.** Only subcutaneous fat from 10 people was analysed.

[← Fig. 6](06-fig6-cell-communication.md) · [README](../README.md) · [Next: Fig. 8 →](08-fig8-clinical.md)
