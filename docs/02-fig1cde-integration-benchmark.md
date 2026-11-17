# 02 · Fig. 1c–e: tissue composition and the choice of integration method

[← Fig. 1a–b](01-fig1ab-jaccard-network.md) · [README](../README.md) · [Next: Fig. 2–3 →](03-fig2-3-immune-vascular.md)

**Short answer.** The benchmark in Fig. 1e scores batch mixing only. scVI was chosen because it led on one of three metrics (kBET) and an outside benchmark recommended it for complex data. It did not win across the board.

## Contents

- [Fig. 1c–d: how much of the tissue is each class?](#fig-1cd-how-much-of-the-tissue-is-each-class)
- [Does Fig. 1e include all the single-cell data?](#does-fig-1e-include-all-the-single-cell-data)
- [Were the datasets annotated before integration?](#were-the-datasets-annotated-before-integration)
- [How was the benchmark run?](#how-was-the-benchmark-run)
- [The five versions of the data](#the-five-versions-of-the-data)
- [How are the three metrics computed, and from what data?](#how-are-the-three-metrics-computed-and-from-what-data)
- [What the benchmark does not measure](#what-the-benchmark-does-not-measure)
- [Result](#result)
- [Where the evidence stops](#where-the-evidence-stops)

## Fig. 1c–d: how much of the tissue is each class?

**Input.** After per-dataset clustering, every cluster was assigned by hand to one of four classes. `Source Data/Figure_1/Annotation_Massier3.txt` shows what that looks like for one dataset: cluster number, annotation, class. The percentage of each class is then computed per dataset. The result is `Figure_1c_and_1d.txt` with the columns `cohort`, `tissue`, `group`, `sum_group` (the percentage) and `method`. For example, in Massier #1 adipocytes are 33.3% and vascular cells 29.9% of all nuclei.

**Plot.** Datasets are grouped by method (snSeq, n = 8, against scSeq, n = 6) or by depot, and drawn as mean ± SD. The test is a two-sided Mann–Whitney U.

**Conclusion.** FAPs are the largest class, about 40% of cells by both methods. Adipocytes and immune cells appear in roughly equal numbers and together make up another 40%. Adipocytes can only be counted in snSeq data, because scSeq does not capture them.

**Limits.**

- One data point is one dataset, not one person, so n is 10 at most. "No difference between depots" means none was detected at that sample size.
- There is no code for these panels in the repository. Statistics and plots were probably made in GraphPad Prism, which Methods lists next to R.

The depot labels used throughout the paper and the scripts:

| Label | Depot | Location |
|---|---|---|
| `sc` | Subcutaneous | Under the skin, usually abdominal |
| `om` | Omental (a visceral depot) | The greater omentum inside the abdomen |
| `pv` or `pvat` | Perivascular | The fat layer around blood vessels |

Each cell carries one of these labels in the `tissue` column of the Seurat object. Note that `sc` in "scSeq" means single-cell, not subcutaneous.

## Does Fig. 1e include all the single-cell data?

No. It uses 15 sn/sc datasets. Two datasets were left out, and the rest were trimmed.

| Left out or trimmed | Reason |
|---|---|
| The spatial dataset (Bäckdahl) | A spot is not a single cell and needs deconvolution. It is analysed separately |
| Hildreth | It failed to mix with the other datasets in a trial integration with scANVI (Supplementary Fig. 1c). It is excluded from everything after that |
| Cell numbers per dataset | Each remaining dataset was capped to a similar number of objects, so that large datasets would not dominate (Supplementary Fig. 2c–d) |
| Adipocytes from scSeq | They do not exist. The adipocyte class comes only from the seven snSeq datasets: 4 subcutaneous, 2 omental, 1 perivascular |

## Were the datasets annotated before integration?

Yes, but only to the four classes.

- Each dataset was clustered on its own, and every cluster was assigned to adipocyte, FAP, vascular or immune. The network in Fig. 1b was used to check that assignment.
- Fine sub-clusters, such as the individual macrophage types, did not exist yet. They come from re-clustering after integration.

The four-class labels were used twice:

1. As input labels for the scANVI trial. scANVI is a semi-supervised version of scVI that uses known labels to guide integration.
2. To split the cells, so that each class could be integrated separately.

The scANVI trial (`scripts/Integration/scANVI.Rmd`) first trains scVI with `n_layers = 2`, `n_latent = 30` and a negative binomial likelihood, then builds scANVI on top of it with `labels_key = "cluster_high"` and trains for 20 epochs with 100 samples per label. With Hildreth removed, the scANVI embedding reproduced the four classes of the network (Supplementary Fig. 1d).

## How was the benchmark run?

1. **Merge** the 15 datasets into one Seurat object.
2. **Run everything twice:** on all cells together, and on each of the four classes separately.
3. **Define the batch.** The batch variable is `orig.ident`. Different depots or methods within one study count as separate batches, so the three Emont datasets are three batches.
4. **Produce five versions** of the data: uncorrected, rPCA, Harmony, BBKNN and scVI.
5. **Score** each version against three confounders (method, depot, cohort) with three metrics.

Fig. 1e shows those scores.

## The five versions of the data

| Method | Script | Key settings |
|---|---|---|
| Uncorrected | — | The merged data without correction |
| Harmony | `Harmony.R` | `RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT", reduction = "ica")`. UMAP on 30 dimensions with `n.neighbors = 50`. Clustering at resolution 1 |
| rPCA | `rPCA.R` | Split by `orig.ident`, `SCTransform(method = "glmGamPoi")` on each part, 3,000 integration features, `FindIntegrationAnchors(reduction = "rpca", k.anchor = 20, dims = 1:30)`, then `IntegrateData` |
| BBKNN | `BBKNN_integration.ipynb` | Convert to h5ad, `bbknn.bbknn(batch_key = 'orig.ident')`, Leiden clustering at 0.3, `ridge_regression(confounder_key = 'leiden')`, new PCA, BBKNN again |
| scVI | `scVI.R` | 2,000 variable genes, counts as input, `batch_key = 'orig.ident'`, default model with 10 latent dimensions |

The scVI script in full, lightly condensed:

```r
seurat   <- readRDS("Seurat_Merged.rds")
top2000  <- head(VariableFeatures(seurat), 2000)
seurat_top <- seurat[top2000]
adata <- convertFormat(seurat_top, from = "seurat", to = "anndata",
                       main_layer = "counts", drop_single_values = FALSE)

scvi$model$SCVI$setup_anndata(adata, batch_key = 'orig.ident')
# script comment: further covariate keys are possible, did not improve integration for this project
model  <- scvi$model$SCVI(adata)     # all defaults
model$train()
latent <- as.matrix(model$get_latent_representation())
rownames(latent) <- colnames(seurat_top)

seurat[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_",
                                         assay = DefaultAssay(seurat_top))
seurat <- RunUMAP(seurat, dims = 1:10, reduction = "scvi", n.components = 2)
```

What scVI gives back is a table with one row per cell and 10 columns. That table is stored in the Seurat object as a reduction called `scvi`, and every later UMAP and clustering runs on it.

## How are the three metrics computed, and from what data?

None of them uses the expression matrix or cell-type labels. They use the **output of integration** plus the **batch label** of each cell. The batch label is one of three things, each scored separately: sequencing method, depot or cohort.

| Metric | R package | Input | Question it asks | Better when |
|---|---|---|---|---|
| ARI | `pdfCluster` | Cluster labels after integration, batch labels | Do the clusters follow the batches? | Lower |
| LISI | `lisi` | UMAP coordinates, batch labels | How many batches surround each cell? | Higher |
| kBET | `kBET` | The integrated low-dimensional embedding (for example the scVI latent space), batch labels | Does the local batch mix match the overall mix? | Higher acceptance rate |

The code that computes them is not in the repository. The descriptions below follow Methods and the algorithms of the three packages.

### ARI (adjusted Rand index)

The same cells carry two groupings: cluster and batch. ARI measures how well the two groupings agree. It is 1 when they are identical and about 0 when they are unrelated. Here a value near 0 is wanted.

The calculation counts pairs of cells. With $n_{ij}$ cells in cluster $i$ and batch $j$, row totals $a_i$, column totals $b_j$ and $n$ cells in all:

$$\mathrm{ARI}=\frac{\sum_{ij}\binom{n_{ij}}{2}-\dfrac{\sum_i\binom{a_i}{2}\sum_j\binom{b_j}{2}}{\binom{n}{2}}}{\dfrac{1}{2}\left[\sum_i\binom{a_i}{2}+\sum_j\binom{b_j}{2}\right]-\dfrac{\sum_i\binom{a_i}{2}\sum_j\binom{b_j}{2}}{\binom{n}{2}}}$$

**Illustrative example.** Six cells come from batches A A A B B B.

| Case | Cluster labels | ARI | Meaning |
|---|---|---|---|
| Failed integration | 1 1 1 2 2 2 | 1 | Clusters copy the batches |
| Good mixing | 1 2 1 2 1 2 | −0.11 | Clusters are unrelated to batch |

For the second case: cluster 1 holds 2 A and 1 B, cluster 2 holds 1 A and 2 B. Pairs that share both cluster and batch: 1 + 0 + 0 + 1 = 2. Each cluster and each batch has 3 cells, so both sums of pairs are 3 + 3 = 6. There are 15 pairs in total. The expected value is 6 × 6 / 15 = 2.4 and the maximum is 6. ARI = (2 − 2.4) / (6 − 2.4) = −0.11.

### LISI (local inverse Simpson index)

For each cell, look at its neighbours on the UMAP and take the proportion $p_b$ of each batch $b$ among them.

$$\mathrm{LISI}=\frac{1}{\sum_b p_b^2}$$

The values of all cells are averaged. The number can be read as "how many batches are effectively present around a cell".

| Batch mix among neighbours | Calculation | LISI |
|---|---|---|
| Half A, half B | 1 / (0.25 + 0.25) | 2 |
| 90% A, 10% B | 1 / (0.81 + 0.01) | 1.22 |
| All A | 1 / 1 | 1 |

The maximum equals the number of batches. It is 2 for method, 3 for depot and much higher for cohort, so LISI values for different confounders cannot be compared with each other.

### kBET (k-nearest-neighbour batch effect test)

A sample of cells is tested. For each one, take its $k$ nearest neighbours in the integrated embedding and compare their batch composition with the composition of the whole dataset, using a chi-square test. A significant deviation (p < 0.05) counts as a rejection. The figure shows the acceptance rate, which is 1 minus the share of rejected cells.

**Illustrative example.** The whole dataset is 50% A and 50% B. With 20 neighbours, 10 of each are expected.

| Neighbours observed | Chi-square | p | Outcome |
|---|---|---|---|
| 11 A, 9 B | (1² / 10) + (1² / 10) = 0.2 | 0.65 | Accepted |
| 18 A, 2 B | (8² / 10) + (8² / 10) = 12.8 | 0.0004 | Rejected |

### How the three differ

- ARI works at the level of clusters and depends on the clustering resolution.
- LISI is computed on the UMAP, which distorts distances.
- kBET is computed on the integrated embedding and is the strictest of the three. That is why its absolute values are low for every method.

## What the benchmark does not measure

It measures whether cells from **different sources are mixed**. It does not measure whether cells of the **same type stay together**.

Suppose the immune class holds macrophages and T cells from two datasets, A and B.

- **A good result:** macrophages of A and B mix, T cells of A and B mix, and macrophages stay apart from T cells.
- **An over-corrected result:** all four groups are blended into one cloud.

Both score well on all three metrics, because in both the batches are evenly mixed.

The paper covers the missing half indirectly:

1. Cells are split into four classes before integration. Within one class the cells are assumed to be similar enough that over-mixing does less harm.
2. The marker genes of the resulting clusters are inspected. Methods mentions this in one sentence and reports no number.

Integration benchmarks usually report both kinds of score: batch mixing, and conservation of biology (for example LISI or ARI computed on cell-type labels). This paper reports the first kind only. So "scVI is best" should be read as "scVI removes batch structure most thoroughly by kBET", not as "scVI preserves cell-type differences best".

## Result

Values from the source data file (`Source Data/Figure_1/Figure_1e.xlsx`). Integration is scored by cohort unless stated.

| Metric | Uncorrected | Harmony | BBKNN | rPCA | scVI |
|---|---|---|---|---|---|
| ARI, all cells (lower is better) | 0.387 | 0.151 | 0.500 | **0.078** | 0.147 |
| LISI, all cells (higher is better) | 1.53 | 3.45 | 2.53 | **4.42** | 3.14 |
| kBET acceptance, immune cells | 0.0005 | **0.044** | 0 | 0.012 | 0.0045 |
| kBET acceptance by depot, all cells | 0.008 | 0.035 | 0.065 | 0.002 | **0.304** |

- rPCA leads on ARI and LISI in most columns.
- scVI leads clearly on kBET by method and by depot.
- Splitting by class helps. LISI by cohort for scVI rises from 3.14 on all cells to 4.30 on vascular cells.
- Adipocytes integrate badly whatever the method. Their best kBET acceptance by cohort is 0.023 (scVI), and it is 0 for the other methods. This is the starting point of [chapter 05](05-fig5-adipocytes.md).

The paper's wording is that Harmony did well for FAPs, that scVI and rPCA gave the best results across scores, and that kBET rejection rates were lowest for scVI. Together with a published benchmark that favours scVI for complex data, this led the authors to apply scVI to each class separately.

## Where the evidence stops

- **kBET acceptance is low everywhere,** 0.475 at best. "Integrated well" is relative to the uncorrected data.
- **Duplicate columns in the source data.** In the three ARI sheets, the all-cells and FAP columns are identical to many decimals, and so are the vascular and immune columns. Whether that is a coincidence or a copying error cannot be told from the file. The numbers should be checked with the authors before being reused.
- **The metric code is missing.** The clustering resolution behind ARI, and the neighbourhood size and sampling fraction for kBET, are unknown. All of them change the values, so Fig. 1e cannot be reproduced exactly.
- **Merging, splitting and capping are missing too.** `scVI.R` starts from an already merged `Seurat_Merged.rds`. How many cells each dataset was capped to can only be read from the supplementary figure. A file name in another script (`fibroblasts_8000.rds`) hints that the cap for FAPs may have been 8,000. That is a guess.
- **Which counts went into scVI?** `convertFormat(..., main_layer = "counts")` is called without naming an assay. If the default assay of the object was SCT, scVI received SCT-corrected counts and not raw UMI counts. The script does not show which it was.
- **No robustness check.** The paper does not show whether the sub-clusters would look the same under rPCA, the method that scored better on two of the three metrics.

[← Fig. 1a–b](01-fig1ab-jaccard-network.md) · [README](../README.md) · [Next: Fig. 2–3 →](03-fig2-3-immune-vascular.md)
