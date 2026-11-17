# 01 · Fig. 1a–b: can the datasets be combined?

[← README](../README.md) · [Next: Fig. 1c–e →](02-fig1cde-integration-benchmark.md)

**Short answer.** Yes, for 16 of the 17 datasets. Before any batch correction, clusters from different studies share enough marker genes to link up into four classes. Hildreth is the exception. This is the bet the whole paper rests on: if independently defined clusters did not match, "cluster first, then integrate by class" would have no footing.

Both panels come from one script: `scripts/Network/Jaccard_and_Network_Fig1.R`.

## Contents

- [What goes in](#what-goes-in)
- [Is every cluster compared with every cluster of the other datasets?](#is-every-cluster-compared-with-every-cluster-of-the-other-datasets)
- [How is the Jaccard index computed?](#how-is-the-jaccard-index-computed)
- [How does one box per dataset come out?](#how-does-one-box-per-dataset-come-out)
- [How is the network built?](#how-is-the-network-built)
- [What do node size, colour and edges mean?](#what-do-node-size-colour-and-edges-mean)
- [What is a `*_Counts.txt` file?](#what-is-a-_countstxt-file)
- [Are edges thresholded? Is five the maximum?](#are-edges-thresholded-is-five-the-maximum)
- [What are the two sides of an overlap?](#what-are-the-two-sides-of-an-overlap)
- [What the network shows](#what-the-network-shows)
- [Where the evidence stops](#where-the-evidence-stops)

## What goes in

Each dataset was first processed on its own with Seurat. That gives, per dataset, a set of clusters and one marker gene table per cluster, saved as one Excel sheet per cluster. The script reads these Excel files for 17 datasets, lettered A to Q (the full list is in the [README](../README.md#the-17-datasets)).

A cluster is then named by letter and number, for example `A_3` is cluster 3 of dataset A.

For each cluster the script keeps only significantly up-regulated genes (L47–89):

```r
FC <- 0   # comment in the script: different FC cutoffs were tested
markers <- markers[markers$avg_log2FC > FC & markers$p_val_adj < 0.05, ]
```

The left part of Fig. 1a (subjects, sex, age and BMI ranges, object counts) is not computed. It is a table copied from the cohort descriptions.

## Is every cluster compared with every cluster of the other datasets?

Yes.

- 17 datasets give 17 × 16 / 2 = **136 dataset pairs**.
- The 397 clusters give about **73,700 cross-dataset cluster pairs**.
- Clusters from the same dataset are never compared with each other.

The comparison is written out by hand as roughly 1,700 lines of nested loops (L100–1793), one block per dataset pair. For each cluster pair it records three numbers:

| Column | Meaning |
|---|---|
| `Number` | How many marker genes the two clusters share |
| `Percent_1` | `Number` as a percentage of the first cluster's marker list |
| `Percent_2` | `Number` as a percentage of the second cluster's marker list |

## How is the Jaccard index computed?

$$J(A,B)=\frac{|A\cap B|}{|A|+|B|-|A\cap B|}$$

$A$ and $B$ are the marker gene sets of two clusters. The index is 0 when they share nothing and 1 when they are identical.

**Worked example.** Cluster `A_3` has 200 markers, cluster `B_2` has 150, and they share 60.

$$J = \frac{60}{200 + 150 - 60} = \frac{60}{290} \approx 0.21$$

The script does not keep the marker totals $|A|$ and $|B|$. It recovers them from the stored percentages (L1835–1850):

```r
output3$total   <- output3$Number / output3$Percent_1 * 100   # = |A|
output3$total_2 <- output3$Number / output3$Percent_2 * 100   # = |B|
jaccard2$jaccard <- jaccard2$Number /
                    (jaccard2$total + jaccard2$total_2 - jaccard2$Number)
```

In the example, 60 / 30% × 100 = 200 and 60 / 40% × 100 = 150, which gives back the two totals. When two clusters share no gene, the division is 0 / 0. To cover those cases, the script then sets the total in every row of a cluster to the largest value seen for that cluster (`max(temp$total, na.rm = T)`), before the line that computes the index.

## How does one box per dataset come out?

Every cluster has a Jaccard value against every cluster of the other 16 datasets, which is hundreds of values. The script reduces them in two steps (L1853–1863):

```r
j_1 <- jaccard2[, c(1, 8)]; colnames(j_1) <- c("cluster", "jaccard")   # pair seen from cluster 1
j_2 <- jaccard2[, c(2, 8)]; colnames(j_2) <- c("cluster", "jaccard")   # pair seen from cluster 2
jac_plot <- rbind(j_1, j_2)
jac_plot <- jac_plot %>% group_by(cluster) %>% slice_max(jaccard, n = 5)
jac_plot$cohort <- substr(jac_plot$cluster, 1, 1)                      # the dataset letter
ggplot(data = jac_plot, aes(x = cohort, y = jaccard)) +
  geom_boxplot(outlier.shape = NA) + theme_classic() + coord_flip()
```

1. Each cluster keeps its **five highest** Jaccard values, whichever datasets they come from.
2. The values are grouped by the dataset the cluster belongs to and drawn as a box.

A dataset with 20 clusters therefore contributes 20 × 5 = 100 values. The source data file `Figure_1a_Jaccard.txt` has 1,985 rows, which is 397 clusters × 5.

**What the box answers.** How similar is the best counterpart that this dataset's clusters find elsewhere? A high box means the clusters of this dataset have close look-alikes in other datasets.

**Result.** Dataset medians run from 0.40 (Emont, subcutaneous snSeq) down to 0.16 (the spatial dataset). Hildreth sits at 0.04. That gap is the first piece of evidence for dropping Hildreth. The per-dataset medians are in the [README table](../README.md#the-17-datasets).

A detail in the plotting code: the `geom_jitter` line that would overlay the single points is not attached to the plot, because the `+` before it is missing.

## How is the network built?

The second half of the script (L1866–1927) turns the same overlap table into a network.

```r
output3$avgPERCENT <- (output3$Percent_1 + output3$Percent_2) / 2

# remove Hildreth dataset as previous mapping showed very low overlap
output3 <- output3 %>% filter(!grepl("G_\\d", output3$Cohort_1))
output3 <- output3 %>% filter(!grepl("G_\\d", output3$Cohort_2))
nodes2  <- unique(c(output3$Cohort_1, output3$Cohort_2))

output_PERCENT <- output3[NULL, ]
for (i in 1:length(nodes2)) {
  # every row in which this node appears, on either side
  output_temp <- output3[which(output3$Cohort_1 == nodes2[[i]] | output3$Cohort_2 == nodes2[[i]]), ]
  output_temp <- output_temp %>% top_n(5, avgPERCENT)       # its five best partners
  output_PERCENT <- rbind(output_PERCENT, output_temp)
}
output_PERCENT <- output_PERCENT[which(output_PERCENT$Percent_1 > 15 | output_PERCENT$Percent_2 > 15), ]
output_PERCENT <- output_PERCENT[which(output_PERCENT$Percent_1 > 5 & output_PERCENT$Percent_2 > 5), ]
output_PERCENT <- output_PERCENT[!duplicated(output_PERCENT), ]

# edges: weight = Number (shared gene count)
routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
cytoscape <- createNetworkFromIgraph(routes_igraph,
                                     title = "TOP5 PERCENT NO HILDRETH FILTERED",
                                     collection = "Meta New")
# final visual adjustments and network clustering using cytoscape
```

R builds the graph and pushes it to Cytoscape. Layout, colouring and grouping are then done by hand in Cytoscape, as the last comment says. The finished session is `Source Data/Figure_1/Figure_1b.cys`.

Comments in the script record what was tried and rejected:

- `### NUMBER not good selection`: ranking edges by the raw shared-gene count was abandoned.
- "If no cutoff, network is overloaded with very weak connections": different percentage cutoffs were tested before settling on 15% and 5%.

## What do node size, colour and edges mean?

A `.cys` file is a zip archive. Unpacking it gives the node table, the edge table and the visual style, so the mappings can be read directly.

| Element | Based on | Detail |
|---|---|---|
| Node | One cluster of one dataset | Hildreth and the spatial dataset are left out, which leaves 15 datasets |
| Node size | `cells_percent`: the cluster's share of cells within its own dataset | Continuous mapping: 0% → size 100, 40% and above → size 500 |
| Node fill colour | `color`: a hex code assigned per cluster in advance | Passed straight through. Hue marks the class, shade separates clusters of one dataset |
| Node border | `annot_group`: the class | Discrete mapping |
| Edge | A marker overlap that passed the selection rules | `weight` is the shared gene count. In the style used by the main view, edge width is not mapped to anything |
| Layout | Cytoscape | Arranged by hand |

The colour ramps, one per class:

| Class | Light end | Dark end |
|---|---|---|
| Immune (blue) | `#C6DBEF` | `#075A84` |
| Vascular (green) | `#A6DBBB` | `#359566` |
| FAP (orange) | `#F3E55C` | `#E8602D` |
| Adipocyte (purple) | `#ECD9F1` | `#967BCE` |

So a large node is a cluster that makes up a big part of its dataset, and nodes of the same hue were judged by the original annotation to be the same class.

## What is a `*_Counts.txt` file?

Node attributes are read from one such file per dataset (16 files, none for the spatial dataset). The files are not in the repository. Their content was recovered from the node table inside the Cytoscape session.

Each row is one cluster. For Massier #3 the reconstructed file starts and ends like this:

| seurat_clusters | Freq | percent | cluster_anno | cluster_group | cluster_color | cohort |
|---|---|---|---|---|---|---|
| 0 | 2307 | 27.23 | Endothelial | Vascular/Endo | `#A6DBBB` | Lipidlab_1S |
| … | … | … | … | … | … | … |
| 14 | 52 | 0.61 | Bcell | Immune | `#075A84` | Lipidlab_1S |

| Column | Meaning | How it is obtained |
|---|---|---|
| `seurat_clusters` | Cluster number | Per-dataset clustering |
| `Freq` | Number of cells in the cluster | A count table of cluster labels. For Massier #3 it sums to 8,471, the dataset's total |
| `percent` | `Freq` as a percentage of the dataset | 2307 / 8471 = 27.23% |
| `cluster_anno` | The original annotation | Assigned by hand from marker genes |
| `cluster_group` | One of the four classes | Assigned by hand |
| `cluster_color` | Hex colour | Taken from the class ramp above |
| `cohort` | Dataset name | — |

The matching hand annotation for this dataset is in the repository as `Source Data/Figure_1/Annotation_Massier3.txt` (cluster number, annotation, class).

## Are edges thresholded? Is five the maximum?

The script applies three rules, in this order:

1. **Top five.** Each node keeps its five best partners, ranked by `avgPERCENT`, the average of the two overlap percentages.
2. **One side above 15%.** At least one of `Percent_1` and `Percent_2` must exceed 15%.
3. **Both sides above 5%.** Both percentages must exceed 5%.

Duplicates are then removed, because the same pair can be picked from both ends.

Five is **not** a maximum per node, and it is not a minimum either.

**Illustrative example.**

- `A_3` ranks all its partners and keeps the best five. All five pass rules 2 and 3, so `A_3` has five edges of its own choosing.
- `N_6` is a cluster of another dataset. It is not among `A_3`'s five best partners, but `A_3` is among `N_6`'s five best. That adds a sixth edge to `A_3`. A node that many others pick ends up with many edges.
- `F_7` comes from a small dataset with short marker lists. Two of its five best partners have 12% and 3% on the two sides. Rule 2 removes them, so `F_7` ends up with three edges.

In the saved network, the number of edges per node runs from 5 to 32, with a median of 6.

## What are the two sides of an overlap?

The shared gene count is one number. Each cluster has its own marker total, so the count turns into two percentages.

| | Calculation | Result |
|---|---|---|
| Side of `A_3` | 60 / 200 | 30% |
| Side of `B_2` | 60 / 150 | 40% |
| Average (`avgPERCENT`) | (30 + 40) / 2 | 35% |

Read it as: 30% of `A_3`'s markers are also markers of `B_2`, and 40% of `B_2`'s markers are also markers of `A_3`. The pair passes rule 2 (one side above 15%) and rule 3 (both above 5%).

**Why both sides matter.** Marker lists differ a lot in length. A small cluster with 20 markers and a large one with 1,000 markers may share 20 genes.

| | Calculation | Result |
|---|---|---|
| Side of the small cluster | 20 / 20 | 100% |
| Side of the large cluster | 20 / 1,000 | 2% |

Seen from the small cluster, the match looks perfect. Seen from the large one, those 20 genes are a sliver of its profile, and the two are probably not the same cell type. Rule 3 removes this one-sided match.

## What the network shows

Numbers from the saved Cytoscape network that was examined:

| Property | Value |
|---|---|
| Nodes | 353 (397 clusters minus 25 Hildreth and 19 spatial) |
| Edges | 1,303 |
| Nodes by class | Immune 120, FAP 116, Vascular 78, Adipocyte 39 |
| Edge weight (shared genes) | 6 to 2,039, median 225 |
| Edges per node | 5 to 32, median 6 |

The network separates into four groups that follow the four classes. No batch correction was applied at this stage. This is the basis for integrating each class separately, and it is the figure the rest of the paper depends on most.

The immune annotations were also checked with CellTypist (Supplementary Fig. 1b). There is no code for that step in the repository.

## Where the evidence stops

- **Edge ranking.** Methods says edges are created from the five highest Jaccard values. The script ranks by average overlap percentage. The two usually agree, but they are different quantities.
- **Marker threshold.** Methods gives log2FC > 0.5 for comparing datasets. This script uses 0. In addition, the filter loop covers 14 datasets. Massier #1, #2 and #4 are not passed through it in the script as published.
- **Which network is the figure?** The Cytoscape session holds several network versions (top 5, top 5 without Hildreth, a second filtered version, top 3). The one examined contains edges whose smaller side is 0.4%, which rule 3 would have removed, and its lowest average overlap is 13.5%. Which version and which thresholds match the published figure cannot be told from the repository.
- **Class labels are manual.** Each cluster was assigned to a class by hand before the network was drawn. The network supports the assignment. It did not produce it.
- **Jaccard depends on list length.** Datasets with few cells give short marker lists, which lowers the index for reasons unrelated to biology. The paper does not separate the two.

[← README](../README.md) · [Next: Fig. 1c–e →](02-fig1cde-integration-benchmark.md)
