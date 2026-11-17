# 06 · Fig. 6: cell–cell communication

[← Fig. 5](05-fig5-adipocytes.md) · [README](../README.md) · [Next: Fig. 7 →](07-fig7-spatial.md)

**Short answer.** "FAPs send, M2-like macrophages receive" was not a starting hypothesis. CellChat scores every combination of sender, receiver and ligand–receptor pair, and the conclusion is read off the totals.

Script: `scripts/Cell_Cell_Communications.R` (CellChat v1.4.0). Subcutaneous fat is analysed in L14–229 and omental fat in L281–499, with identical steps.

## Contents

- [Is every combination explored?](#is-every-combination-explored)
- [How does the conclusion come out of the matrix?](#how-does-the-conclusion-come-out-of-the-matrix)
- [What each panel draws](#what-each-panel-draws)
- [What is concluded](#what-is-concluded)
- [What is computed and what is chosen by hand](#what-is-computed-and-what-is-chosen-by-hand)
- [Where the evidence stops](#where-the-evidence-stops)

## Is every combination explored?

Yes.

```r
WAT.s.sub <- readRDS("WAT_scRNA_MTX_sc_20220921.RDS")
WAT.s.sub <- NormalizeData(WAT.s.sub)

cellchat.sc <- createCellChat(object = WAT.s.sub, assay = "RNA")   # groups = annotated sub-clusters
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat.sc@DB <- CellChatDB.use

cellchat.sc <- subsetData(cellchat.sc)
cellchat.sc <- identifyOverExpressedGenes(cellchat.sc)
cellchat.sc <- identifyOverExpressedInteractions(cellchat.sc)
cellchat.sc <- computeCommunProb(cellchat.sc)
cellchat.sc <- filterCommunication(cellchat.sc, min.cells = 10)
cellchat.sc <- computeCommunProbPathway(cellchat.sc)
cellchat.sc <- aggregateNet(cellchat.sc)
```

| Step | Function | What it does |
|---|---|---|
| 1. Input | `createCellChat` | All annotated sub-clusters of one depot go in together, several dozen groups with names such as `myC15_Mox`. Adipocytes enter as one group, `adipose` |
| 2. Database | `subsetDB` | The human CellChat database of ligand–receptor pairs, restricted to secreted signalling |
| 3. Candidate genes | `identifyOverExpressedGenes`, `identifyOverExpressedInteractions` | Keeps ligands and receptors that are over-expressed in at least one group |
| 4. Exhaustive scoring | `computeCommunProb` | For every sender group, receiver group and ligand–receptor pair, a communication probability is computed from the average ligand expression in the sender and the average receptor expression in the receiver. Group labels are then shuffled to test whether the value is higher than chance |
| 5. Filter | `filterCommunication` | Groups with fewer than 10 cells are dropped |
| 6. Pathway level | `computeCommunProbPathway` | Ligand–receptor pairs that belong to one pathway (for example CCL) are summed |
| 7. Summary | `aggregateNet` | For each sender–receiver pair, the probabilities of all significant ligand–receptor pairs are added up. The result is a sender × receiver matrix, stored in `@net$weight` |

With $N$ groups there are $N \times N$ ordered sender–receiver pairs, a group signalling to itself included, and each is scored for every candidate ligand–receptor pair. Nothing is selected in advance.

## How does the conclusion come out of the matrix?

**Illustrative matrix with three groups.**

| Sender (rows) → receiver (columns) | FAP | M2 | T cell | Outgoing total |
|---|---|---|---|---|
| FAP | 0.5 | 3.0 | 0.5 | **4.0** |
| M2 | 0.3 | 0.4 | 0.3 | 1.0 |
| T cell | 0.1 | 0.6 | 0.1 | 0.8 |
| **Incoming total** | 0.9 | **4.0** | 0.9 | |

- Each row sum is what a group sends. Each column sum is what it receives.
- FAPs send 4.0 and receive 0.9, so they are senders.
- M2 macrophages send 1.0 and receive 4.0, so they are receivers.
- The largest single entry is FAP → M2 (3.0), so the main direction of communication is from FAPs to M2 macrophages.

The real analysis does the same with dozens of groups. Fig. 6a plots the two totals of every group against each other.

## What each panel draws

| Panel | Function | Content |
|---|---|---|
| 6a | `netAnalysis_computeCentrality(slot.name = "netP")`, then `netAnalysis_signalingRole_scatter()` | One dot per sub-cluster, outgoing strength on x and incoming on y, with a diagonal added by `geom_abline(slope = 1)`. FAPs fall below the diagonal, M2-like macrophages above it |
| Supp. 6a | `pheatmap(cellchat.sc@net$weight, clustering_method = "ward.D2", cutree_rows = 3, cutree_cols = 3)` | The whole sender × receiver matrix as a clustered heat map. The FAP rows by M2 columns block stands out. Omental fat uses `ward.D` |
| 6b | `netVisual_circle()` on a sub-matrix | A circle plot for a hand-picked list of sub-clusters. Line width is strength, colour is the sender |
| 6c | `Immune_Cells.R` L326–327 | Stacked violins on the myeloid object for TNF, IL1B, CXCL1/2/3/8 and CCL2/3/4 |
| 6d | `netAnalysis_contribution(cellchat.sc, signaling = pathway)` | Within each of ten pathways, the contribution of every ligand–receptor pair as a bar. The pathways are TNF, KIT, CXCL, CCL, IL16, IL2, CSF, GRN, GALECTIN and PAR. The blue colouring of Mox-driven bars was added by hand afterwards |
| 6e | `netVisual(..., layout = "circle", out.format = "pdf")` | Circle plots per pathway and per ligand–receptor pair |
| Supp. 6b | `computeNetSimilarity(type = "functional")`, `netEmbedding`, `netClustering`, `netVisual_embedding` | Pathways grouped by functional similarity, that is by having similar senders and receivers |

The sub-clusters drawn in 6b for subcutaneous fat are written out in the script (L132–146): sfC16, adipose, sfC06, sfC10, sfC11, sfC15, myC02 (LAM), myC15 (Mox), myC14, vC11, myC00, myC01 and myC13.

## What is concluded

1. **FAPs relay, M2-like macrophages receive.** In both depots, FAPs mainly send signals, and those signals are mainly received by M2-like macrophages. Pathways that link the two include complement, chemerin and IL16 (Supplementary Fig. 6b).
2. **Mox macrophages are a hub in subcutaneous fat.** Mox (myC15) occur only in subcutaneous fat and express very high levels of pro-inflammatory cytokines and chemokines (6c). They both send and receive CCL, CXCL and TNF signals:

| Direction | Partner | Ligand → receptor |
|---|---|---|
| Mox receive | CD8 T cells (lyC02, lyC03) | CCL5 → CCR1 |
| Mox send | Venous endothelial cells (vC01, vC11) | CCL2, CXCL2, CXCL8 → ACKR1 |
| Mox send | Adipocytes and myeloid cells | TNF → TNFRSF1A |
| Mox send | CD8 T and NK cells (lyC03, lyC06, lyC09), endothelial cells (vC01–03) | TNF → TNFRSF1B |

## What is computed and what is chosen by hand

| Computed | Chosen by hand |
|---|---|
| Every communication probability | The dozen or so sub-clusters drawn in 6b |
| The outgoing and incoming totals of each group | The ten pathways shown in 6d |
| The clustering of the matrix in Supp. 6a | The focus on Mox, and the blue highlight in 6d |

## Where the evidence stops

- **Inference from expression.** A ligand high in one group and its receptor high in another counts as communication. Whether the two cell types sit next to each other in the tissue is not considered, and there is no protein-level evidence.
- **Secreted signalling only.** Extracellular-matrix and cell-contact pairs are excluded by the `subsetDB` call. Methods says only that CellChatDB was used. FAPs produce much of the extracellular matrix, so this restriction is not neutral for them.
- **Depots are not comparable.** The axes of 6a go up to 9.1 for subcutaneous and 3.7 for omental fat. Strengths depend on the number and size of groups in each run and cannot be compared between depots.
- **The authors' check.** Two known facts are reproduced: the KIT receptor appears only in mast cells, and periostin is sent only by pericytes (Supplementary Fig. 6c). That shows the pipeline is not broken. It does not prove the new predictions.
- **Mox come from one depot.** The script has a separate block (L329–331 of `Immune_Cells.R`) that plots Mox by dataset, to check that they come from more than one cohort. That plot is not in the main figures.
- **A general property of the method.** Cell types that express many secreted ligands tend to come out as senders, and cell types with many receptors as receivers. This is my caution, not a statement from the paper.

[← Fig. 5](05-fig5-adipocytes.md) · [README](../README.md) · [Next: Fig. 7 →](07-fig7-spatial.md)
