# 05 · Fig. 5: why adipocytes resist subtyping

[← Fig. 4](04-fig4-fap-pseudotime.md) · [README](../README.md) · [Next: Fig. 6 →](06-fig6-cell-communication.md)

**Short answer.** Fig. 5 explains a failure. Adipocyte sub-clusters from different single-nucleus studies do not match, the cause points to the sequencing platform, and so all adipocytes were pooled into one class.

None of the code for this figure is in the repository. Everything below follows the main text, Methods and the figure legend.

## Contents

- [What went wrong?](#what-went-wrong)
- [What does each panel show?](#what-does-each-panel-show)
- [How the argument runs](#how-the-argument-runs)
- [What follows from it](#what-follows-from-it)
- [Where the evidence stops](#where-the-evidence-stops)

## What went wrong?

Immune, vascular and FAP cells all gave stable sub-clusters after integration. Adipocytes did not.

- Only snSeq captures mature adipocytes, because the cells are too large and fragile for droplet-based single-cell sequencing. The class therefore comes from seven datasets: 4 subcutaneous, 2 omental, 1 perivascular.
- Integrated together, adipocytes scored poorly on ARI, kBET and LISI (see the result table in [chapter 02](02-fig1cde-integration-benchmark.md#result)).
- Integrating each depot separately did not help. Clusters still followed the study they came from and not shared biology (Supplementary Fig. 5a–b).

A possible reason for incomplete integration is that the studies share few common features. To test this, the authors stepped back and analysed the studies one by one.

## What does each panel show?

| Panel | What was done | What it shows |
|---|---|---|
| 5a, left | In the two largest subcutaneous snSeq datasets (Emont and Massier #1), the top 50 markers of every cluster were compared across the four classes | Adipocyte sub-clusters have lower fold changes than the other classes. They differ less from each other |
| 5a, right | Cell-type labels were transferred from one dataset to the other, and marker overlap was measured | Overlap is high for the other classes and low for adipocyte sub-clusters. The sub-clusters of the two studies do not match |
| 5b | Adipocyte markers that do replicate between the two studies | A short list of lipid-metabolism genes: ABCD2, ACACB, CD36, DGAT2, GPAM, HACD2, LPL |
| 5c, top | Pseudobulk profiles of the snSeq datasets and of the spatial data were clustered together with bulk RNA-seq of isolated mature adipocytes | The spatial data sit close to bulk. All snSeq datasets form their own branch |
| 5c, bottom | The same with snSeq, spatial and bulk data from one person (Massier #3) | The same picture |
| 5d | 218 transcripts enriched at least 50-fold in adipocytes against all other samples of the FANTOM5 atlas were compared across platforms as z-scores | The same adipocyte genes are expressed very differently depending on platform |
| 5e | Genes whose z-score differs by more than 5 between snSeq and spatial data | LEP, PLIN4, SAA1 and RBP4 are high in spatial and bulk data and low or absent in snSeq. PPARG, WDPCP and PDE3B are higher in snSeq |

**Illustrative reading of 5e.** If LEP has a z-score of 6.2 in the spatial data and 0.4 in the integrated snSeq data, the difference of 5.8 exceeds 5 and the gene is flagged as platform-discordant.

Methods says only that each dataset was z-score normalised. I read this as standardising within each dataset.

## How the argument runs

1. **Adipocyte sub-clusters differ little to begin with** (5a, left).
2. **The sub-clusters found by different studies do not match** (5a right, 5b).
3. **The cause points to the platform.** The snSeq profile of adipocytes is further from isolated mature adipocytes than the spatial profile is (5c).
4. **The discordant genes are important ones.** They include the very markers that an earlier spatial study (Bäckdahl et al.) used to define three adipocyte subtypes, named after LEP, PLIN and SAA (5d, 5e).

Point 4 closes the loop. The markers that define the known adipocyte subtypes are measured poorly by snSeq, so snSeq cannot recover those subtypes.

The mechanism offered is taken from earlier literature: snSeq captures nuclear RNA and is biased towards nascent and long transcripts. The paper's own evidence is the comparison with bulk data in 5c.

## What follows from it

- In all later analyses, adipocytes from snSeq are **one pooled class**. In the communication analysis they enter as a single group called `adipose`, and in the spatial reference too.
- In Fig. 8, adipocyte subtypes are scored with the **spatial markers**, not with snSeq markers (Supplementary Fig. 8a).
- The authors' conclusion is modest. Platform-specific biases "may contribute" to the low reproducibility, and several platforms are needed to define adipocyte subtypes reliably.

## Where the evidence stops

- **Bias, not absence.** The figure supports "snSeq is biased for adipocytes". It does not show that adipocytes have no subtypes, and the authors do not claim that.
- **Bulk as the reference.** Bulk RNA-seq of isolated adipocytes is treated as the truth. The paper does not discuss contamination by other cell types or expression changes caused by the isolation procedure.
- **Mixtures resemble mixtures.** A spatial spot and a bulk sample both average many cells. That shared property may explain part of their similarity. The paper does not rule this out.
- **Whole cell against nucleus.** Bulk and spatial data measure whole-cell RNA, and snSeq measures nuclear RNA. Part of the difference is expected for any cell type. The paper does not show the same comparison for a class that integrates well, such as FAPs, which would tell how much of the gap is specific to adipocytes.
- **Methods gaps.** The method used for label transfer in 5a, and the distance and gene set behind the dendrogram in 5c, are not stated.

For what it teaches about method, this is the most useful figure in the paper. It shows how to use bulk and spatial data as outside references to judge whether heterogeneity seen in single-nucleus data can be trusted.

[← Fig. 4](04-fig4-fap-pseudotime.md) · [README](../README.md) · [Next: Fig. 6 →](06-fig6-cell-communication.md)
