# Spider Embryo scRNA-seq — Drop-seq Pipeline & Cell Atlas

Single-cell RNA-seq pipeline and analysis for *Parasteatoda tepidariorum* (common house spider) embryos, covering four developmental-stage samples from raw reads through to an integrated, clustered cell atlas with marker genes, plus a per-sample isoform usage analysis.

## Background

*P. tepidariorum* is a key chelicerate model for studying arthropod body-plan evolution and segmentation, thanks to its accessible embryology and an established genome/gene-annotation resource. Chelicerates diverged from the insect/crustacean lineage early enough that comparing their developmental gene expression programs to *Drosophila* and other arthropods helps resolve which patterning mechanisms (e.g. segmentation clock genes, appendage identity genes) are ancestral versus lineage-specific.

Single-cell RNA-seq lets us move beyond bulk, whole-embryo expression profiles to ask which of these patterning genes are expressed in which cell populations — germ layers, segment addition zone, appendage primordia, nervous system, hemocytes, etc. — at a given developmental stage. Layering an isoform-level view on top of that lets us ask a further question: whether individual genes switch which transcript isoform they predominantly express as development proceeds, independent of overall expression level. This repository reprocesses raw Drop-seq reads from four embryonic-stage samples into a labeled, integrated cell atlas, cross-checks the resulting cluster identities against a previously published cell-type labeling, and separately screens for genes with sample-to-sample isoform switching.

## Methods

### scRNA-seq pipeline

The pipeline (see `Snakefile`) follows the standard Drop-seq-tools workflow:

1. **QC** (`fastqc`) — raw read quality per sample/mate.
2. **Barcode tagging** — cell barcode (XC) and UMI (XM) are extracted from read 2 and attached to read 1 as BAM tags (`TagBamWithReadSequenceExtended`), followed by quality-based read filtering (`FilterBam`).
3. **Adapter/poly-A trimming** — template-switch oligo and poly-A tail sequences are trimmed from the cDNA read (`TrimStartingSequence`, `PolyATrimmer`).
4. **Alignment** — trimmed cDNA reads are aligned to the *P. tepidariorum* genome with STAR, then merged back with the original barcode/UMI tags (`MergeBamAlignment`).
5. **Gene tagging & quantification** — reads are annotated with overlapping gene/exon/intron info (`TagReadWithGeneFunction`) and collapsed into a cell-by-gene digital expression matrix (`DigitalExpression`), counting coding + UTR + intronic reads (to capture nascent/unspliced signal).
6. **Seurat object construction** — one raw `Seurat` object per sample is built directly from each DGE matrix.
7. **Per-sample QC + doublet removal** — cells are filtered on feature count, UMI count, and percent-mitochondrial content (stage-specific UMI floors where warranted), then putative doublets are removed with `DoubletFinder` (parameter-swept `pK`, `SCTransform`-normalized counts).
8. **Integration & clustering** — the four QC'd, SCT-normalized samples are integrated with Seurat's anchor-based `IntegrateData` (reciprocal PCA), then clustered (Louvain, shared-nearest-neighbor graph) and visualized with UMAP.
9. **Marker gene identification** — `FindAllMarkers` (Wilcoxon rank-sum) identifies genes enriched in each cluster relative to all others.

### Isoform usage analysis

Per-sample, per-transcript read counts (from the genome annotation's alternative isoforms) were converted to within-gene isoform proportions for each of the four samples, then screened for genes whose isoform proportions vary substantially across samples (`isoform_switching_summary.tsv`: max range in isoform proportion across samples, and the standard deviation of that proportion across samples, per gene).

## Results

### QC and doublet removal

| Sample | Raw cells (barcodes) | After QC filtering | After doublet removal | QC pass rate | Doublet rate (of QC-passed) |
|---|---:|---:|---:|---:|---:|
| ERR9871746 | 12,489 | 6,650 | 6,099 | 53.2% | 8.3% |
| ERR9871747 | 12,089 | 3,071 | 2,920 | 25.4% | 4.9% |
| ERR9871748 | 15,865 | 10,215 | 9,125 | 64.4% | 10.7% |
| ERR9871749 | 15,236 | 5,198 | 4,864 | 34.1% | 6.4% |
| **Total** | **55,679** | **25,134** | **23,008** | — | — |

QC-pass rates vary considerably by sample (25–64%), consistent with differences in raw library quality/ambient RNA content rather than a fixed cutoff issue:

![ERR9871746 QC](results/qc/ERR9871746_raw_qc_violin.png)
![ERR9871747 QC](results/qc/ERR9871747_raw_qc_violin.png)
![ERR9871748 QC](results/qc/ERR9871748_raw_qc_violin.png)
![ERR9871749 QC](results/qc/ERR9871749_raw_qc_violin.png)

### Clustering

Integration of the four QC'd samples yielded **22 clusters** across 23,008 cells.

![Integrated UMAP by cluster](results/clustering/integrated_umap_clusters.png)

### Cluster identity

Comparing this re-run's clustering against a previously published cell-type labeling for the same clusters shows the same 22-cluster topology, labeled with recognizable tissue/cell-type categories: segment addition zone (SAZ, clusters 2/7), mesoderm (3/9), leg (4/6), stripe/segment-mediating zone (SMZ, 5), peripheral nervous system (PNS, 8), posterior compartment (10), endoderm (11), pedipalp (12/18), hemocytes (13/17/20), central nervous system (CNS, 14), and posterior/anterior "Pc" populations (15/19, 16/21).

![Cluster identities — this re-run](results/clustering/UMAP_Clusters.png)
![Cluster identities — published labeling](results/clustering/umap_paper_labels.png)

Two clusters left as **"Undetermined"** in the published labeling (clusters 0 and 1) resolve to specific marker-gene identities in this re-run: cluster 0's top differential markers align with **netrin-1**, and cluster 1's with an **ecdysone receptor**-associated gene — both plausible given their roles in neural guidance and molting/developmental timing, respectively, though this identification is based on marker genes alone and hasn't been validated further (e.g. in situ).

### Marker genes

`FindAllMarkers` returned **1,809 significant marker genes** (adjusted p ≤ 0.05) across the 22 clusters (full table: `results/markers/all_markers.tsv`), ranging from 22 markers (cluster 21) to 224 markers (cluster 20, hemocytes) per cluster. Top markers per cluster are strongly cluster-specific — e.g. cluster 7 (SAZ) and cluster 20 (hemocytes) both show marker genes with log2FC > 7 and near-exclusive detection (`pct.1` > 0.8, `pct.2` < 0.02), indicating clean separation rather than a continuum.

### Isoform usage across developmental stages

Of **21,593 genes** with isoform-level counts, **6,523 (30.2%)** have two or more annotated transcripts. Screening these for cross-sample variation in isoform proportion (`results/isoforms/isoform_switching_summary.tsv`) identified **2,656 genes** with detectable isoform-usage variability across the four samples, of which **1,143** show a proportion range greater than 0.9 (near-complete switching) and **1,100** show a full switch (range ≈ 1.0 — an isoform going from 100% to 0% usage, or vice versa, between at least two samples).

Three representative examples (proportions from `results/isoforms/variable_isoforms.tsv`):

- **LOC107439500** — a clean two-isoform switch between the two latest-stage samples: XM_043054336.2 goes from ~50–57% (ERR9871746/47) to 100% in ERR9871748, then to 0% in ERR9871749, where XM_043054335.2 (0% in ERR9871748) becomes exclusive.
- **LOC107436484** — XM_043040168.2 is the exclusive isoform in ERR9871746/47 (100%), drops to 19.5% in ERR9871748, and disappears entirely (0%) in ERR9871749, where two additional isoforms (XM_071179417.1, XM_071179418.1, both entirely absent in the first two samples) each account for 50%.
- **LOC107457063** — a three-way switch: XM_071180825.1 dominates ERR9871746 (100%) and largely disappears elsewhere; XM_043044375.2 dominates ERR9871747 only (96.1%); XM_043044377.2 dominates both ERR9871748 and ERR9871749 (81–82%) but is entirely absent from the first two samples.

![LOC107436484](results/isoforms/isoform_LOC107436484.png)
![LOC107438554](results/isoforms/isoform_LOC107438554.png)
![LOC107439500](results/isoforms/isoform_LOC107439500.png)
![LOC107442251](results/isoforms/isoform_LOC107442251.png)
![LOC107442412](results/isoforms/isoform_LOC107442412.png)
![LOC107446349](results/isoforms/isoform_LOC107446349.png)
![LOC107449306](results/isoforms/isoform_LOC107449306.png)
![LOC107452495](results/isoforms/isoform_LOC107452495.png)
![LOC107455433](results/isoforms/isoform_LOC107455433.png)
![LOC107457063](results/isoforms/isoform_LOC107457063.png)

These are per-sample, whole-embryo isoform proportions rather than cluster-resolved — a natural extension would be checking whether any of these switches track a specific cell population from the atlas above rather than the whole embryo uniformly.

### Repository layout for results

```
results/
  qc/
    ERR9871746_raw_qc_violin.png   # per-sample nFeature/nCount distributions
    ERR9871747_raw_qc_violin.png
    ERR9871748_raw_qc_violin.png
    ERR9871749_raw_qc_violin.png
    qc_summary.tsv                 # cell counts before/after QC/doublet removal
  clustering/
    UMAP_Clusters.png              # this re-run's cluster labels
    umap_paper_labels.png          # published reference labeling, same clusters
    integrated_umap_clusters.png   # integrated UMAP colored by cluster number
  markers/
    all_markers.tsv                # FindAllMarkers output, all 22 clusters
  isoforms/
    isoform_LOC*.png               # per-gene isoform usage bar plots (10 examples)
    isoform_switching_summary.tsv  # per-gene proportion range/SD across samples
    transcripts_per_gene.tsv       # annotated transcript count per gene
    variable_isoforms.tsv          # per-transcript proportions for variable genes
```
