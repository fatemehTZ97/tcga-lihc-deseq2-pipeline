# TCGA-LIHC Apoptosis/Inflammation Candidate Gene Panel

## What this pipeline does

Differential expression analysis of TCGA-LIHC (liver hepatocellular
carcinoma) tumor vs. TCGA's own adjacent-normal liver tissue, focused on a
literature-selected panel of 9 apoptosis and inflammation genes (IL6, TNF,
IL10, BAX, BCL2, CASP3, CASP9, TLR4, TLR2). The panel was chosen a priori
based on mechanistic relevance to a wet-lab study of *Akkermansia
muciniphila* effects on HepG2 hepatocellular carcinoma cells, and this
pipeline evaluates the clinical/human-tumor relevance of that panel before
functional validation.

## Data source

- **Cohort:** TCGA-LIHC (Liver Hepatocellular Carcinoma), accessed via the
  [GDC Data Portal](https://portal.gdc.cancer.gov/) using
  [`TCGAbiolinks`](https://bioconductor.org/packages/TCGAbiolinks/) (v2.32.0).
- **Data type:** Gene Expression Quantification, STAR - Counts workflow
  (raw counts), reference genome hg38.
- **Samples:** 424 total — 371 Primary Tumor, 3 Recurrent Tumor, 50 Solid
  Tissue Normal (adjacent non-tumor liver tissue from the same patients).
  Primary and Recurrent Tumor samples were pooled into a single Tumor group
  (n = 374) for all analyses.
- **Retrieved:** 2026 (see `sessionInfo()` output in script for exact
  package/Bioconductor versions used).

## Method summary

1. **Data access & QC** — TCGA-LIHC raw counts and sample metadata retrieved
   via `TCGAbiolinks::GDCquery`/`GDCdownload`/`GDCprepare`. Sample counts
   cross-checked at three independent stages (project summary, query
   results, prepared `SummarizedExperiment`) before any analysis.
2. **Sample/expression alignment** — TCGA barcodes retained as row/column
   names throughout (never renamed by position), with `identical()` checks
   confirming expression-matrix and metadata alignment at every step.
3. **Differential expression** — [`DESeq2`](https://bioconductor.org/packages/DESeq2/)
   (v1.44.0), binary Tumor-vs-Normal design (`~ sample_type_binary`,
   `Normal` as reference level), default independent filtering and
   Benjamini-Hochberg FDR correction.
4. **Validation** — Every reported gene's direction (higher/lower in tumor)
   was independently confirmed by comparing raw group means of
   DESeq2-normalized counts against the sign of `log2FoldChange`, before
   being used in any figure.
5. **Visualization** — Heatmap (`pheatmap`, VST-transformed, per-gene
   z-scored, samples and genes explicitly ordered and alignment-checked
   before plotting), volcano plot (`EnhancedVolcano`, genome-wide with
   candidate panel highlighted), PCA (`DESeq2::plotPCA`).
6. **Functional enrichment** — GO (Biological Process) and KEGG pathway
   enrichment (`clusterProfiler`) on the full set of significant DEGs
   (padj < 0.05, |log2FC| > 1), as a genome-wide characterization —
   **not** used to select the candidate gene panel, which was defined a
   priori from the literature.

## Key findings

- **7 of 9 candidate genes were significantly dysregulated** (padj < 0.05):
  IL6, TNF, IL10, TLR4, TLR2, BAX, CASP3. BCL2 and CASP9 trended in the
  expected direction but did not reach significance.
- **All 5 inflammation/TLR genes (IL6, TNF, IL10, TLR4, TLR2) were *lower*
  in tumor than in TCGA's adjacent-normal tissue.** This was initially
  unexpected (it runs counter to some external references, e.g. GEPIA-style
  tumor-vs-pooled-normal comparisons), but was consistent and internally
  validated across every gene checked. The most likely explanation, and the
  central limitation of this analysis, is discussed below.
- **BAX and CASP3 (pro-apoptotic) were significantly higher in tumor**,
  consistent with expected HCC biology.
- **KEGG "Apoptosis" (hsa04210) and "Toll-like receptor signaling pathway"
  were NOT significantly enriched genome-wide** (padj = 1.0), despite
  individual panel genes in those pathways being strongly significant.
  This likely reflects reduced statistical power for pathway-level
  enrichment when a very large fraction of the transcriptome (~30% here) is
  already flagged as differentially expressed, diluting any pathway-specific
  enrichment signal — not necessarily an absence of real pathway
  involvement.

## How to run it

1. Install R (this pipeline was developed and tested on R 4.4.2) and the
   packages listed in the `SETUP` section at the top of
   `tcga_lihc_deseq2_pipeline.R`.
2. Set your working directory (edit the `setwd()` line near the top of the
   script, or open the script inside an RStudio Project so relative paths
   work automatically).
3. Run the script top to bottom, in order — each section depends on objects
   created in the previous one (there is no separate script-splitting; all
   steps are in a single file, organized into numbered sections with
   checkpoint comments).
4. **Expected runtime:** the initial GDC download (~424 files, ~1.8 GB) is
   the slowest step and can take anywhere from 10 minutes to significantly
   longer depending on connection speed; it only needs to run once; on
   subsequent runs `GDCdownload()` will detect already-downloaded files.
   `DESeq()` itself takes a few minutes on this sample size.
5. **Outputs:** three PNG figures are written to the working directory —
   `heatmap_9gene_panel_grouped.png`, `volcano_9gene_highlighted.png`,
   `pca_tumor_normal.png` (plus `pca_detailed_sampletype.png`).

No external input files are required beyond what `TCGAbiolinks` downloads
directly from GDC.

## Known limitations

- **TCGA-LIHC's "Solid Tissue Normal" is adjacent tissue, not a healthy
  baseline.** The great majority of TCGA-LIHC patients have underlying
  HBV/HCV-driven cirrhosis, and adjacent non-tumor liver tissue in this
  cohort is frequently itself chronically inflamed/cirrhotic. This is the
  most plausible explanation for why inflammatory cytokines (IL6, TNF, IL10)
  and TLR2/TLR4 read *lower* in tumor than in "normal" here — the
  comparison is better understood as *tumor vs. chronically inflamed
  adjacent tissue*, not *tumor vs. healthy liver*. A planned follow-up
  pipeline incorporates GTEx healthy-donor liver (via the UCSC Xena Toil
  harmonized TCGA/GTEx recompute) to establish a true healthy baseline and
  resolve this ambiguity.
- **Ensembl gene ID version suffixes are resolved dynamically** (via pattern
  matching on the unversioned ID) rather than hardcoded, but the underlying
  gene set/annotation could still shift slightly across future GDC data
  releases or `org.Hs.eg.db` versions.
- **Sample size imbalance** (Tumor n = 374 vs. Normal n = 50) gives high
  statistical power to detect even small-magnitude expression differences,
  which is part of why a large fraction (~30%) of the tested transcriptome
  reaches statistical significance — significance here should be
  interpreted alongside effect size (log2FoldChange), not in isolation.
- **Recurrent Tumor samples (n = 3) are pooled with Primary Tumor.** This
  is a deliberate, documented simplification for statistical power; these
  3 samples were not analyzed as a separate group.
- **GO/KEGG enrichment universe correction is applied**, but pathway
  databases (KEGG, GO) are themselves imperfect and only partially
  overlapping — the "not significantly enriched" finding for the apoptosis
  pathway should be read as a property of this specific enrichment test on
  this specific gene list, not as evidence that apoptosis is uninvolved in
  HCC (which is well-established in the broader literature).
