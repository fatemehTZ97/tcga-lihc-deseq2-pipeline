## ============================================================================
## TCGA-LIHC RNA-seq Pipeline: Apoptosis/Inflammation Candidate Gene Panel
## ============================================================================
## Purpose: Differential expression analysis of TCGA-LIHC tumor vs. TCGA's own
## adjacent-normal liver tissue, focused on a literature/mechanism-driven
## candidate gene panel (IL6, TNF, IL10, BAX, BCL2, CASP3, CASP9, TLR4, TLR2)
## selected to support a wet-lab study on Akkermansia muciniphila effects on
## HepG2 hepatocellular carcinoma cells.
##
## IMPORTANT SCOPE NOTE: this pipeline compares TCGA tumor vs. TCGA's own
## "Solid Tissue Normal" samples only. As documented in the README, this
## adjacent-normal tissue is frequently cirrhotic/chronically inflamed rather
## than truly healthy liver -- a key limitation discussed in the README and
## addressed in a separate, planned TCGA-vs-GTEx follow-up pipeline.
## ============================================================================


## ----------------------------------------------------------------------------
## 0. SETUP: packages and environment
## ----------------------------------------------------------------------------

# One-time installs (uncomment if not already installed)
# devtools::install_github("kevinblighe/EnhancedVolcano")
# install.packages(c("ggplot2", "ggrepel", "msigdbr"))
# BiocManager::install(c("DESeq2", "TCGAbiolinks", "clusterProfiler",
#                         "org.Hs.eg.db", "pheatmap"))

library(ggplot2)
library(ggrepel)
library(DESeq2)
library(SummarizedExperiment)
library(tidyverse)
library(TCGAbiolinks)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Record package versions for reproducibility -- always worth capturing
# alongside your results, since TCGAbiolinks/DESeq2 APIs have changed
# across versions in the past.
sessionInfo()

# NOTE: set this to your own local working directory before running.
# The original analysis used a hardcoded Windows path (not portable);
# replace the line below with your own path, or better, use an RStudio
# Project (.Rproj) and a relative path so the script runs on any machine.
setwd(".")  # <-- CHANGE THIS to your working directory


## ----------------------------------------------------------------------------
## 1. DATA ACCESS: TCGA-LIHC RNA-seq via GDC / TCGAbiolinks
## ----------------------------------------------------------------------------

# Full list of GDC projects (reference only, not used downstream)
gdcprojects <- getGDCprojects()

# Project-level summary for TCGA-LIHC specifically -- used as an independent
# baseline case count to cross-check against the query results below.
projectsummary <- getProjectSummary('TCGA-LIHC')
projectsummary$data_categories
# Expect ~376-377 cases with "Transcriptome Profiling" data available.

# Build the query: STAR-Counts gene expression quantification for TCGA-LIHC.
# This only queries the GDC API -- no files are downloaded at this step.
queryLIHC <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# CHECKPOINT: confirm total sample count and tumor/normal breakdown BEFORE
# downloading. TCGA-LIHC is expected to have ~371 Primary Tumor,
# a small number (~3) Recurrent Tumor, and ~50 Solid Tissue Normal samples.
nrow(queryLIHC$results[[1]])
table(queryLIHC$results[[1]]$sample_type)

# Download the actual per-sample files (only needs to run once; re-running
# will just verify checksums of already-downloaded files).
GDCdownload(queryLIHC, method = 'api', files.per.chunk = 20)

# Assemble all downloaded files into a single SummarizedExperiment object
# (genes x samples), with sample metadata (incl. sample_type) in colData().
TCGAData <- GDCprepare(queryLIHC, summarizedExperiment = TRUE)

# CHECKPOINT: confirm dimensions and sample_type breakdown match the query
# above -- this is the first place a silent GDCprepare data-loss issue
# would show up.
dim(TCGAData)
table(colData(TCGAData)$sample_type)
assayNames(TCGAData)  # confirm "unstranded" (raw counts) is present


## ----------------------------------------------------------------------------
## 2. BUILD EXPRESSION MATRIX AND SAMPLE METADATA (barcode-matched throughout)
## ----------------------------------------------------------------------------
## Design choice: TCGA barcodes are kept as row/column names throughout this
## entire pipeline (never renamed to generic labels like "Tumor1"/"Normal1").
## This is deliberate -- renaming samples by position is a common source of
## sample/annotation mismatch bugs, and keeping real barcodes lets every
## alignment be checked explicitly with identical() rather than assumed.

# Raw counts matrix (required by DESeq2; NOT normalized TPM/FPKM)
expr_mat <- assay(TCGAData, "unstranded")

# Sample metadata, indexed by the same real TCGA barcodes
sample_info <- as.data.frame(colData(TCGAData))
sample_type <- factor(sample_info$sample_type)
names(sample_type) <- rownames(sample_info)

# CHECKPOINT: expression matrix columns and sample_type names must match
# exactly, in the same order.
identical(colnames(expr_mat), names(sample_type))
dim(expr_mat)
length(sample_type)

# Build coldata for DESeq2, again keeping real barcodes as rownames
coldata <- data.frame(
  sample_type = sample_type,
  row.names = names(sample_type)
)
table(coldata$sample_type)
identical(colnames(expr_mat), rownames(coldata))  # must be TRUE


## ----------------------------------------------------------------------------
## 3. DEFINE TUMOR VS. NORMAL GROUPING AND RUN DESeq2
## ----------------------------------------------------------------------------
## Design decision: Primary Tumor and Recurrent Tumor samples are pooled into
## a single "Tumor" group (n=374) and compared against "Solid Tissue Normal"
## (n=50). This is deliberately a BINARY factor, not the raw 3-level
## sample_type -- using a multi-level factor with an unspecified contrast is
## a known pitfall with DESeq2's results(): results() silently returns only
## the LAST coefficient in the model by default, which for a 3-level factor
## means comparing the smallest group (here, n=3 Recurrent Tumor) against the
## reference, NOT the intended pooled Tumor-vs-Normal comparison. Collapsing
## to a binary factor up front avoids this ambiguity entirely.

coldata$sample_type_binary <- ifelse(
  coldata$sample_type %in% c("Primary Tumor", "Recurrent Tumor"),
  "Tumor",
  "Normal"
)
coldata$sample_type_binary <- relevel(factor(coldata$sample_type_binary), ref = "Normal")

table(coldata$sample_type_binary)   # expect Normal=50, Tumor=374
levels(coldata$sample_type_binary)  # "Normal" must be first/reference level,
                                     # so that positive log2FC means "higher in Tumor"

# STAR counts can occasionally be non-integer due to internal GDC processing;
# DESeq2 requires integer counts.
expr_mat_int <- round(expr_mat)

dds <- DESeqDataSetFromMatrix(
  countData = expr_mat_int,
  colData = coldata,
  design = ~ sample_type_binary
)

# Filter out genes with negligible expression across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dim(dds)  # expect ~49,000 genes retained out of ~60,660

dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Confirm which contrast this actually represents -- should read
# "sample_type_binary_Tumor_vs_Normal", not an ambiguous multi-level contrast.
resultsNames(dds)

norm_counts <- counts(dds, normalized = TRUE)


## ----------------------------------------------------------------------------
## 4. VALIDATE: manual mean-comparison vs. DESeq2 log2FoldChange sign
## ----------------------------------------------------------------------------
## For every gene reported below, we independently compute the mean
## normalized expression in Tumor vs. Normal and confirm its direction
## agrees with the SIGN of DESeq2's log2FoldChange. This check is run for
## every gene in the candidate panel before any figure is built -- it is
## what allows us to distinguish genuine biological findings from pipeline
## bugs (e.g., a mis-specified factor level or column/row mismatch).

candidate_genes <- c(
  IL6   = "ENSG00000136244", TNF  = "ENSG00000232810", IL10  = "ENSG00000136634",
  BAX   = "ENSG00000087088", BCL2 = "ENSG00000171791", CASP3 = "ENSG00000164305",
  CASP9 = "ENSG00000132906", TLR4 = "ENSG00000136869", TLR2  = "ENSG00000137462"
)

# Ensembl IDs in the count matrix carry version suffixes (e.g. ".12") that can
# change between GDC data releases -- so we look up the current suffix by
# pattern match rather than hardcoding it, for a small amount of future-
# proofing (see README "Known limitations").
candidate_ids <- sapply(candidate_genes, function(x) grep(paste0("^", x), rownames(dds), value = TRUE))
candidate_ids

for (gene in names(candidate_ids)) {
  id <- candidate_ids[[gene]]
  if (length(id) == 0) { cat(gene, ": NOT FOUND in filtered dds\n\n"); next }
  mt <- mean(norm_counts[id, coldata$sample_type_binary == "Tumor"])
  mn <- mean(norm_counts[id, coldata$sample_type_binary == "Normal"])
  lfc <- res[id, "log2FoldChange"]
  padj <- res[id, "padj"]
  cat(sprintf("%s (%s)\n  mean_tumor=%.2f  mean_normal=%.2f  higher_in_tumor=%s  log2FC=%.3f  padj=%.3g\n\n",
              gene, id, mt, mn, mt > mn, lfc, padj))
}

## RESULT SUMMARY (see README for full discussion):
## 7 of 9 candidate genes are significant at padj < 0.05.
## All 5 inflammation/TLR genes (IL6, TNF, IL10, TLR4, TLR2) are LOWER in
## tumor than in TCGA's adjacent-normal tissue -- consistent with TCGA-LIHC's
## "Solid Tissue Normal" representing chronically inflamed/cirrhotic adjacent
## tissue rather than a healthy baseline (see README limitations).
## BAX and CASP3 (pro-apoptotic) are significantly HIGHER in tumor, as
## expected from HCC biology. BCL2 and CASP9 trend in the expected
## anti-apoptotic direction but do not reach significance.


## ----------------------------------------------------------------------------
## 5. [OPTIONAL / EXPLORATORY] Extended receptor and pathway-partner checks
## ----------------------------------------------------------------------------
## These checks were used during development to inform gene-panel scoping
## decisions (e.g., whether to include additional TLR-family receptors, and
## whether nearby apoptosis-pathway genes showed stronger effect sizes than
## the core panel). They are NOT part of the final 9-gene candidate panel or
## any published figure, but are kept here for transparency into the
## decision-making process. Safe to skip/remove if not of interest.

# Other pattern-recognition receptors related to TLR2/TLR4 signaling
prr_genes <- c(
  TLR1 = "ENSG00000174125", TLR3 = "ENSG00000164342", TLR5 = "ENSG00000187554",
  TLR6 = "ENSG00000174130", TLR7 = "ENSG00000196664", TLR8 = "ENSG00000101916",
  TLR9 = "ENSG00000239887", NOD1 = "ENSG00000106100", NOD2 = "ENSG00000167207",
  NLRP3 = "ENSG00000162711", MYD88 = "ENSG00000172936", AGER = "ENSG00000204305",
  AHR = "ENSG00000106546"
)
prr_ids <- sapply(prr_genes, function(x) grep(paste0("^", x), rownames(dds), value = TRUE))
for (gene in names(prr_ids)) {
  id <- prr_ids[[gene]]
  if (length(id) == 0) { cat(gene, ": NOT FOUND in filtered dds\n\n"); next }
  mt <- mean(norm_counts[id, coldata$sample_type_binary == "Tumor"])
  mn <- mean(norm_counts[id, coldata$sample_type_binary == "Normal"])
  lfc <- res[id, "log2FoldChange"]; padj <- res[id, "padj"]
  cat(sprintf("%s (%s)\n  mean_tumor=%.2f  mean_normal=%.2f  higher_in_tumor=%s  log2FC=%.3f  padj=%.3g\n\n",
              gene, id, mt, mn, mt > mn, lfc, padj))
}

# Direct pathway partners of the core panel (intrinsic apoptosis cascade;
# TLR/NF-kB signaling cascade)
pathway_partners <- c(
  APAF1 = "ENSG00000120868", CYCS = "ENSG00000172115", BAK1 = "ENSG00000030110",
  BID = "ENSG00000015475", BAD = "ENSG00000002330", MCL1 = "ENSG00000143384",
  XIAP = "ENSG00000101966", DIABLO = "ENSG00000184047", CASP7 = "ENSG00000165806",
  TRAF6 = "ENSG00000175104", IRAK1 = "ENSG00000184216", IRAK4 = "ENSG00000198001",
  NFKB1 = "ENSG00000109320", RELA = "ENSG00000173039", IKBKB = "ENSG00000104365",
  TNFAIP3 = "ENSG00000118503", IL1R1 = "ENSG00000115594", CD14 = "ENSG00000170458"
)
partner_ids <- sapply(pathway_partners, function(x) grep(paste0("^", x), rownames(dds), value = TRUE))
for (gene in names(partner_ids)) {
  id <- partner_ids[[gene]]
  if (length(id) == 0) { cat(gene, ": NOT FOUND in filtered dds\n\n"); next }
  mt <- mean(norm_counts[id, coldata$sample_type_binary == "Tumor"])
  mn <- mean(norm_counts[id, coldata$sample_type_binary == "Normal"])
  lfc <- res[id, "log2FoldChange"]; padj <- res[id, "padj"]
  cat(sprintf("%s (%s)\n  mean_tumor=%.2f  mean_normal=%.2f  higher_in_tumor=%s  log2FC=%.3f  padj=%.3g\n\n",
              gene, id, mt, mn, mt > mn, lfc, padj))
}


## ----------------------------------------------------------------------------
## 6. HEATMAP: 9-gene candidate panel, grouped by pathway and sample type
## ----------------------------------------------------------------------------

# Full versioned Ensembl IDs as they exist in dds (re-derived from the
# pattern-matched IDs above, so this stays correct even if GDC version
# suffixes change on a future data pull)
panel_genes <- candidate_ids
names(panel_genes) <- names(candidate_genes)

# Variance-stabilizing transform -- better suited to heatmap visualization
# than raw normalized counts, since it removes the mean-variance relationship
# inherent to count data.
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)[panel_genes, ]
rownames(vsd_mat) <- names(panel_genes)

# Row-wise (per-gene) z-scoring: t(scale(t(x))) is required here because
# base R's scale() normalizes by COLUMN by default. Since our matrix is
# genes (rows) x samples (columns), a plain scale(x) would incorrectly
# z-score each SAMPLE across genes rather than each GENE across samples --
# a common and easy-to-miss heatmap bug. Transposing twice ensures scaling
# happens along the correct margin.
vsd_scaled <- t(scale(t(vsd_mat)))

# CHECKPOINT: every gene (row) should now have mean ~0, sd ~1 across samples
round(rowMeans(vsd_scaled), 5)
round(apply(vsd_scaled, 1, sd), 5)

# Sample-type annotation, built directly from coldata and matched by real
# barcode (never by position) -- this is the second common heatmap bug this
# pipeline explicitly guards against.
annotation_col <- data.frame(
  sample_type = coldata$sample_type_binary,
  row.names = rownames(coldata)
)
identical(colnames(vsd_scaled), rownames(annotation_col))  # must be TRUE

# Gene-level pathway grouping (manually assigned, not clustering-derived) so
# the heatmap visually separates the apoptosis axis from the inflammation
# axis rather than relying on hierarchical clustering to (maybe) find it.
gene_groups <- data.frame(
  pathway = c(
    BAX = "Apoptosis", CASP3 = "Apoptosis", CASP9 = "Apoptosis", BCL2 = "Apoptosis",
    IL6 = "Inflammation", TNF = "Inflammation", IL10 = "Inflammation",
    TLR4 = "Inflammation", TLR2 = "Inflammation"
  ),
  row.names = c("BAX", "CASP3", "CASP9", "BCL2", "IL6", "TNF", "IL10", "TLR4", "TLR2")
)

# Reorder BOTH columns (by sample_type) and rows (by pathway group) on the
# same matrix object before plotting -- reordering only one dimension while
# leaving the other in original order is a subtle bug that will silently
# desynchronize the annotation bars from the data (encountered and fixed
# during development of this pipeline).
sample_order <- order(annotation_col$sample_type)
vsd_scaled_ordered <- vsd_scaled[, sample_order]
annotation_col_ordered <- annotation_col[sample_order, , drop = FALSE]

gene_order <- rownames(gene_groups)
vsd_scaled_grouped <- vsd_scaled_ordered[gene_order, ]

# CHECKPOINT: both alignments must hold before plotting
identical(colnames(vsd_scaled_grouped), rownames(annotation_col_ordered))
identical(rownames(vsd_scaled_grouped), rownames(gene_groups))

ann_colors <- list(
  sample_type = c(Normal = "cyan", Tumor = "pink"),
  pathway = c(Apoptosis = "darkgreen", Inflammation = "orange")
)

pheatmap(
  vsd_scaled_grouped,
  annotation_col = annotation_col_ordered,
  annotation_row = gene_groups,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,   # preserve manual Apoptosis/Inflammation block order
  cluster_cols = FALSE,   # preserve manual Normal/Tumor block order
  show_colnames = FALSE,  # 424 individual sample labels would be unreadable
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Apoptosis/Inflammation Panel: TCGA-LIHC Tumor vs Normal",
  filename = "heatmap_9gene_panel_grouped.png",
  width = 10, height = 5
)


## ----------------------------------------------------------------------------
## 7. VOLCANO PLOT: genome-wide DE results, candidate panel highlighted
## ----------------------------------------------------------------------------

res_full_df <- as.data.frame(res)
res_full_df$ensembl_id <- gsub("\\.[0-9]+$", "", rownames(res_full_df))
res_full_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_full_df$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

panel_symbols <- names(panel_genes)

# CHECKPOINT: confirm the 9 candidate genes' stats match what was validated
# manually in Section 4 above.
res_full_df[res_full_df$symbol %in% panel_symbols, c("symbol", "log2FoldChange", "padj")]

EnhancedVolcano(
  res_full_df,
  lab = res_full_df$symbol,
  x = 'log2FoldChange',
  y = 'padj',
  selectLab = panel_symbols,   # label only the 9 candidate genes, though all
                                # ~49,000 tested genes are plotted as points
  title = 'TCGA-LIHC: Tumor vs Normal',
  subtitle = 'Apoptosis/Inflammation panel highlighted',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 1.5,
  labSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  legendPosition = 'right'
)
ggsave("volcano_9gene_highlighted.png", width = 10, height = 8, dpi = 300)


## ----------------------------------------------------------------------------
## 8. GO / KEGG ENRICHMENT (genome-wide characterization, not panel-building)
## ----------------------------------------------------------------------------
## NOTE ON SCOPE: this enrichment is run on the full set of significant DEGs
## as a general characterization of TCGA-LIHC tumor biology -- it is
## deliberately NOT used to select or justify the 9-gene candidate panel
## above, which was chosen a priori from literature/mechanism. Running
## enrichment on ~9,000+ DEGs (a large fraction of the transcriptome, a
## consequence of this dataset's large sample size) tends to surface very
## broad, non-specific categories; see README for the specific finding that
## the KEGG apoptosis pathway is NOT significantly enriched genome-wide here,
## despite individual apoptosis genes (BAX, CASP3) being strongly significant.

deg_for_enrichment <- res_full_df[!is.na(res_full_df$padj) &
                                     res_full_df$padj < 0.05 &
                                     abs(res_full_df$log2FoldChange) > 1, ]
nrow(deg_for_enrichment)

deg_for_enrichment$entrez <- mapIds(
  org.Hs.eg.db,
  keys = deg_for_enrichment$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
entrez_ids <- unique(na.omit(deg_for_enrichment$entrez))
length(entrez_ids)

# Background universe = all genes actually tested (post-filtering), not the
# whole genome -- standard practice for correct enrichment background
# correction.
universe_entrez <- unique(na.omit(mapIds(
  org.Hs.eg.db,
  keys = res_full_df$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)))
length(universe_entrez)

go_enrich <- enrichGO(
  gene = entrez_ids,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

kegg_enrich <- enrichKEGG(
  gene = entrez_ids,
  universe = universe_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

dim(go_enrich@result)
dim(kegg_enrich@result)
head(go_enrich@result[, c("Description", "GeneRatio", "p.adjust")], 10)
head(kegg_enrich@result[, c("Description", "GeneRatio", "p.adjust")], 10)

# Subset specifically to apoptosis/inflammation-related terms, since these
# are diluted by very broad, non-specific top hits in this large gene list
# (see README).
go_results <- go_enrich@result
apop_inflam_go <- go_results[grepl("apoptosis|apoptotic|inflammat|cytokine|immune response|necrosis",
                                     go_results$Description, ignore.case = TRUE), ]
nrow(apop_inflam_go)
head(apop_inflam_go[, c("Description", "GeneRatio", "p.adjust")], 20)

kegg_results <- kegg_enrich@result
apop_inflam_kegg <- kegg_results[grepl("apoptosis|inflammat|immune|cytokine|NF-kappa|Toll-like",
                                         kegg_results$Description, ignore.case = TRUE), ]
nrow(apop_inflam_kegg)
head(apop_inflam_kegg[, c("Description", "GeneRatio", "p.adjust")], 20)

## RESULT: KEGG "Apoptosis" (hsa04210) and "Toll-like receptor signaling
## pathway" are NOT significantly enriched genome-wide (padj = 1.0) in this
## dataset, despite BAX/CASP3 and TLR2/TLR4/IL6/TNF/IL10 being individually
## significant. This is documented and discussed as a limitation in the
## README rather than omitted.


## ----------------------------------------------------------------------------
## 9. PCA: sample-level overview of tumor/normal separation
## ----------------------------------------------------------------------------

pca_data <- plotPCA(vsd, intgroup = "sample_type_binary", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
head(pca_data)
percentVar

ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_type_binary)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(Normal = "cyan3", Tumor = "deeppink3")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: TCGA-LIHC Tumor vs Normal") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())
ggsave("pca_tumor_normal.png", width = 8, height = 6, dpi = 300)

# Optional: same PCA with the detailed 3-level sample_type (Primary Tumor /
# Recurrent Tumor / Solid Tissue Normal) for visual inspection only -- this
# does NOT change or re-run any analysis, it is a display-only overlay.
pca_data$sample_type_detailed <- coldata$sample_type[match(pca_data$name, rownames(coldata))]
table(pca_data$sample_type_detailed)

ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_type_detailed)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(
    "Solid Tissue Normal" = "cyan3",
    "Primary Tumor" = "deeppink3",
    "Recurrent Tumor" = "darkorange"
  )) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: TCGA-LIHC by Detailed Sample Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.title = element_blank())
ggsave("pca_detailed_sampletype.png", width = 8, height = 6, dpi = 300)

## ============================================================================
## END OF PIPELINE
## ============================================================================
