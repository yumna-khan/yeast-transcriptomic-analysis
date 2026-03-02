# Yumna Khan
# Feb 18, 2026

# Overview ----
# The purpose of this script is to perform differential expression analysis on Saccharomyces cerevisiae RNA-seq data across three biofilm developmental stages: Early, Thin, and Mature. Transcriptional shifts are explored using PCA, volcano plots, and Over-Representation Analysis (ORA), followed by expression profiling of key functional genes to track replicate consistency over time.

# --------------- 1. Load libraries ---------------
library(tximport)
library(DESeq2)
library(ashr)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)
library(pheatmap)
library(tidyverse)


# --------------- 2. Load Salmon quantification ---------------
# Recursive search because each sample has its own folder
files <- list.files("../data/quants", pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
files
sampleNames <- basename(dirname(files))

# Define condition (stage)
condition <- factor(c(
  "Mature", "Mature", "Mature",
  "Thin", "Thin", "Thin",
  "Early", "Early", "Early"
))

coldata <- data.frame(row.names = sampleNames, condition = condition)


# --------------- 3. Load GTF and create tx2gene ---------------
gtf <- rtracklayer::import("../data/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz")
gtf_tx <- gtf[gtf$type == "transcript"]

# Create mapping
tx2gene <- data.frame(
  TXNAME = gtf_tx$transcript_id,
  GENEID = gtf_tx$gene_id
)

head(tx2gene)


# --------------- 4. Import Salmon counts ---------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Sanity check salmon names = tx names
dim(txi$counts)

# Get the lines of one Salmon file
salmon_names <- read.table(files[1], header = TRUE, sep = "\t")

common_ids <- intersect(salmon_names$Name, tx2gene$TXNAME)
print(paste("Matching IDs:", length(common_ids)))
print(paste("Salmon IDs not in tx2gene:", length(setdiff(salmon_names$Name, tx2gene$TXNAME))))


# --------------- 5. Create DESeq2 dataset ---------------
dds <- DESeqDataSetFromTximport(txi,
  colData = coldata,
  design = ~condition
)
# Run DESeq2
dds <- DESeq(dds)


# --------------- 6. Shrink LFC for interpretation ---------------
# Thin vs Early
res_Thin_vs_Early <- lfcShrink(
  dds,
  contrast = c("condition", "Thin", "Early"),
  type = "ashr"
)

# Mature vs Early
res_Mature_vs_Early <- lfcShrink(
  dds,
  contrast = c("condition", "Mature", "Early"),
  type = "ashr"
)

# Mature vs Thin
res_Mature_vs_Thin <- lfcShrink(
  dds,
  contrast = c("condition", "Mature", "Thin"),
  type = "ashr"
)

# Function for total DEGs, up and down regulated
count_DEGs <- function(res, padj_cutoff=0.05, log2FC_cutoff=1){
  DEGs <- res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > log2FC_cutoff), ]
  n_up <- nrow(DEGs[DEGs$log2FoldChange > log2FC_cutoff, ])
  n_down <- nrow(DEGs[DEGs$log2FoldChange < -log2FC_cutoff, ])
  n_total <- n_up + n_down
  return(list(total=n_total, up=n_up, down=n_down))
}

# Call function
count_DEGs(res_Thin_vs_Early)
count_DEGs(res_Mature_vs_Early)
count_DEGs(res_Mature_vs_Thin)


# --------------- 7. PCA (overall structure) ---------------
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA of Velum Stages")


# --------------- 8. Heatmap of top 50 DEGs ---------------
# Function for heatmap of 50 DEGs
top50_genes <- function(res_stage, vsd_stage, coldata_stage, title_text) {
  res_ordered <- res_stage[order(res_stage$padj), ]
  top50_genes <- rownames(head(res_ordered, 50))
  plot_data <- assay(vsd_stage)[top50_genes, ]
  
  # Z-score normalization
  plot_data <- t(scale(t(plot_data)))

  pheatmap(plot_data,
    annotation_col = coldata_stage,
    main = title_text,
    show_colnames = FALSE
  )
}

top50_genes(res_Thin_vs_Early, vsd, coldata, "Top 50 DEGs: Thin vs Early")
top50_genes(res_Mature_vs_Early, vsd, coldata, "Top 50 DEGs: Mature vs Early")
top50_genes(res_Mature_vs_Thin, vsd, coldata, "Top 50 DEGs: Mature vs Thin")


# --------------- 9. Volcano plot ---------------
# Create dataframe for all velum stages and add the gene column
res_df_tve <- res_Thin_vs_Early %>%
  as.data.frame() %>%
  rownames_to_column("gene")

res_df_mve <- res_Mature_vs_Early %>%
  as.data.frame() %>%
  rownames_to_column("gene")

res_df_mvt <- res_Mature_vs_Thin %>%
  as.data.frame() %>%
  rownames_to_column("gene")

# Function for volcano plot
res_df_volcano <- function(res_df_velum, title_text) {
  res_df_velum$significant <- ifelse(res_df_velum$padj < 0.05 & abs(res_df_velum$log2FoldChange) > 1, ifelse(res_df_velum$log2FoldChange > 0, "Up", "Down"), "Not Sig")

  res_df_velum <- na.omit(res_df_velum)

  ggplot(res_df_velum, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
    theme_minimal() +
    labs(
      title = title_text,
      x = "Log2 Fold Change",
      y = "-log10 Adjusted p-value"
    )
}

res_df_volcano(res_df_tve, "Volcano plot: Thin vs Early")
res_df_volcano(res_df_mve, "Volcano plot: Mature vs Early")
res_df_volcano(res_df_mvt, "Volcano plot: Mature vs Thin")


# --------------- 10. ORA Functional Analysis ---------------
# Function to run GO analysis
go_analysis <- function(res_df, plot_title) {
  # Define significant genes
  sig_genes <- res_df %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  # Define background genes (the universe)
  all_genes <- res_df %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  # Run GO Analysis
  ego <- enrichGO(
    gene          = sig_genes,
    universe      = all_genes,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )

  # Plot
  print(dotplot(ego) + ggtitle(plot_title))

  # Return the ego object to use it later
  return(ego)
}

ego_tve <- go_analysis(res_df_tve, "GO Enrichment (ORA): Thin vs Early")
ego_mve <- go_analysis(res_df_mve, "GO Enrichment (ORA): Mature vs Early")
ego_mvt <- go_analysis(res_df_mvt, "GO Enrichment (ORA): Mature vs Thin")


# --------------- 12. Gene Functional Annotation ---------------
# Function to retrieve the top gene ID for top 3 metabolic pathways
get_top_gene <- function(ego_obj, desc) {
  as.data.frame(ego_obj) %>%
    filter(Description == desc) %>%
    pull(geneID) %>%
    strsplit("/") %>%
    unlist() %>%
    head(1)
}

# Call function for each metabolic pathway
top_oxo_gene <- get_top_gene(ego_tve, "oxoacid metabolic process")
top_trans_gene <- get_top_gene(ego_mve, "transmembrane transport")
top_bio_gene <- get_top_gene(ego_mvt, "small molecule biosynthetic process")


# Identify IDs and Map Symbols
target_orfs <- c(top_oxo_gene, top_trans_gene, top_bio_gene)
gene_symbols <- mapIds(
  org.Sc.sgd.db,
  keys = target_orfs,
  column = "GENENAME",
  keytype = "ORF"
)

# Extract and reshape data
plot_data <- as.data.frame(assay(vsd)[target_orfs, ])
plot_data$Symbol <- gene_symbols

plot_data_long <- plot_data %>%
  pivot_longer(cols = -Symbol, names_to = "sample", values_to = "expression") %>%
  left_join(as.data.frame(coldata) %>%
    mutate(sample = rownames(coldata)), by = "sample")

# Assign replicate numbers (1, 2, 3)
plot_data_long <- plot_data_long %>%
  group_by(Symbol, condition) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

# Chronological order
plot_data_long$condition <- factor(plot_data_long$condition, levels = c("Early", "Thin", "Mature"))

# Loop to create and print 3 separate plots
for (gene in gene_symbols) {
  # Filter data for just one gene
  single_gene_data <- filter(plot_data_long, Symbol == gene)

  p <- ggplot(single_gene_data, aes(
    x = condition, y = expression,
    group = replicate, color = as.factor(replicate)
  )) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = paste("Expression Trajectory:", gene),
      x = "Biofilm Stage",
      y = "Normalized Expression (VST)",
      color = "Replicate"
    )

  print(p)
}
