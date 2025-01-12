# Load required libraries
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(DESeq2)
library(pheatmap)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(parallel)
library(cluster)
library(RColorBrewer)
library(org.Hs.eg.db)

# Define your gene sets
resistant_markers <- c("JUND", "ZNF683", "FCER1G", "DENND5A", "CD27", "CD81", "TSC22D3", "PECAM1")
various_markers <- c(
  "IL2RB", "IL7R", "CXCR3", "CCR7",
  "CD27", "CD28", "TNFRSF4", "TNFRSF9", "ICOS",
  "ENTPD1", "CD101", "CTLA4", "LAG3", "CD244", "BTLA", "CD274", "PDCD1", "HAVCR2", "TIGIT",
  "IL2RA", "CD38", "MKI67", "PRDM1", "EOMES", "IKZF2", "TBX21", "B3GAT1"
)

# Combine both gene sets and remove duplicates
all_genes <- unique(c(resistant_markers, various_markers))

# Define sample table (needed for DESeq2)
peak_files <- c(
  "GSM5171828_Day0-CD4-CM-1_peaks.narrowPeak", "GSM5171829_Day0-CD4-CM-2_peaks.narrowPeak", 
  "GSM5171830_Day0-CD4-N-1_peaks.narrowPeak", "GSM5171831_Day0-CD4-N-2_peaks.narrowPeak", 
  "GSM5171836_Day7-CD4-CM-CD19-1_peaks.narrowPeak", "GSM5171837_Day7-CD4-CM-CD19-2_peaks.narrowPeak", 
  "GSM5171840_Day7-CD4-N-CD19-1_peaks.narrowPeak", "GSM5171841_Day7-CD4-N-CD19-2_peaks.narrowPeak", 
  "GSM5171852_Day14-CD4-CM-CD19-1_peaks.narrowPeak", "GSM5171853_Day14-CD4-CM-CD19-2_peaks.narrowPeak", 
  "GSM5171856_Day14-CD4-N-CD19-1_peaks.narrowPeak", "GSM5171857_Day14-CD4-N-CD19-2_peaks.narrowPeak"
)

sample_table <- data.frame(
  File = peak_files,
  Sample = c(
    "Day0-CD4-CM-1", "Day0-CD4-CM-2", "Day0-CD4-N-1", "Day0-CD4-N-2",
    "Day7-CD4-CM-CD19-1", "Day7-CD4-CM-CD19-2", "Day7-CD4-N-CD19-1", "Day7-CD4-N-CD19-2",
    "Day14-CD4-CM-CD19-1", "Day14-CD4-CM-CD19-2", "Day14-CD4-N-CD19-1", "Day14-CD4-N-CD19-2"
  ),
  Condition = c(
    "Day0-CD4-CM", "Day0-CD4-CM", "Day0-CD4-N", "Day0-CD4-N",
    "Day7-CD4-CM-CD19", "Day7-CD4-CM-CD19", "Day7-CD4-N-CD19", "Day7-CD4-N-CD19",
    "Day14-CD4-CM-CD19", "Day14-CD4-CM-CD19", "Day14-CD4-N-CD19", "Day14-CD4-N-CD19"
  )
)

# Load the saved data
filtered_peaks <- readRDS("filtered_peaks.rds")
peak_annotation <- readRDS("peak_annotation.rds")
peak_annotation_df <- readRDS("peak_annotation_df.rds")
count_matrix <- readRDS("count_matrix.rds")

# DESeq2 analysis
sample_table$Condition <- as.factor(sample_table$Condition)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_table, design = ~ Condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)
res <- results(dds)

# Clustering using silhouette method
rld <- rlog(dds, blind = TRUE)

# Select the top 5000 most variable genes
topVar <- head(order(rowVars(assay(rld)), decreasing = TRUE), 5000)
mat <- assay(rld)[topVar, ]
rownames(mat) <- peak_annotation_df$SYMBOL[topVar]  # Add gene symbols to rows

# Determine maximum number of clusters
max_clusters <- min(15, ncol(mat) - 1)

# Calculate silhouette width for each cluster count
sil_width <- sapply(2:max_clusters, function(k) {
  kmeans_result <- kmeans(t(mat), centers = k)
  mean(silhouette(kmeans_result$cluster, dist(t(mat)))[, 3])
})

# Determine the optimal number of clusters
optimal_k <- which.max(sil_width) + 1
set.seed(123)
clusters <- kmeans(t(mat), centers = optimal_k)

# Create annotation for the heatmap
sampleClusters <- data.frame(
  Sample_ID = colnames(mat),
  Cluster = factor(clusters$cluster),
  Condition = sample_table$Condition
)
rownames(sampleClusters) <- colnames(mat)

# Create gene type annotation
gene_type <- rep("Other", length(rownames(mat)))
gene_type[rownames(mat) %in% resistant_markers] <- "Resistant"
gene_type[rownames(mat) %in% various_markers] <- "Various"
gene_type[rownames(mat) %in% intersect(resistant_markers, various_markers)] <- "Both"

# Create annotation colors
# Ensure optimal_k is at least 3 for color palette
if (optimal_k < 3) {
  optimal_k <- 3
}

# Adjust your color palette accordingly
ann_colors <- list(
  Cluster = setNames(brewer.pal(optimal_k, "Set1"), 1:optimal_k),
  Condition = setNames(brewer.pal(length(unique(sample_table$Condition)), "Set2"), 
                       unique(sample_table$Condition))
)


# Generate heatmap
mat_scaled <- t(scale(t(mat)))
# Determine whether to show rownames based on gene types
show_rownames_flag <- any(gene_type != "Other")
# Subset the matrix to include only marker genes
mat_scaled_subset <- mat_scaled[rownames(mat_scaled) %in% all_genes, ]
# Generate heatmap
pheatmap(mat_scaled_subset,
         annotation_col = sampleClusters[, c("Cluster", "Condition"), drop = FALSE],
         annotation_colors = ann_colors,
         show_rownames = show_rownames_flag,  # Use flag here
         show_colnames = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         main = "Global Clustering with Marker Genes Highlighted",
         fontsize = 8,
         filename = "global_clustering_markers.pdf",
         width = 12,
         height = 8)


# Extract and print information about marker genes
marker_genes <- rownames(mat)[rownames(mat) %in% all_genes]
marker_clusters <- data.frame(
  Gene = marker_genes,
  Cluster = cutree(hclust(dist(mat[marker_genes,])), k = optimal_k),
  Type = ifelse(marker_genes %in% resistant_markers, 
                "Resistant", 
                ifelse(marker_genes %in% various_markers, "Various", "Both"))
)

# Print summary of marker genes clustering
print("Marker Genes Cluster Distribution:")
for(i in 1:optimal_k) {
  cat(sprintf("\nCluster %d:\n", i))
  cluster_genes <- marker_clusters[marker_clusters$Cluster == i, ]
  print(cluster_genes)
}