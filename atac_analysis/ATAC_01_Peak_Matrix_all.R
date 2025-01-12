# Load required libraries
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(DESeq2)
library(pheatmap)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(parallel)
library(org.Hs.eg.db)

# Step 1: Read in the narrowPeak files using rtracklayer::import()
peak_files <- c(
  "GSM5171828_Day0-CD4-CM-1_peaks.narrowPeak",
  "GSM5171829_Day0-CD4-CM-2_peaks.narrowPeak", "GSM5171830_Day0-CD4-N-1_peaks.narrowPeak",
  "GSM5171831_Day0-CD4-N-2_peaks.narrowPeak", "GSM5171832_Day0-CD8-CM-1_peaks.narrowPeak",
  "GSM5171833_Day0-CD8-CM-2_peaks.narrowPeak", "GSM5171834_Day0-CD8-N-1_peaks.narrowPeak",
  "GSM5171835_Day0-CD8-N-2_peaks.narrowPeak", "GSM5171836_Day7-CD4-CM-CD19-1_peaks.narrowPeak",
  "GSM5171837_Day7-CD4-CM-CD19-2_peaks.narrowPeak", "GSM5171840_Day7-CD4-N-CD19-1_peaks.narrowPeak",
  "GSM5171841_Day7-CD4-N-CD19-2_peaks.narrowPeak", "GSM5171844_Day7-CD8-CM-CD19-1_peaks.narrowPeak",
  "GSM5171845_Day7-CD8-CM-CD19-2_peaks.narrowPeak", "GSM5171848_Day7-CD8-N-CD19-1_peaks.narrowPeak",
  "GSM5171849_Day7-CD8-N-CD19-2_peaks.narrowPeak", "GSM5171852_Day14-CD4-CM-CD19-1_peaks.narrowPeak",
  "GSM5171853_Day14-CD4-CM-CD19-2_peaks.narrowPeak", "GSM5171856_Day14-CD4-N-CD19-1_peaks.narrowPeak",
  "GSM5171857_Day14-CD4-N-CD19-2_peaks.narrowPeak", "GSM5171860_Day14-CD8-CM-CD19-1_peaks.narrowPeak",
  "GSM5171861_Day14-CD8-CM-CD19-2_peaks.narrowPeak", "GSM5171864_Day14-CD8-N-CD19-1_peaks.narrowPeak",
  "GSM5171865_Day14-CD8-N-CD19-2_peaks.narrowPeak"
)

# Create a sample table
sample_table <- data.frame(
  File = peak_files,
  Sample = c(
    "Day0-CD4-CM-1", "Day0-CD4-CM-2", "Day0-CD4-N-1", "Day0-CD4-N-2",
    "Day0-CD8-CM-1", "Day0-CD8-CM-2", "Day0-CD8-N-1", "Day0-CD8-N-2",
    "Day7-CD4-CM-CD19-1", "Day7-CD4-CM-CD19-2", "Day7-CD4-N-CD19-1", "Day7-CD4-N-CD19-2",
    "Day7-CD8-CM-CD19-1", "Day7-CD8-CM-CD19-2", "Day7-CD8-N-CD19-1", "Day7-CD8-N-CD19-2",
    "Day14-CD4-CM-CD19-1", "Day14-CD4-CM-CD19-2", "Day14-CD4-N-CD19-1", "Day14-CD4-N-CD19-2",
    "Day14-CD8-CM-CD19-1", "Day14-CD8-CM-CD19-2", "Day14-CD8-N-CD19-1", "Day14-CD8-N-CD19-2"
  ),
  Condition = c(
    "Day0-CD4-CM", "Day0-CD4-CM", "Day0-CD4-N", "Day0-CD4-N",
    "Day0-CD8-CM", "Day0-CD8-CM", "Day0-CD8-N", "Day0-CD8-N",
    "Day7-CD4-CM-CD19", "Day7-CD4-CM-CD19", "Day7-CD4-N-CD19", "Day7-CD4-N-CD19",
    "Day7-CD8-CM-CD19", "Day7-CD8-CM-CD19", "Day7-CD8-N-CD19", "Day7-CD8-N-CD19",
    "Day14-CD4-CM-CD19", "Day14-CD4-CM-CD19", "Day14-CD4-N-CD19", "Day14-CD4-N-CD19",
    "Day14-CD8-CM-CD19", "Day14-CD8-CM-CD19", "Day14-CD8-N-CD19", "Day14-CD8-N-CD19"
  )
)

# Step 2: Read in narrowPeak files in parallel and create GRanges objects
peak_ranges_list <- mclapply(peak_files, function(file) {
  peaks <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(peaks) >= 5) {
    GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3), strand = "*", score = peaks$V5)
  } else {
    warning(paste("Skipping file:", file, "- insufficient columns"))
    NULL
  }
}, mc.cores = detectCores() - 1)

peak_ranges_list <- Filter(Negate(is.null), peak_ranges_list)

# Step 3: Extend peaks to 500bp around summits
extended_peaks_list <- lapply(peak_ranges_list, function(peaks) {
  resize(peaks, width = 500, fix = 'center')
})

# Step 4: Merge peaks into union peak set and select highest scores
merged_peaks <- reduce(unlist(GRangesList(extended_peaks_list)), with.revmap = TRUE)
max_scores <- sapply(mcols(merged_peaks)$revmap, function(cluster) {
  max(sapply(cluster, function(i) {
    unlist(lapply(extended_peaks_list, function(gr) mcols(gr)$score[i]))
  }))
})

mcols(merged_peaks)$score <- max_scores

# Step 5: Filter peaks by ENCODE hg19 blacklist regions
blacklist <- import('hg19-blacklist.v2.bed', format = 'bed')
filtered_peaks <- subsetByOverlaps(merged_peaks, blacklist, invert = TRUE)
rm(merged_peaks, blacklist)
gc()

# Step 6: Peak annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peak_annotation <- annotatePeak(filtered_peaks, TxDb = txdb, tssRegion = c(-3000, 3000), level = 'gene', annoDb = 'org.Hs.eg.db')
peak_annotation_df <- as.data.frame(peak_annotation)

saveRDS(filtered_peaks, 'filtered_peaks_all.rds')
saveRDS(peak_annotation, 'peak_annotation_all.rds')
saveRDS(peak_annotation_df, 'peak_annotation_df_all.rds')

# Step 7: Generate count matrix
count_matrix <- matrix(0, nrow = length(filtered_peaks), ncol = length(peak_files))
rownames(count_matrix) <- paste0('peak_', 1:length(filtered_peaks))
colnames(count_matrix) <- sample_table$Sample

# Fill count matrix in parallel
count_matrix <- mclapply(1:length(peak_files), function(i) {
  peaks <- import(peak_files[i], format = 'narrowPeak')
  countOverlaps(filtered_peaks, peaks)
}, mc.cores = detectCores() - 1)

count_matrix <- do.call(cbind, count_matrix)
rownames(count_matrix) <- row_names
colnames(count_matrix) <- col_names

saveRDS(count_matrix, 'count_matrix_all.rds')
