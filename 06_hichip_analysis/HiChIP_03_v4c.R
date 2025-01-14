# Load required libraries
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

# Set file paths for your data
day0_files <- c(
  "/Users/weiwu/Desktop/2nd project/Chip_Seq/GSM5171868_CD8_D0_HiChIP_H3K27ac_r1.hic",
  "/Users/weiwu/Desktop/2nd project/Chip_Seq/GSM5171869_CD8_D0_HiChIP_H3K27ac_r2.hic"
)
day10_files <- c(
  "/Users/weiwu/Desktop/2nd project/Chip_Seq/GSM5171870_CD8_D10N19_HiChIP_H3K27ac_r1.hic",
  "/Users/weiwu/Desktop/2nd project/Chip_Seq/GSM5171871_CD8_D10N19_HiChIP_H3K27ac_r2.hic"
)

# Function to run Juicebox dump command from R
run_juicebox_dump <- function(hic_file, chr, resolution, output_dir) {
  # Create output filename based on input file
  base_name <- basename(hic_file)
  output_file <- file.path(output_dir, paste0(base_name, "_dump.txt"))
  
  # Construct Juicebox command
  command <- sprintf(
    "java -jar juicer_tools.jar dump observed NONE %s %s %s BP %d %s",
    hic_file, chr, chr, resolution, output_file
  )
  
  # Run command
  system(command)
  
  return(output_file)
}

# Function to read and process Juicebox dump data
read_juicebox_dump <- function(file_path) {
  data <- read.table(file_path, header = FALSE)
  colnames(data) <- c("bin1", "bin2", "count")
  return(data)
}

# Function to process v4C data for a specific anchor
process_v4c <- function(dump_files, anchor_bin, resolution = 5000) {
  # Read and combine replicates
  replicate_data <- map(dump_files, read_juicebox_dump) %>%
    bind_rows(.id = "replicate")
  
  # Filter for anchor bin and normalize
  v4c_profile <- replicate_data %>%
    filter(bin1 == anchor_bin | bin2 == anchor_bin) %>%
    mutate(
      position = ifelse(bin1 == anchor_bin, bin2, bin1),
      distance = abs(position - anchor_bin) * resolution
    ) %>%
    group_by(replicate, position, distance) %>%
    summarize(count = sum(count), .groups = "drop")
  
  return(v4c_profile)
}

# Function to plot v4C profile
plot_v4c <- function(v4c_data, title = "V4C Profile", window_size = 1000000) {
  # Calculate mean and standard deviation
  summary_stats <- v4c_data %>%
    group_by(position, distance) %>%
    summarize(
      mean_count = mean(count),
      sd_count = sd(count),
      .groups = "drop"
    )
  
  # Create plot
  p <- ggplot(summary_stats, aes(x = distance)) +
    geom_ribbon(aes(ymin = mean_count - sd_count, 
                    ymax = mean_count + sd_count),
                alpha = 0.2) +
    geom_line(aes(y = mean_count)) +
    scale_x_continuous(limits = c(0, window_size),
                       labels = function(x) paste0(x/1e6, "Mb")) +
    labs(title = title,
         x = "Distance from anchor",
         y = "Interaction frequency") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  return(p)
}

# Example usage with complete workflow
main <- function() {
  # Set parameters
  chr <- "chr1"              # Replace with your chromosome of interest
  resolution <- 5000         # 5kb resolution
  anchor_bin <- 1000        # Replace with your bin of interest
  output_dir <- "dumps"      # Directory for dump files
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE)
  
  # Process Day 0 files
  day0_dumps <- map(day0_files, function(file) {
    run_juicebox_dump(file, chr, resolution, output_dir)
  })
  
  # Process Day 10 files
  day10_dumps <- map(day10_files, function(file) {
    run_juicebox_dump(file, chr, resolution, output_dir)
  })
  
  # Generate v4C data
  day0_data <- process_v4c(day0_dumps, anchor_bin, resolution)
  day10_data <- process_v4c(day10_dumps, anchor_bin, resolution)
  
  # Create plots
  p1 <- plot_v4c(day0_data, title = "Day 0 V4C Profile")
  p2 <- plot_v4c(day10_data, title = "Day 10 V4C Profile")
  
  # Save plots
  ggsave("day0_v4c.pdf", p1, width = 8, height = 4)
  ggsave("day10_v4c.pdf", p2, width = 8, height = 4)
}

# Run the main function
main()