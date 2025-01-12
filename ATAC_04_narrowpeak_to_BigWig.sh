#!/bin/bash

# Exit immediately if any command fails
set -e

# Define input files and parameters
CHROM_SIZES="hg19.chrom.sizes"
NARROWPEAK_FILES=(
    #"GSM5171834_Day0-CD8-N-1_peaks.narrowPeak" 
    #"GSM5171835_Day0-CD8-N-2_peaks.narrowPeak"
    #"GSM5171848_Day7-CD8-N-CD19-1_peaks.narrowPeak"
    #"GSM5171849_Day7-CD8-N-CD19-2_peaks.narrowPeak"
    "GSM4058193_Day10_CD4_CM_CD19_1_peaks.narrowPeak"
    "GSM4058194_Day10_CD4_CM_CD19_2_peaks.narrowPeak"
    #"GSM5171864_Day14-CD8-N-CD19-1_peaks.narrowPeak" 
    #"GSM5171865_Day14-CD8-N-CD19-2_peaks.narrowPeak"
)
OUTPUT_BW_FILES=(
    #"day0_CD8_N_signal.bw"
    #"day7_CD8_N_signal.bw"
    "day10_CD4_CM_signal.bw"
    #"day14_CD8_N_signal.bw"
)

# Create temporary directory and set cleanup
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "Working in temporary directory: $TMPDIR"

# Step 1: Create genomic bins
echo "Creating 100bp bins across the genome"
bedtools makewindows -g "$CHROM_SIZES" -w 100 > "$TMPDIR/genome_100bp_bins.bed"
sort -k1,1 -k2,2n "$TMPDIR/genome_100bp_bins.bed" -o "$TMPDIR/genome_100bp_bins.bed"

# Step 2: Process each input file
TOTAL_SIGNAL=0
SIGNAL_FILES=()

for FILE in "${NARROWPEAK_FILES[@]}"; do
    if [ ! -f "$FILE" ]; then
        echo "Error: Input file $FILE not found"
        exit 1
    fi
    
    OUTPUT_FILE="$TMPDIR/${FILE%.narrowPeak}_signal_coverage.bedGraph"
    
    # Aggregate signalValue within each bin using bedtools map
    bedtools map -a "$TMPDIR/genome_100bp_bins.bed" -b "$FILE" -c 7 -o sum > "$OUTPUT_FILE"
    SIGNAL_FILES+=("$OUTPUT_FILE")
    
    # Calculate total signal
    FILE_SIGNAL=$(awk '{sum += $4} END {printf "%.0f", sum}' "$OUTPUT_FILE")
    TOTAL_SIGNAL=$(awk -v curr="$TOTAL_SIGNAL" -v new="$FILE_SIGNAL" 'BEGIN {printf "%.0f", curr + new}')
    
    echo "  Total signal value for $FILE: $FILE_SIGNAL"
done

echo "Total signal across all files: $TOTAL_SIGNAL"

# Step 3: Normalize each file to CPM
NORMALIZED_FILES=()

for FILE in "${SIGNAL_FILES[@]}"; do
    NORMALIZED_FILE="${FILE%_signal_coverage.bedGraph}_cpm.bedGraph"
    
    # Normalize signal values to CPM
    awk -v total="$TOTAL_SIGNAL" \
        'BEGIN {OFS="\t"} 
         {if(total > 0) print $1, $2, $3, ($4 / total) * 1000000; else print $1, $2, $3, 0}' \
        "$FILE" > "$NORMALIZED_FILE"
    
    NORMALIZED_FILES+=("$NORMALIZED_FILE")
    
    # Calculate statistics
    awk '{
        sum+=$4;
        sumsq+=$4^2;
        if($4>max) max=$4;
        n++
    } END {
        if(n > 0) {
            mean=sum/n;
            std=sqrt(sumsq/n - (sum/n)^2);
            print "  Mean CPM:", mean;
            print "  Max CPM:", max;
            print "  CV%:", (std/mean)*100
        }
    }' "$NORMALIZED_FILE"
done

# Step 4: Process each timepoint
echo "Processing timepoints and calculating correlations"
for i in {0..2}; do
    timepoint=$((i * 2))
    day_num=$((i * 7))  # 0, 7, 14
    
    echo "Processing Day $day_num"
    
    # Get the correct files from the normalized files array
    file1="${NORMALIZED_FILES[$timepoint]}"
    file2="${NORMALIZED_FILES[$timepoint+1]}"
    
    if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
        echo "Error: Missing normalized files for Day $day_num"
        continue
    fi
    
    echo "Day $day_num replicate correlation:"
    paste "$file1" "$file2" | \
    awk '{
        sum1+=$4; sum2+=$8;
        sum12+=$4*$8;
        sum11+=$4*$4;
        sum22+=$8*$8;
        n++
    } END {
        if(n > 0) {
            denom=sqrt((n*sum11-sum1*sum1)*(n*sum22-sum2*sum2));
            if(denom > 0) {
                correlation=(n*sum12-sum1*sum2)/denom;
                print "  Correlation:", correlation
            } else {
                print "  Correlation: NA"
            }
        }
    }'
    
    # Average the replicates
    paste "$file1" "$file2" | \
        awk 'BEGIN {OFS="\t"} 
             {print $1, $2, $3, ($4 + $8) / 2}' \
        > "$TMPDIR/averaged_day${day_num}_cpm.bedGraph"
    
    # Sort the averaged file
    LC_COLLATE=C sort -k1,1 -k2,2n "$TMPDIR/averaged_day${day_num}_cpm.bedGraph" \
        > "$TMPDIR/averaged_day${day_num}_cpm.sorted.bedGraph"
done

# Step 5: Create final bigWig files
echo "Creating final bigWig files"
for i in {0..2}; do
    day_num=$((i * 7))
    input_file="$TMPDIR/averaged_day${day_num}_cpm.sorted.bedGraph"
    output_file="${OUTPUT_BW_FILES[$i]}"
    
    if [ -f "$input_file" ]; then
        echo "Converting $input_file to $output_file"
        bedGraphToBigWig "$input_file" "$CHROM_SIZES" "$output_file"
    else
        echo "Error: Missing input file $input_file"
    fi
done

echo "Processing complete! BigWig files created for Day 0, Day 7, and Day 14"