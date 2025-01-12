import numpy as np
import pandas as pd
import subprocess
import os

# File paths
day0_files = [
    "/home/emma/data/CART/GSM5171868_CD8_D0_HiChIP_H3K27ac_r1.hic",
    "/home/emma/data/CART/GSM5171869_CD8_D0_HiChIP_H3K27ac_r2.hic"
]

day10_files = [
    "/home/emma/data/CART/GSM5171870_CD8_D10N19_HiChIP_H3K27ac_r1.hic",
    "/home/emma/data/CART/GSM5171871_CD8_D10N19_HiChIP_H3K27ac_r2.hic"
]

def check_juicer_tools():
    """Check if juicer_tools is installed and download if needed"""
    try:
        subprocess.run(['which', 'juicer_tools'], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        print("juicer_tools not found. Downloading...")
        download_cmd = "wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar"
        subprocess.run(download_cmd, shell=True)
        # Make executable and create alias
        os.system('chmod +x juicer_tools_1.22.01.jar')
        return 'juicer_tools_1.22.01.jar'

def call_loops_hiccups(hic_file, output_prefix):
    """
    Call loops using HiCCUPS from juicer tools with parameters from the paper
    """
    # Create output directory if it doesn't exist
    os.makedirs(f"{output_prefix}_hiccups", exist_ok=True)
    
    # Run HiCCUPS with paper parameters and ignore-sparsity flag
    if os.path.exists('juicer_tools_1.22.01.jar'):
        cmd = f"java -jar juicer_tools_1.22.01.jar hiccups "
        cmd += f"-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 --ignore-sparsity "
        cmd += f"{hic_file} {output_prefix}_hiccups"
    else:
        cmd = f"juicer_tools hiccups "
        cmd += f"-m 500 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 --ignore-sparsity "
        cmd += f"{hic_file} {output_prefix}_hiccups"
    
    print(f"Running HiCCUPS command: {cmd}")
    subprocess.run(cmd, shell=True)
    
    # Check for resolution-specific output files
    possible_files = {
        5000: [f"{output_prefix}_hiccups/postprocessed_pixels_5000.bedpe",
               f"{output_prefix}_hiccups/5000_loops.bedpe"],
        10000: [f"{output_prefix}_hiccups/postprocessed_pixels_10000.bedpe",
                f"{output_prefix}_hiccups/10000_loops.bedpe"]
    }
    
    loops_by_resolution = {}
    
    for resolution, filenames in possible_files.items():
        for file in filenames:
            if os.path.exists(file):
                print(f"Found loops file for {resolution}bp resolution: {file}")
                try:
                    loops_df = pd.read_csv(file, sep='\t')
                    # Save to resolution-specific BEDPE file
                    output_file = f"{output_prefix}_loops_{resolution}bp.bedpe"
                    loops_df.to_csv(output_file, sep='\t', index=False)
                    loops_by_resolution[resolution] = loops_df
                    print(f"Saved {resolution}bp resolution loops to {output_file}")
                    break
                except pd.errors.EmptyDataError:
                    print(f"Warning: Loops file {file} is empty")
                except Exception as e:
                    print(f"Error reading loops file: {e}")
    
    return loops_by_resolution

def process_replicates(hic_files, output_prefix):
    """
    Process multiple replicates and merge loops by resolution
    """
    all_loops_by_resolution = {5000: [], 10000: []}
    
    # First, call loops for each replicate
    for i, hic_file in enumerate(hic_files):
        rep_prefix = f"{output_prefix}_rep{i+1}"
        loops_by_resolution = call_loops_hiccups(hic_file, rep_prefix)
        
        # Read the generated BEDPE files
        for resolution in [5000, 10000]:
            bedpe_file = f"{rep_prefix}_loops_{resolution}bp.bedpe"
            if os.path.exists(bedpe_file):
                print(f"Reading file: {bedpe_file}")
                loops_df = pd.read_csv(bedpe_file, sep='\t')
                
                # Convert problematic columns to numeric, replacing errors with NaN
                numeric_cols = ['score', 'observed', 'expectedBL', 'expectedDonut', 'expectedH', 'expectedV', 
                              'fdrBL', 'fdrDonut', 'fdrH', 'fdrV', 'numCollapsed', 'centroid1', 'centroid2', 'radius']
                
                for col in numeric_cols:
                    if col in loops_df.columns:
                        loops_df[col] = pd.to_numeric(loops_df[col], errors='coerce')
                
                print(f"Columns in file: {loops_df.columns.tolist()}")
                all_loops_by_resolution[resolution].append(loops_df)
    
    # Merge replicates for each resolution
    for resolution, loops_list in all_loops_by_resolution.items():
        if loops_list:
            # Merge replicates
            merged_loops = pd.concat(loops_list)
            
            # Define columns
            coordinate_cols = ['#chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']
            string_cols = ['name', 'strand1', 'strand2', 'color']
            
            # Create aggregation dictionary
            agg_dict = {}
            for col in merged_loops.columns:
                if col in coordinate_cols:
                    continue
                elif col in string_cols:
                    agg_dict[col] = 'first'
                else:
                    agg_dict[col] = 'mean'
            
            print(f"\nGrouping by coordinates: {coordinate_cols}")
            print(f"Aggregation dictionary: {agg_dict}")
            
            # Group and aggregate
            merged_loops = merged_loops.groupby(coordinate_cols).agg(agg_dict).reset_index()
            
            # Save merged loops
            output_file = f"{output_prefix}_merged_loops_{resolution}bp.bedpe"
            merged_loops.to_csv(output_file, sep='\t', index=False)
            print(f"\nSaved merged loops to {output_file}")
            
            # Debug: print sample of merged data
            print("\nSample of merged data:")
            print(merged_loops.head())
            print("\nMerged data types:")
            print(merged_loops.dtypes)

def main():
    juicer_path = check_juicer_tools()
    
    # Process Day 0 samples
    process_replicates(
        day0_files,
        output_prefix="CD8_D0_HiChIP_H3K27ac"
    )
    
    # Process Day 10 samples
    process_replicates(
        day10_files,
        output_prefix="CD8_D10_HiChIP_H3K27ac"
    )

if __name__ == "__main__":
    main()
