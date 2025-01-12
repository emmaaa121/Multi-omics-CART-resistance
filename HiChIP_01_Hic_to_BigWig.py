import numpy as np
import pandas as pd
import subprocess
import os

# File paths remain the same
day0_files = [
    "/home/emma/data/CART/GSM5171868_CD8_D0_HiChIP_H3K27ac_r1.hic",
    "/home/emma/data/CART/GSM5171869_CD8_D0_HiChIP_H3K27ac_r2.hic"
]

day10_files = [
    "/home/emma/data/CART/GSM5171870_CD8_D10N19_HiChIP_H3K27ac_r1.hic",
    "/home/emma/data/CART/GSM5171871_CD8_D10N19_HiChIP_H3K27ac_r2.hic"
]

def check_juicer_tools():
    try:
        subprocess.run(['which', 'juicer_tools'], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        print("juicer_tools not found. Downloading...")
        download_cmd = "wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar"
        subprocess.run(download_cmd, shell=True)
        os.system('chmod +x juicer_tools_1.22.01.jar')
        return 'juicer_tools_1.22.01.jar'

def extract_hic_tracks(hic_file, resolution, chromosome, output_prefix):
    if os.path.exists('juicer_tools_1.22.01.jar'):
        cmd = f"java -jar juicer_tools_1.22.01.jar dump norm KR {hic_file} {chromosome} BP {resolution} {output_prefix}_norm.txt"
    else:
        cmd = f"juicer_tools dump norm KR {hic_file} {chromosome} BP {resolution} {output_prefix}_norm.txt"
    
    print(f"Running command: {cmd}")
    subprocess.run(cmd, shell=True)
    
    try:
        CHROM_SIZES = "/home/emma/data/CART/hg19.chrom.sizes"
        chrom_size = None
        with open(CHROM_SIZES) as f:
            for line in f:
                chrom, size = line.strip().split()
                if chrom == chromosome:
                    chrom_size = int(size)
                    break
        
        norm_factors = pd.read_csv(f"{output_prefix}_norm.txt", sep='\t', header=None)
        norm_factors.columns = ['norm_factor']
        
        # Create bins with staggered extensions to avoid overlaps
        bins = np.arange(len(norm_factors)) * resolution
        track_df = pd.DataFrame({
            'chromosome': chromosome,
            'start': bins,
            'end': bins + resolution,
            'coverage': norm_factors['norm_factor']
        })
        
        # Apply extensions without creating overlaps
        track_df['start'] = track_df['start'].apply(lambda x: max(0, x - 75))
        track_df['end'] = track_df.apply(lambda row: min(row['end'] + 75, 
                                                        track_df.loc[track_df.index > row.name, 'start'].iloc[0] if row.name < len(track_df)-1 else chrom_size), 
                                       axis=1)
        
        if chrom_size:
            track_df = track_df[track_df['end'] <= chrom_size]
        
        track_df = track_df.replace([np.inf, -np.inf], 0)
        track_df = track_df.fillna(0)
        
        track_df.to_csv(f"{output_prefix}_5K_KR_norm_75extend_coverage.bedgraph", sep='\t', header=False, index=False)
        
        return track_df
    except Exception as e:
        print(f"Warning: Could not process normalization factors for {hic_file}: {str(e)}")
        return None

# Rest of the functions remain the same
def process_replicates(hic_files, resolution, chromosome, output_prefix):
    coverage_tracks = []
    for i, hic_file in enumerate(hic_files):
        rep_prefix = f"{output_prefix}_rep{i+1}"
        track_df = extract_hic_tracks(hic_file, resolution, chromosome, rep_prefix)
        if track_df is not None:
            coverage_tracks.append(track_df)
    
    if coverage_tracks:
        merged_track = coverage_tracks[0].copy()
        merged_track['coverage'] = np.mean([track['coverage'] for track in coverage_tracks], axis=0)
        merged_track.to_csv(f"{output_prefix}_merged_5K_KR_norm_75extend_coverage.bedgraph", sep='\t', header=False, index=False)

def prepare_washu_track(bedgraph_file, track_name):
    sorted_bg = bedgraph_file.replace('.bedgraph', '.sorted.bedgraph')
    cmd = f"sort -k1,1 -k2,2n {bedgraph_file} > {sorted_bg}"
    subprocess.run(cmd, shell=True)
    
    output_bw = bedgraph_file.replace('.bedgraph', '.bw')
    CHROM_SIZES="/home/emma/data/CART/hg19.chrom.sizes"
    cmd = f"bedGraphToBigWig {sorted_bg} {CHROM_SIZES} {output_bw}"
    subprocess.run(cmd, shell=True)
    
    os.remove(sorted_bg)

def main():
    juicer_path = check_juicer_tools()
    resolution = 5000
    chromosome = "chr11"
    
    process_replicates(day0_files, resolution, chromosome, "CD8_D0_HiChIP_H3K27ac")
    process_replicates(day10_files, resolution, chromosome, "CD8_D10_HiChIP_H3K27ac")
    
    prepare_washu_track(
        "CD8_D0_HiChIP_H3K27ac_merged_5K_KR_norm_75extend_coverage.bedgraph",
        "CD8_D0_HiChIP_H3K27ac_5kb_KR_norm_75extend_track.bw"
    )
    prepare_washu_track(
        "CD8_D10_HiChIP_H3K27ac_merged_5K_KR_norm_75extend_coverage.bedgraph",
        "CD8_D10_HiChIP_H3K27ac_5kb_KR_norm_75_extend_track.bw"
    )

if __name__ == "__main__":
    main()