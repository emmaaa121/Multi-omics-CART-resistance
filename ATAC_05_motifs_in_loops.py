import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

def find_motifs_in_region(motif_file, chrom, start, end, extend=5000):
    """Extract motifs within a specified genomic region (with optional extension)."""
    motifs_df = pd.read_csv(
        motif_file, 
        sep='\t',
        header=None,
        names=['chr', 'start', 'end', 'motif_name', 'score', 'strand']
    )
    
    region_motifs = motifs_df[
        (motifs_df['chr'] == chrom) &
        (motifs_df['start'] >= start - extend) &
        (motifs_df['end'] <= end + extend)
    ]
    
    # Add core_motif_name column
    region_motifs['core_motif_name'] = region_motifs['motif_name'].str.split('/').str[0]
    
    return region_motifs

def process_motif_scores(region_motifs):
    """Process motif scores by taking the mean score for each motif type."""
    motif_stats = region_motifs.groupby('core_motif_name').agg({
        'score': ['mean', 'count']
    }).reset_index()
    
    motif_stats.columns = ['motif', 'score', 'count']
    return motif_stats

def assign_motif_family(motif_name):
    """Assign motif family based on name patterns."""
    motif_families = {
        'CTCF': ['CTCF', 'BORIS'],
        'ETS': ['ETS', 'ELK', 'ETV'],
        'RUNT': ['RUNX', 'RUNT'],
        'Other': []  # Default category
    }
    
    for family, patterns in motif_families.items():
        if any(pattern in motif_name.upper() for pattern in patterns):
            return family
    return 'Other'

def visualize_publication_style(motif_stats, output_file):
    """Create a publication-style visualization of motif scores."""
    if len(motif_stats) == 0:
        print("No motif statistics to plot.")
        return
    
    # Sort by score and add rank
    motif_stats = motif_stats.sort_values('score', ascending=False)
    motif_stats['rank'] = range(len(motif_stats))
    
    # Assign families and colors
    motif_stats['family'] = motif_stats['motif'].apply(assign_motif_family)
    family_colors = {
        'CTCF': 'purple',
        'ETS': 'red',
        'RUNT': 'green',
        'Other': 'gray'
    }
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot points colored by family
    texts = []
    for family in family_colors:
        family_data = motif_stats[motif_stats['family'] == family]
        if not family_data.empty:
            ax.scatter(family_data['rank'], family_data['score'],
                      color=family_colors[family], label=family, s=50)
            
            # Add labels for top motifs in each family
            for idx, row in family_data.head(min(3, len(family_data))).iterrows():
                texts.append(ax.text(row['rank'], row['score'], row['motif'],
                           fontsize=8, ha='center', va='bottom'))
    
    # Adjust text positions to avoid overlap
    adjust_text(texts, 
               arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
               expand_points=(1.5, 1.5),
               force_text=(0.5, 0.5))
    
    # Customize the plot
    ax.set_xlabel('Score Rank')
    ax.set_ylabel('Motif Score')
    ax.set_title('Motifs in ATAC-seq peaks')
    
    # Only add legend if we have multiple families
    if len(motif_stats['family'].unique()) > 1:
        ax.legend(title='TF family', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# Main execution
def main():
    motif_file = "/Users/weiwu/Desktop/2nd_project/ATAC_seq_data/all_motif_locations.bed"
    
    # Define regions
    regions = {
        'CD81_Anchor1': ('chr11', 2320000, 2330000),
        'Gene_Body_Region': ('chr11', 2397407, 2418649),
        'Downstream_Anchor1': ('chr11', 2430000, 2440000),
        'Downstream_Anchor2': ('chr11', 2440000, 2450000)
    }
    
    # Process each region
    for region_name, (chrom, start, end) in regions.items():
        print(f"Processing {region_name}...")
        
        # Get motifs in region
        region_motifs = find_motifs_in_region(motif_file, chrom, start, end)
        
        if len(region_motifs) == 0:
            print(f"No motifs found in {region_name}")
            continue
        
        # Process motif scores
        motif_stats = process_motif_scores(region_motifs)
        
        if len(motif_stats) == 0:
            print(f"No significant motifs found in {region_name}")
            continue
            
        # Create visualization
        output_file = f"/Users/weiwu/Desktop/2nd_project/ATAC_seq_data/motifs_in_loop/{region_name}_motif_scores.png"
        visualize_publication_style(motif_stats, output_file)
        print(f"Score plot saved for {region_name}: {output_file}")

if __name__ == "__main__":
    main()