import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scrublet as scr

adata_harvard = sc.read_h5ad('/home/emma/data/CART/Harvard_Yescarta_infusion_D7sorted_raw.h5ad')
adata_stanford = sc.read_h5ad('/home/emma/data/CART/Stanford_Yescarta_D7sorted_raw.h5ad')
adata_stanford.obs['timepoint'] = 'D7-CART'

adata_combined = adata_harvard.concatenate(adata_stanford, batch_key = 'organization')
batch_map = {'0': 'Harvard', '1': 'Stanford'}
adata_combined.obs['organization'] = adata_combined.obs['organization'].apply(lambda x: batch_map.get(x, x))

# Mitochondrial genes
adata_combined.var['mt'] = adata_combined.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Preprocess the combined data
sc.pp.filter_cells(adata_combined, min_genes=200)
sc.pp.filter_genes(adata_combined, min_cells=3)

# Filter based on quality control metrics
upper_lim = np.quantile(adata_combined.obs.n_genes_by_counts.values, .98)
lower_lim = np.quantile(adata_combined.obs.n_genes_by_counts.values, .02)
adata_combined = adata_combined[(adata_combined.obs.n_genes_by_counts < upper_lim) & (adata_combined.obs.n_genes_by_counts > lower_lim)]
adata_combined = adata_combined[adata_combined.obs.pct_counts_mt < 15]

# Doublet detection with Scrublet
scrub = scr.Scrublet(adata_combined.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata_combined.obs['doublet_scores'] = doublet_scores
adata_combined.obs['predicted_doublets'] = predicted_doublets
adata_combined = adata_combined[~adata_combined.obs['predicted_doublets']]

adata_combined.raw = adata_combined

if 'CAR' in adata_combined.obs.columns:
    adata_combined.obs['CAR'] = adata_combined.obs['CAR'].astype(str)

adata_combined.write_h5ad('/home/emma/data/CART/Harvard_Stanford_Yescarta_infusion_D7sorted_raw.h5ad')