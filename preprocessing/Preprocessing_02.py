import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import scanpy as sc
#os.environ['CUDA_VISIBLE_DEVICES'] = '1' 
import genomap.genoVis as gpv 
import tensorflow as tf 
import os
from combat.pycombat import pycombat

# os.environ['CUDA_VISIBLE_DEVICES'] = '1' 
#os.environ['XLA_FLAGS'] = '--xla_gpu_cuda_data_dir=/usr/lib/cuda'
# print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
adata = sc.read_h5ad('/home/emma/data/CART/Harvard_Stanford_Yescarta_infusion_D7sorted_raw.h5ad')
cd3_mask = adata.obs_vector('CD3E') > 0
adata = adata[cd3_mask]
cd4_mask = adata.obs_vector('CD4') > 0
cd8a_mask = adata.obs_vector('CD8A') > 0
#both_mask = cd4_mask & cd8a_mask
either_mask = cd4_mask | cd8a_mask
adata = adata[either_mask, :]
adata.raw = adata
adata.write_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_raw.h5ad')

#(3+1)
adata_raw = adata.copy()
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4) 
sc.pp.log1p(adata) 
sc.pp.highly_variable_genes(adata, min_mean=0.005, max_mean=3, min_disp=0.1)
#print("Number of highly variable genes after adjustment:", adata.var['highly_variable'].sum())
adata.var.loc[['XIST','RPS4Y1','RPS4Y2'],'highly_variable'] = False
adata.var.loc[adata.var.index.str.match('TR.V|IG.V'),'highly_variable'] = False
adata_raw = adata_raw[:, adata.var['highly_variable']]
adata_raw.write_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_raw.h5ad')

#Generate adata_log.h5ad, for sc.rank_genes_clusters which requires log transformation but not scaled, since it needs highly expressed genes to caculate foldchange, if you scale the data, fold change will be NaN.
adata = adata[:, adata.var['highly_variable']]
adata.write_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_log.h5ad')
adata_combat = adata.copy()
log1p_df = pd.DataFrame(adata_combat.X.toarray().T, index=adata_combat.var_names, columns=adata_combat.obs_names)
batch_labels = adata.obs['channel'].tolist()  
corrected_df = pycombat(log1p_df, batch_labels)
adata_combat.X = corrected_df.T
adata_combat.write_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_combat.h5ad')

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10, zero_center=True)
sc.tl.pca(adata, svd_solver='arpack')
sc.external.pp.harmony_integrate(adata,'channel')
adata.obsm['X_pca_harmony_sub'] = adata.obsm['X_pca_harmony'][:,0:35]
sc.pp.neighbors(adata,use_rep='X_pca_harmony_sub')

# Run UMAP
sc.tl.umap(adata)
sc.tl.tsne(adata)
sc.tl.leiden(adata,resolution=0.9,key_added='leiden_0.9')
sc.tl.leiden(adata,resolution=0.8,key_added='leiden_0.8')
sc.tl.leiden(adata,resolution=0.7,key_added='leiden_0.7')
sc.tl.leiden(adata,resolution=0.6,key_added='leiden_0.6')
sc.tl.leiden(adata,resolution=0.5,key_added='leiden_0.5')
sc.tl.leiden(adata,resolution=0.4,key_added='leiden_0.4')
sc.tl.leiden(adata,resolution=0.3,key_added='leiden_0.3')
sc.tl.leiden(adata,resolution=1,key_added='leiden_1')
sc.tl.leiden(adata,resolution=1.1,key_added='leiden_1.1')
sc.tl.leiden(adata,resolution=1.2,key_added='leiden_1.2')
sc.tl.leiden(adata,resolution=1.3,key_added='leiden_1.3')
sc.tl.leiden(adata,resolution=1.4,key_added='leiden_1.4')
sc.tl.leiden(adata,resolution=1.5,key_added='leiden_1.5')
sc.tl.louvain(adata)
adata.write_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_harmony_umap.h5ad')

# Run Genovis
pca_results = adata.obsm['X_pca_harmony_sub']
resVis1 = gpv.genoVis(pca_results, n_clusters=15, colNum=32, rowNum=32)  

adata.obs['cluster_15'] = resVis1[1]
adata.obsm['X_genovis'] = resVis1[0]
adata.write_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_harmony_genovis.h5ad')
print("Files saved for each group.")