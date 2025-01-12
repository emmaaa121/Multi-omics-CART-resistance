import scanpy as sc
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import joblib
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import shap
import matplotlib.pyplot as plt
import h5py
import os

# GPU configuration
os.environ['CUDA_VISIBLE_DEVICES'] = '2' 
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Neural Network Model definition (must match the saved model architecture)
class NeuralNet(nn.Module):
    def __init__(self, input_size, hidden_layers, dropout_rates):
        super(NeuralNet, self).__init__()
        layers = []
        current_size = input_size
        
        for hidden_size, dropout_rate in zip(hidden_layers, dropout_rates):
            layers.extend([
                nn.Linear(current_size, hidden_size),
                nn.BatchNorm1d(hidden_size),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            current_size = hidden_size
        
        self.feature_layers = nn.Sequential(*layers)
        self.final_layer = nn.Linear(current_size, 2)

    def forward(self, x):
        features = self.feature_layers(x)
        logits = self.final_layer(features)
        if self.training or not torch.is_grad_enabled():
            return torch.sigmoid(logits[:, 1:2])
        return logits

# Load data
print("Loading data...")
adata = sc.read_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_combat_CR_NR.h5ad')
adata = adata[adata.obs['cell_type'] == 'CD8'].copy()
adata = adata[adata.obs['leiden_0.9'].isin(['0', '1', '4', '7', '9', '10', '11', '13', '14'])].copy()

gene_names = pd.read_csv('/home/emma/result/CART/CD3E_CD4_high_confidence_clusters_CART_response_diff_result_last.csv')['names']
genes_to_use = [gene for gene in gene_names if gene in adata.var_names and gene not in ['CD8A', 'CD8B']]  # Exclude CD8A and CD8B


X = adata[:, genes_to_use].X
X = X.toarray() if not isinstance(X, np.ndarray) else X
y = np.where(adata.obs['response'].values == 'CR', 0, 1)

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Load the saved model and scaler
print("Loading model and scaler...")
checkpoint = torch.load('/home/emma/data/CART/CD8_MLP_complete_model_logits.pth')
model_architecture = checkpoint['model_architecture']

# Create and load model
model = NeuralNet(
    input_size=X.shape[1],
    hidden_layers=model_architecture['hidden_layers'],
    dropout_rates=model_architecture['dropout_rates']
).to(device)
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()

# Load and apply scaler
scaler = joblib.load('/home/emma/data/CART/CD8_MLP_scaler_logits.pkl')
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Convert to tensors
X_train_tensor = torch.tensor(X_train_scaled, dtype=torch.float32).to(device)
X_test_tensor = torch.tensor(X_test_scaled, dtype=torch.float32).to(device)

# SHAP Analysis
print("Starting SHAP analysis...")
background_size = 1000
rng = np.random.RandomState(42)
background_indices = rng.choice(len(X_train_tensor), background_size, replace=False)
background_data = X_train_tensor[background_indices]

# Create explainer
explainer = shap.DeepExplainer(model, background_data)

try:
    print("Calculating SHAP values...")
    shap_values = explainer.shap_values(X_test_tensor, check_additivity=False)
    
    if isinstance(shap_values, list):
        shap_values_final = shap_values[1]
    else:
        shap_values_final = shap_values
        
    if len(shap_values_final.shape) == 1:
        shap_values_final = shap_values_final.reshape(len(X_test_tensor), -1)
        
except RuntimeError as e:
    print("Memory limit exceeded, switching to batch processing...")
    shap_values_list = []
    batch_size = 100
    
    for i in range(0, len(X_test_tensor), batch_size):
        batch = X_test_tensor[i:i + batch_size]
        batch_shap_values = explainer.shap_values(batch, check_additivity=False)
        if isinstance(batch_shap_values, list):
            batch_shap_values = batch_shap_values[1]
        shap_values_list.append(batch_shap_values)
    
    shap_values_final = np.concatenate(shap_values_list, axis=0)

# Save SHAP values
print("Saving results...")
shap_values_path = '/home/emma/result/CART/CD8_MLP_shap_values.h5'
with h5py.File(shap_values_path, 'w') as hf:
    hf.create_dataset('shap_values', data=shap_values_final)
    hf.create_dataset('feature_names', data=np.array(genes_to_use, dtype='S'))

# Generate plots
print("Generating plots...")
plt.figure(figsize=(12, 8))

# Ensure feature names match shap_values
genes_to_use = list(genes_to_use)
assert len(genes_to_use) == shap_values_final.shape[1], "Mismatch in feature names and SHAP values dimensions"

# Debugging output
print("Type of genes_to_use:", type(genes_to_use))
print("Number of features:", len(genes_to_use))
print("SHAP values shape:", shap_values_final.shape)

# Check and fix SHAP values shape
if len(shap_values_final.shape) == 3 and shap_values_final.shape[-1] == 1:
    shap_values_final = shap_values_final.squeeze(-1)  # Remove the last dimension

# Debugging output
print("Corrected SHAP values shape:", shap_values_final.shape)


shap.summary_plot(
    shap_values_final,
    X_test_scaled,
    feature_names=genes_to_use,
    plot_type="bar",
    show=False
)
plt.tight_layout()
plt.savefig('/home/emma/result/CART/CD8_MLP_shap_bar.png', dpi=300, bbox_inches='tight')
plt.close()

plt.figure(figsize=(12, 8))
shap.summary_plot(
    shap_values_final,
    X_test_scaled,
    feature_names=genes_to_use,
    plot_type="violin",
    show=False
)
plt.tight_layout()
plt.savefig('/home/emma/result/CART/CD8_MLP_shap_violin.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nAnalysis complete! Files saved:")
print(f"- SHAP values: {shap_values_path}")
print("- Bar plot: /home/emma/result/CART/CD8_MLP_shap_bar.png")
print("- Violin plot: /home/emma/result/CART/CD8_MLP_shap_violin.png")