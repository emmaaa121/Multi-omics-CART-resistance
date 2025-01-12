import scanpy as sc
import pandas as pd
import numpy as np
import joblib
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (accuracy_score, classification_report, balanced_accuracy_score,
                           roc_auc_score, average_precision_score, precision_recall_curve)
import optuna
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0' 
os.environ['XLA_FLAGS'] = '--xla_gpu_cuda_data_dir=/usr/lib/cuda'

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Set random seeds for reproducibility
torch.manual_seed(42)
np.random.seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed(42)
    torch.backends.cudnn.deterministic = True

# Data loading and preprocessing
adata = sc.read_h5ad('/home/emma/data/CART/Harvard_Stanford_infusion_D7sorted_CD3E_CD4_CD8A_highly_variable_combat_CR_NR.h5ad')
adata = adata[adata.obs['cell_type'] == 'CD8'].copy()
adata = adata[adata.obs['leiden_0.9'].isin(['0', '1', '4', '7', '9', '10', '11', '13', '14'])].copy()

gene_names = pd.read_csv('/home/emma/result/CART/CD3E_CD4_high_confidence_clusters_CART_response_diff_result_last.csv')['names']
genes_to_use = [gene for gene in gene_names if gene in adata.var_names and gene not in ['CD8A', 'CD8B']]  # Exclude CD8A and CD8B

# Extract matrix and labels
X = adata[:, genes_to_use].X
X = X.toarray() if not isinstance(X, np.ndarray) else X
y = np.where(adata.obs['response'].values == 'CR', 0, 1)

# Print dataset dimensions
print("\nDataset Summary:")
print(f"Total samples: {X.shape[0]}")
print(f"Number of features: {X.shape[1]}")
print(f"Class distribution: {np.bincount(y)}")

# First split into train_val and test
X_train_val, X_test, y_train_val, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Create final scaler for the complete training set
final_scaler = StandardScaler()
X_train_val_scaled = final_scaler.fit_transform(X_train_val)
X_test_scaled = final_scaler.transform(X_test)  # Transform test set using full training set's parameters

# Calculate class weights using training data
class_counts = np.bincount(y_train_val)
class_weights = torch.FloatTensor(len(class_counts) / (class_counts * len(class_counts))).to(device)
print(f"Class weights: {class_weights.cpu().numpy()}")

# Convert to PyTorch tensors for final training
X_train_val_tensor = torch.tensor(X_train_val_scaled, dtype=torch.float32).to(device)
y_train_val_tensor = torch.tensor(y_train_val, dtype=torch.float32).view(-1, 1).to(device)
X_test_tensor = torch.tensor(X_test_scaled, dtype=torch.float32).to(device)
y_test_tensor = torch.tensor(y_test, dtype=torch.float32).view(-1, 1).to(device)

class EarlyStopping:
    def __init__(self, patience=5, min_delta=1e-3):
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.best_loss = None
        self.best_model = None
        self.early_stop = False

    def __call__(self, model, val_loss):
        if self.best_loss is None:
            self.best_loss = val_loss
            self.save_checkpoint(model)
        elif val_loss > self.best_loss - self.min_delta:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_loss = val_loss
            self.save_checkpoint(model)
            self.counter = 0

    def save_checkpoint(self, model):
        self.best_model = {k: v.clone() for k, v in model.state_dict().items()}

class NeuralNet(nn.Module):
    def __init__(self, input_size, hidden_layers, dropout_rates):
        super(NeuralNet, self).__init__()
        # Feature layers
        layers = []
        current_size = input_size
        
        for hidden_size, dropout_rate in zip(hidden_layers, dropout_rates):
            linear = nn.Linear(current_size, hidden_size)
            nn.init.xavier_uniform_(linear.weight, gain=nn.init.calculate_gain('relu'))
            nn.init.zeros_(linear.bias)
            
            layers.extend([
                linear,
                nn.BatchNorm1d(hidden_size),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            current_size = hidden_size
        
        self.feature_layers = nn.Sequential(*layers)
        
        # Final layer
        self.final_layer = nn.Linear(current_size, 2)  # Output 2 logits
        nn.init.xavier_uniform_(self.final_layer.weight, gain=1.0)
        nn.init.zeros_(self.final_layer.bias)

    def forward(self, x):
        features = self.feature_layers(x)
        logits = self.final_layer(features)
        if self.training or not torch.is_grad_enabled():
            return torch.sigmoid(logits[:, 1:2])  # Return positive class probability during training
        return logits  # Return logits for SHAP analysis

def evaluate_model(y_true, y_pred, y_pred_proba):
    return {
        'accuracy': accuracy_score(y_true, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_true, y_pred),
        'auc_roc': roc_auc_score(y_true, y_pred_proba),
        'avg_precision': average_precision_score(y_true, y_pred_proba)
    }

def objective(trial):
    # Hyperparameters
    lr = trial.suggest_float('lr', 1e-5, 1e-3, log=True)
    batch_size = trial.suggest_categorical('batch_size', [256, 512, 1024, 2048])
    num_epochs = trial.suggest_categorical('num_epochs', [20, 30, 40])
    
    # Architecture options
    architecture_options = [
        (2048, 1024),
        (2048, 1024, 512),
        (2048, 1024, 512, 256),
        (4096, 2048, 1024),
        (4096, 2048, 1024, 512),
        (4096, 2048, 1024, 512, 256)
    ]
    
    selected_architecture = trial.suggest_categorical('architecture', list(range(len(architecture_options))))
    hidden_layers = architecture_options[selected_architecture]
    
    # Dropout rates
    dropout_rates = []
    for i, layer_size in enumerate(hidden_layers):
        if i == 0:
            dropout_rates.append(trial.suggest_float(f'dropout_first', 0.3, 0.5))
        elif i == len(hidden_layers) - 1:
            dropout_rates.append(trial.suggest_float(f'dropout_last', 0.2, 0.4))
        else:
            dropout_rates.append(trial.suggest_float(f'dropout_{i}', 0.3, 0.6))
    
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_metrics = []
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X_train_val, y_train_val)):
        # Get the training and validation splits for this fold
        X_train_fold = X_train_val[train_idx]
        y_train_fold = y_train_val[train_idx]
        X_val_fold = X_train_val[val_idx]
        y_val_fold = y_train_val[val_idx]
        
        # Scale the data using only training fold
        scaler_fold = StandardScaler()
        X_train_fold_scaled = scaler_fold.fit_transform(X_train_fold)
        X_val_fold_scaled = scaler_fold.transform(X_val_fold)  # Transform validation using training fold's parameters
        
        # Convert to tensors
        X_train = torch.tensor(X_train_fold_scaled, dtype=torch.float32).to(device)
        y_train = torch.tensor(y_train_fold, dtype=torch.float32).view(-1, 1).to(device)
        X_val = torch.tensor(X_val_fold_scaled, dtype=torch.float32).to(device)
        y_val = torch.tensor(y_val_fold, dtype=torch.float32).view(-1, 1).to(device)
        
        model = NeuralNet(
            input_size=X_train_val.shape[1],
            hidden_layers=list(hidden_layers),
            dropout_rates=dropout_rates
        ).to(device)
        
        optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=0.02)
        criterion = nn.BCELoss(weight=None)
        
        early_stopping = EarlyStopping(patience=5)
        best_val_metrics = None
        
        for epoch in range(num_epochs):
            # Training
            model.train()
            indices = torch.randperm(len(X_train))
            total_loss = 0
            
            for i in range(0, len(X_train), batch_size):
                batch_indices = indices[i:i+batch_size]
                batch_X = X_train[batch_indices]
                batch_y = y_train[batch_indices]
                
                outputs = model(batch_X)
                weights = class_weights[batch_y.long().squeeze()]
                loss = criterion(outputs.squeeze(), batch_y.view(-1)) * weights.mean()
                
                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
                optimizer.step()
                
                total_loss += loss.item()
            
            # Validation
            model.eval()
            with torch.no_grad():
                val_pred = model(X_val)
                val_loss = criterion(val_pred.squeeze(), y_val.view(-1))
                
                val_pred_binary = (val_pred > 0.5).float()
                metrics = evaluate_model(
                    y_val.cpu().numpy(),
                    val_pred_binary.cpu().numpy(),
                    val_pred.cpu().numpy()
                )
                
                if best_val_metrics is None or metrics['balanced_accuracy'] > best_val_metrics['balanced_accuracy']:
                    best_val_metrics = metrics
                    early_stopping.save_checkpoint(model)
                
                early_stopping(model, val_loss)
                if early_stopping.early_stop:
                    break
        
        cv_metrics.append(best_val_metrics['balanced_accuracy'])
        model.load_state_dict(early_stopping.best_model)
    
    return np.mean(cv_metrics)

# Create and run the optimization study
print("\nStarting hyperparameter optimization...")
study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=100)

# Print best parameters
print("\nBest trial:")
trial = study.best_trial
print(f"  Value: {trial.value}")
print("  Params: ")
for key, value in trial.params.items():
    print(f"    {key}: {value}")

# Get best parameters and architecture
best_params = study.best_trial.params
architecture_options = [
    (2048, 1024),
    (2048, 1024, 512),
    (2048, 1024, 512, 256),
    (4096, 2048, 1024),
    (4096, 2048, 1024, 512),
    (4096, 2048, 1024, 512, 256)
]
hidden_layers = architecture_options[best_params['architecture']]

# Recreate dropout rates
dropout_rates = []
for i, layer_size in enumerate(hidden_layers):
    if i == 0:
        dropout_rates.append(best_params['dropout_first'])
    elif i == len(hidden_layers) - 1:
        dropout_rates.append(best_params['dropout_last'])
    else:
        dropout_rates.append(best_params[f'dropout_{i}'])

print("\nTraining final model with best parameters on full training set...")

# Create final model
final_model = NeuralNet(
    input_size=X_train_val.shape[1],
    hidden_layers=list(hidden_layers),
    dropout_rates=dropout_rates
).to(device)

optimizer = optim.AdamW(final_model.parameters(), lr=best_params['lr'], weight_decay=0.02)
criterion = nn.BCELoss(weight=None)

# Training loop for final model
train_losses = []
for epoch in range(best_params['num_epochs']):
    final_model.train()
    indices = torch.randperm(len(X_train_val_tensor))
    total_loss = 0
    
    for i in range(0, len(X_train_val_tensor), best_params['batch_size']):
        batch_indices = indices[i:i+best_params['batch_size']]
        batch_X = X_train_val_tensor[batch_indices]
        batch_y = y_train_val_tensor[batch_indices]
        
        outputs = final_model(batch_X)
        weights = class_weights[batch_y.long().squeeze()]
        loss = criterion(outputs.squeeze(), batch_y.view(-1)) * weights.mean()
        
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(final_model.parameters(), max_norm=0.5)
        optimizer.step()
        
        total_loss += loss.item()
    
    avg_loss = total_loss / (len(X_train_val_tensor) // best_params['batch_size'])
    train_losses.append(avg_loss)
    print(f"Epoch {epoch+1}/{best_params['num_epochs']}, Loss: {avg_loss:.4f}")

# Final evaluation on test set
print("\nEvaluating on test set...")
final_model.eval()
with torch.no_grad():
    test_pred = final_model(X_test_tensor)
    test_pred_binary = (test_pred > 0.5).float()
    
    final_metrics = evaluate_model(
        y_test_tensor.cpu().numpy(),
        test_pred_binary.cpu().numpy(),
        test_pred.cpu().numpy()
    )
    
    print("\nFinal Test Set Metrics:")
    for metric_name, value in final_metrics.items():
        print(f"  {metric_name}: {value:.4f}")
    
    print("\nDetailed Classification Report:")
    print(classification_report(y_test_tensor.cpu().numpy(), test_pred_binary.cpu().numpy(), digits=4))

# Save model and artifacts
save_path = '/home/emma/data/CART/CD8_MLP_complete_model_logits.pth'
torch.save({
    'model_state_dict': final_model.state_dict(),
    'optimizer_state_dict': optimizer.state_dict(),
    'scaler_state': final_scaler,  # Using final_scaler here
    'training_history': {
        'train_losses': train_losses,
    },
    'best_parameters': best_params,
    'final_metrics': final_metrics,
    'model_architecture': {
        'hidden_layers': hidden_layers,
        'dropout_rates': dropout_rates
    },
    'class_weights': class_weights
}, save_path)

# Save scaler separately
scaler_path = '/home/emma/data/CART/CD8_MLP_scaler_logits.pkl'
joblib.dump(final_scaler, scaler_path)  # Using final_scaler here

print(f"\nModel and training history saved to {save_path}")
print(f"Scaler saved to {scaler_path}")

# Print training summary
print("\nTraining Summary:")
print(f"Best architecture: {hidden_layers}")
print(f"Final training loss: {train_losses[-1]:.4f}")
print("\nBest Parameters:")
for param, value in best_params.items():
    print(f"  {param}: {value}")

# Save study results
study_path = '/home/emma/result/CART/CD8_MLP_study_results_logits.pkl'
joblib.dump(study, study_path)
print(f"\nStudy results saved to {study_path}")