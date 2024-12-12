"""
This script computes the distribution of attention scores or QK^T / sqrt(d_k) for a given model.
The script takes two arguments:
1. The type of distribution to compute: "attention" or "qk"
2. The model name (optional): The name of the model to use (default: "bert-base-uncased")
"""

import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
from transformers import AutoModel
import sys

# Dimensions of the model
batch_size = 64  # Number of sequences
seq_len = 1024    # Sequence length
num_heads = 1  # Number of attention heads

def main():
    model_name = "bert-base-uncased"
    if len(sys.argv) == 3:
        model_name = str(sys.argv[2])
    elif len(sys.argv) != 2:
        print("Usage: python test.py [attention/qk] [model=bert-base-uncased]")
        sys.exit(1)

    model = AutoModel.from_pretrained(model_name)

    # Dimensions needed for splitting the weight matrix and scaling the dot product
    d_model = model.config.hidden_size # Dimensionality of the model
    d_k = d_model // num_heads  # Dimensionality of the key and query vectors

    if sys.argv[1] == "attention":
        compute_attention_scores(model, d_model, d_k)
    elif sys.argv[1] == "qk":
        compute_qk_scores(model, d_model, d_k)
    else:
        print("Invalid argument")
        sys.exit(1)

def compute_attention_scores(model, d_model, d_k):
    input_ids = torch.randint(0, model.config.vocab_size, (batch_size, seq_len))
    if model.config.model_type == "gpt2":
        input_embeddings = model.wte(input_ids)
    else:
        input_embeddings = model.embeddings.word_embeddings(input_ids)
    
    if model.config.model_type == "gpt2":
        # Access attention weights from 1st layer
        W_qkv = model.h[0].attn.c_attn.weight
        W_q = W_qkv[:, :d_model]
        W_k = W_qkv[:, d_model:2 * d_model]
    else:
        W_q = model.encoder.layer[0].attention.self.query.weight
        W_k = model.encoder.layer[0].attention.self.key.weight
    
    Q = torch.matmul(input_embeddings, W_q.T)
    K = torch.matmul(input_embeddings, W_k.T)
    
    # Compute scaled dot-product QK^T / sqrt(d_k)
    scaling_factor = torch.sqrt(torch.tensor(d_k, dtype=torch.float32))
    QK_T = torch.matmul(Q, K.transpose(-2, -1)) / scaling_factor
    
    attention_scores = torch.softmax(QK_T, dim=-1)  # Normalize over the last dimension
    data = attention_scores.flatten().detach().numpy()
    title = "Distribution of Attention Scores"
    xlabel = "Attention Score"
    mean_value = np.mean(data)
    variance_value = np.var(data)
    print(f"Mean of Attention Scores: {mean_value:.4f}")
    print(f"Variance of Attention Scores: {variance_value:.4f}")
    
    plot(data, title, xlabel, density=True)

def compute_qk_scores(model, d_model, d_k):
    input_ids = torch.randint(0, model.config.vocab_size, (batch_size, seq_len))
    if model.config.model_type == "gpt2":
        input_embeddings = model.wte(input_ids)
    else:
        input_embeddings = model.embeddings.word_embeddings(input_ids)

    if model.config.model_type == "gpt2":
        # Access attention weights from 1st layer
        W_qkv = model.h[0].attn.c_attn.weight
        # Slice weight matrix along second dimension
        W_q = W_qkv[:, :d_model]
        W_k = W_qkv[:, d_model:2 * d_model]
    else:
        W_q = model.encoder.layer[0].attention.self.query.weight
        W_k = model.encoder.layer[0].attention.self.key.weight
    
    Q = torch.matmul(input_embeddings, W_q.T)
    K = torch.matmul(input_embeddings, W_k.T)
    
    # Compute scaled dot-product QK^T / sqrt(d_k)
    scaling_factor = torch.sqrt(torch.tensor(d_k, dtype=torch.float32))
    QK_T = torch.matmul(Q, K.transpose(-2, -1)) / scaling_factor  # Shape: [batch_size, seq_len, seq_len]

    data = QK_T.flatten().detach().numpy()
    title = "Distribution of QK^T / sqrt(d_k)"
    xlabel = "QK^T / sqrt(d_k)"
    mean_value = np.mean(data)
    variance_value = np.var(data)
    print(f"Mean of QK^T / sqrt(d_k): {mean_value:.4f}")
    print(f"Variance of QK^T / sqrt(d_k): {variance_value:.4f}")

    # Proportion of large values
    proportion_large_values = np.mean(np.abs(data) > 3)
    print(f"Proportion of values exceeding |3|: {proportion_large_values:.4f}")

    """
    HYPOTHESIS: W_q and W_k are making the dot product too large
    """
    mean_W_q = torch.mean(W_q).item()
    var_W_q = torch.var(W_q).item()
    mean_W_k = torch.mean(W_k).item()
    var_W_k = torch.var(W_k).item()
    
    print(f"Mean of W_q: {mean_W_q:.4f}, Variance of W_q: {var_W_q:.4f}")
    print(f"Mean of W_k: {mean_W_k:.4f}, Variance of W_k: {var_W_k:.4f}")

    plot(data, title, xlabel, density=True)

def plot(data, title, xlabel, density):
    # Plot histogram
    plt.hist(data, bins=100, density=density, alpha=0.75)
    plt.title(title)
    plt.xlabel(xlabel)
    if density:
        plt.ylabel("Density")
    else:
        plt.ylabel("Frequency")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
