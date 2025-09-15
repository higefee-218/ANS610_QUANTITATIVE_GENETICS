import numpy as np
import pandas as pd

# --- Example: Build the Genomic Relationship Matrix (VanRaden, 2008) ---

# Raw genotype matrix W (individuals x markers)
# Genotypes coded as 0, 1, 2
W = np.array([
    [1, 2, 0],
    [1, 1, 1],
    [2, 2, 0],
    [0, 1, 2]
], dtype=float)

ind_names = ["Ind1", "Ind2", "Ind3", "Ind4"]
marker_names = ["M1", "M2", "M3"]

print("Raw genotype matrix W:")
print(pd.DataFrame(W, index=ind_names, columns=marker_names))

# Step 1. Calculate allele frequencies
n_ind = W.shape[0]
p = W.sum(axis=0) / (2 * n_ind)
print("\nAllele frequencies p:")
print(dict(zip(marker_names, p)))

# Step 2. Centered genotype matrix M = W - 2p
P = np.tile(2 * p, (n_ind, 1))   # expected genotype matrix
M = W - P
print("\nCentered genotype matrix M:")
print(pd.DataFrame(M, index=ind_names, columns=marker_names))

# Step 3. Scaling factor k = 2 * sum(p * (1 - p))
k = 2 * np.sum(p * (1 - p))
print("\nScaling factor k:")
print(k)

# Step 4. Genomic Relationship Matrix G = (M M') / k
G = (M @ M.T) / k
print("\nGenomic Relationship Matrix G:")
print(pd.DataFrame(np.round(G, 3), index=ind_names, columns=ind_names))
