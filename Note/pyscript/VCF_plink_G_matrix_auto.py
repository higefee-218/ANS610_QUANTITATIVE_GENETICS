import numpy as np
import pandas as pd
import allel
from plinkio import plinkfile

def calc_G(W, ind_names=None):
    n_ind = W.shape[0]
    p = W.sum(axis=0) / (2 * n_ind)
    P = np.tile(2 * p, (n_ind, 1))
    M = W - P
    k = 2 * np.sum(p * (1 - p))
    G = (M @ M.T) / k
    if ind_names is None:
        ind_names = [f"Ind{i+1}" for i in range(n_ind)]
    return pd.DataFrame(G, index=ind_names, columns=ind_names)

def calc_G_auto(input_file):
    if input_file.endswith(".vcf"):
        # --- VCF ---
        callset = allel.read_vcf(input_file)
        gt = callset['calldata/GT']
        W = gt.sum(axis=2).T.astype(float)
        W[W < 0] = np.nan
        col_means = np.nanmean(W, axis=0)
        inds = np.where(np.isnan(W))
        W[inds] = np.take(col_means, inds[1])
        ind_names = callset['samples']
        return calc_G(W, ind_names=ind_names)

    elif input_file.endswith(".bed"):
        # --- PLINK ---
        plink = plinkfile.open(input_file)
        samples = [s.iid for s in plink.get_samples()]
        loci = plink.get_loci()
        W = []
        for row in plink:
            geno = np.array(row, dtype=float)
            geno[geno < 0] = np.nan
            W.append(geno)
        plink.close()
        W = np.array(W).T
        col_means = np.nanmean(W, axis=0)
        inds = np.where(np.isnan(W))
        W[inds] = np.take(col_means, inds[1])
        return calc_G(W, ind_names=samples)

    else:
        raise ValueError("Input must be a VCF or PLINK .bed file.")

# Example usage:
# G1 = calc_G_auto("genotypes.vcf")
# G2 = calc_G_auto("genotypes.bed")
