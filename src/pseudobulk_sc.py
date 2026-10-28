"""Pseudobulk HepG2 / Jurkat raw single-cell CRISPRi Perturb-seq (Nadig et al. 2025) with the
Replogle et al. 2022 normalisation: library-size normalisation, per-gem-group z-score against
non-targeting cells, mean per perturbation. Streams the CSR matrix in row blocks (8 GB RAM)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import sys, os
import numpy as np, pandas as pd, h5py, anndata as ad

RAW = str(_cfg.RAW)
name = sys.argv[1]  # hepg2 | jurkat
f = h5py.File(os.path.join(RAW, f"{name}_raw_singlecell.h5ad"), "r")


def read_col(g, k):
    x = g[k]
    if "__categories" in g and k in g["__categories"]:  # legacy anndata categorical
        cats = g["__categories"][k][:].astype(str); return cats[x[:]]
    if isinstance(x, h5py.Group):  # categorical
        cats = x["categories"][:].astype(str); codes = x["codes"][:]
        return cats[codes]
    v = x[:]
    return v.astype(str) if v.dtype.kind in "SO" else v


obs = f["obs"]; var = f["var"]
idx_key = obs.attrs.get("_index", "_index")
cells = pd.Index(obs[idx_key][:].astype(str))
cols = {k: read_col(obs, k) for k in obs.keys() if k != idx_key and k != "__categories"}
O = pd.DataFrame(cols, index=cells)
print(O.columns.tolist()); print(O.head(3).T)
vkey = var.attrs.get("_index", "_index")
gid = var[vkey][:].astype(str)
gname = read_col(var, "gene_name") if "gene_name" in var else gid
X = f["X"]
n_cells, n_genes = len(cells), len(gid)
print("cells", n_cells, "genes", n_genes, "X", X.shape)

gene_col = "gene" if "gene" in O else [c for c in O if "gene" in c.lower()][0]
tgt = O[gene_col].astype(str).values
gem = O["gem_group"].astype(str).values if "gem_group" in O else np.zeros(n_cells).astype(str)
is_nt = tgt == "non-targeting"
print("NT cells", is_nt.sum(), "gem groups", len(set(gem)))

# keep genes present in the Replogle bulk gene space
keepg = set(pd.read_pickle(_cfg.PROC / "ctx_K562.pkl")[0].columns) | \
        set(pd.read_pickle(_cfg.PROC / "ctx_RPE1.pkl")[0].columns)
gmask = np.array([g in keepg for g in gname])
gsel = np.where(gmask)[0]; G = len(gsel)
remap = -np.ones(n_genes, int); remap[gsel] = np.arange(G)

BLOCK = 4000


UMI_OBS = O["UMI_count"].values.astype(float)


def blocks():
    for s in range(0, n_cells, BLOCK):
        e = min(n_cells, s + BLOCK)
        D = X[s:e, :][:, gsel].astype(np.float64)
        yield s, e, D, UMI_OBS[s:e]


umi = UMI_OBS
med = np.median(umi[is_nt])
# pass 1: control mean/sd per gem group on normalised counts
gems = sorted(set(gem)); gi = {g: i for i, g in enumerate(gems)}
S1 = np.zeros((len(gems), G)); S2 = np.zeros((len(gems), G)); N = np.zeros(len(gems))
for s, e, D, tot in blocks():
    Dn = D / np.maximum(umi[s:e, None], 1) * med
    for j in np.where(is_nt[s:e])[0]:
        k = gi[gem[s + j]]; S1[k] += Dn[j]; S2[k] += Dn[j] ** 2; N[k] += 1
mu = S1 / np.maximum(N[:, None], 1)
sd = np.sqrt(np.maximum(S2 / np.maximum(N[:, None], 1) - mu ** 2, 0))
sd[sd == 0] = np.nan
# pass 2: z-score per cell, accumulate per perturbation (gene) and per NT guide pseudo-replicates
perts = sorted(set(tgt[~is_nt])); pi = {p: i for i, p in enumerate(perts)}
Z = np.zeros((len(perts), G)); C = np.zeros(len(perts))
nt_idx = np.where(is_nt)[0]
rng = np.random.default_rng(0)
nt_group = np.full(n_cells, -1); n_nt_groups = 100
nt_group[nt_idx] = rng.integers(0, n_nt_groups, len(nt_idx))  # random NT pseudo-perturbations
ZN = np.zeros((n_nt_groups, G)); CN = np.zeros(n_nt_groups)
for s, e, D, tot in blocks():
    Dn = D / np.maximum(umi[s:e, None], 1) * med
    k = np.array([gi[g] for g in gem[s:e]])
    Zb = (Dn - mu[k]) / sd[k]
    Zb = np.nan_to_num(Zb, nan=0.0)
    for j in range(e - s):
        t = tgt[s + j]
        if t == "non-targeting":
            q = nt_group[s + j]; ZN[q] += Zb[j]; CN[q] += 1
        else:
            Z[pi[t]] += Zb[j]; C[pi[t]] += 1
    print("block", s, flush=True)
keep = C >= 20
Zm = Z[keep] / C[keep, None]; ZNm = ZN / np.maximum(CN[:, None], 1)
names = [f"{i}_{p}_P1P2_NA" for i, p in enumerate(np.array(perts)[keep])] + \
        [f"{100000 + i}_non-targeting_non-targeting_non-targeting" for i in range(n_nt_groups)]
obs_out = pd.DataFrame({"num_cells_filtered": np.r_[C[keep], CN]}, index=names)
ctrl_mean = mu.T @ N / N.sum()  # mean normalised UMI in controls
var_out = pd.DataFrame({"gene_name": gname[gsel], "mean": ctrl_mean}, index=gid[gsel])
A = ad.AnnData(X=np.vstack([Zm, ZNm]).astype(np.float32), obs=obs_out, var=var_out)
A.write_h5ad(os.path.join(RAW, f"{name}_normalized_bulk.h5ad"))
print("wrote", A.shape, "perturbations", keep.sum(), "median cells", np.median(C[keep]))


