"""Load CRISPRi Perturb-seq pseudobulk z-shift matrices per context."""
import os
import numpy as np, pandas as pd, anndata as ad

RAW = r"D:\Study\RejuvBench\data\raw"
PROC = r"D:\Study\RejuvBench\data\proc"

CONTEXTS = {
    "K562": "K562_essential_normalized_bulk_01.h5ad",
    "K562gw": "K562_gwps_normalized_bulk_01.h5ad",
    "RPE1": "rpe1_normalized_bulk_01.h5ad",
    "HepG2": "hepg2_normalized_bulk.h5ad",   # built by pseudobulk_sc.py
    "Jurkat": "jurkat_normalized_bulk.h5ad",
}


def load_context(name):
    """Returns (shift: perturbed-gene x measured-gene DataFrame, controls: NT x gene, gene_mean Series)."""
    cache = os.path.join(PROC, f"ctx_{name}.pkl")
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    a = ad.read_h5ad(os.path.join(RAW, CONTEXTS[name]))
    X = pd.DataFrame(np.asarray(a.X, dtype=np.float32), index=a.obs.index, columns=a.var.gene_name.values)
    X = X.T.groupby(level=0).mean().T  # duplicate gene symbols
    X = X.loc[:, np.isfinite(X.values).all(0)]  # drop genes with non-finite z-scores
    tgt = pd.Series([i.split("_")[1] for i in X.index], index=X.index)
    nt = tgt == "non-targeting"
    controls = X[nt.values]
    shift = X[~nt.values].groupby(tgt[~nt].values).mean()  # average guides per gene
    ncells = a.obs.loc[~nt.values, "num_cells_filtered"].groupby(tgt[~nt].values).sum() \
        if "num_cells_filtered" in a.obs else None
    gmean = pd.Series(a.var["mean"].values, index=a.var.gene_name.values).groupby(level=0).mean()
    out = (shift, controls, gmean, ncells)
    pd.to_pickle(out, cache)
    return out


if __name__ == "__main__":
    import sys
    for c in sys.argv[1:]:
        s, nt, gm, nc = load_context(c)
        print(c, "perturbed genes", s.shape, "NT", nt.shape)

