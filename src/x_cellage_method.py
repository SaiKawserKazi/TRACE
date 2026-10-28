"""Exploratory (Amendment 4): CellAge manipulation type (A9).
The CellAge Build 3 web table (data/raw/cellage_web_table.json; /cells/data/data_cellage.php) lists the method for each gene.
'Senescence effect' describes the gene's own function (e.g. ABCC6 knockdown causes senescence -> 'Inhibits'), so the
labels already state the expected effect of knockdown. Sensitivity: labels restricted to genes with loss-of-function
evidence (knockdown or knockout), the closest match to CRISPRi. Output: results/x_cellage_method.csv"""
import os, sys, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import trace_analysis as T
import rev3_scores as S

RAW, RES = T.RAW, T.RES; CTX = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]


def main():
    W = pd.DataFrame(json.load(open(os.path.join(RAW, "cellage_web_table.json"), encoding="utf-8"))["data"]).set_index("gene_symbol")
    lab = W.senescence_effect.where(W.senescence_effect.isin(["Induces", "Inhibits"]))
    lof = W.method.str.contains("Knockdown|Knockout", case=False)
    rows = []
    for sig, folder in [("Discovery signature", RES), ("GTEx-fibroblast signature", os.path.join(RES, "e18_gtex_discovery"))]:
        rob = pd.read_csv(os.path.join(folder, "robust_scores.csv"), index_col=0)
        sc = {c: pd.read_csv(os.path.join(folder, f"scores_{c}.csv"), index_col=0) for c in CTX}
        raw = pd.concat([sc[c].R_raw for c in CTX], axis=1).mean(axis=1); pro = pd.concat([sc[c].prolif for c in CTX], axis=1).mean(axis=1)
        for subset, keep in [("all CellAge genes", pd.Series(True, index=W.index)), ("loss-of-function evidence", lof), ("overexpression only", ~lof)]:
            l_ = lab[keep]
            for sname, s in [("R_rob (V1)", rob.R_rob), ("mean raw score", raw), ("proliferation shift alone", pro)]:
                x, y = S.evid(s, l_, "Induces", "Inhibits")
                if len(x) >= 5 and len(y) >= 5:
                    rows.append({"signature": sig, "subset": subset, "score": sname, "n_pos": len(x), "n_neg": len(y), "auroc": S.auc(x, y),
                                 "min_detectable_auroc_80pct": S.min_detectable_auc(len(x), len(y))})
    D = pd.DataFrame(rows); D.to_csv(os.path.join(RES, "x_cellage_method.csv"), index=False)
    print(W.method.value_counts().head(6).to_string()); print(D.round(3).to_string())


if __name__ == "__main__":
    main()
