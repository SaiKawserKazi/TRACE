"""E8/E9 with DepMap (exploratory, added in revision): overlap of CellAge and GenAge label classes with DepMap 24Q2
CRISPR-inferred common-essential genes (Broad DepMap, Figshare+ doi:10.25452/figshare.plus.25880521.v1, CC BY 4.0),
next to the Hart CEGv2 overlap. Adds columns to results/rev3_scores_label_essential_overlap.csv and writes a Fisher test per label set."""
import os, sys, json
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import rev3_scores as S3

RAW, RES = S3.RAW, S3.RES
dep = set(pd.read_csv(os.path.join(RAW, "depmap24q2_CRISPRInferredCommonEssentials.csv")).iloc[:, 0].str.split(" ").str[0])
OV = pd.read_csv(os.path.join(RES, "rev3_scores_label_essential_overlap.csv"))
labs = S3.labels(); n_dep, fis = [], []
for _, r in OV.iterrows():
    lab = labs[r["labels"]][0]; g = set(lab[lab == r["class"]].index); n_dep.append(len(g & dep))
OV["in_DepMap_common_essential"] = n_dep; OV["frac_DepMap_common_essential"] = OV.in_DepMap_common_essential / OV.n
OV.to_csv(os.path.join(RES, "rev3_scores_label_essential_overlap.csv"), index=False)
for nm, (lab, pos, neg) in labs.items():
    a = OV[(OV["labels"] == nm) & (OV["class"] == pos)].iloc[0]; b = OV[(OV["labels"] == nm) & (OV["class"] == neg)].iloc[0]
    for col in ["DepMap_common_essential", "CEGv2"]:
        k1, k2 = int(a[f"in_{col}"]), int(b[f"in_{col}"])
        orr, p = stats.fisher_exact([[k1, a.n - k1], [k2, b.n - k2]])
        fis.append({"labels": nm, "reference": col, "positive_class": pos, "frac_pos": k1 / a.n, "negative_class": neg, "frac_neg": k2 / b.n,
                    "odds_ratio": orr, "fisher_p_two_sided": p})
F = pd.DataFrame(fis); F.to_csv(os.path.join(RES, "rev3_depmap_label_fisher.csv"), index=False)
json.dump({"depmap_release": "DepMap 24Q2 Public", "doi": "10.25452/figshare.plus.25880521.v1", "file": "CRISPRInferredCommonEssentials.csv",
           "md5": "76dd2c044e9aa42bdfa344c8681904f7", "n_common_essential": len(dep)}, open(os.path.join(RES, "rev3_depmap.json"), "w"), indent=1)
print(OV[["labels", "class", "n", "frac_CEGv2", "frac_DepMap_common_essential"]].round(3).to_string()); print(F.round(4).to_string())
