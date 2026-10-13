"""Sample metadata for GSE113957 (Fleischer et al. 2018) from the two GEO series matrices:
platform, biopsy site and repository cell ID, matched to FPKM column titles."""
import gzip, os
import pandas as pd

RAW = r"D:\Study\RejuvBench\data\raw"


def load_meta():
    rows = []
    for g, plat in [("GPL16791", "HiSeq 2500"), ("GPL18573", "NextSeq 500")]:
        rec = {}
        for l in gzip.open(os.path.join(RAW, f"GSE113957-{g}_series_matrix.txt.gz"), "rt", errors="replace"):
            f = [x.strip('"') for x in l.rstrip("\n").split("\t")]
            if f[0] == "!Sample_title": rec["title"] = f[1:]
            elif f[0] == "!Sample_source_name_ch1": rec["site"] = [x.replace("Skin; ", "") for x in f[1:]]
            elif f[0] == "!Sample_characteristics_ch1" and f[1].startswith("cell id"):
                rec["cell_id"] = [x.split(": ", 1)[1] for x in f[1:]]
        for t, s, c in zip(rec["title"], rec["site"], rec["cell_id"]):
            rows.append({"title": t, "platform": plat, "site": s, "cell_id": c, "repository": "".join(ch for ch in c if ch.isalpha())})
    return pd.DataFrame(rows).set_index("title")


if __name__ == "__main__":
    import sys; sys.path.insert(0, os.path.dirname(__file__))
    import aging
    e, m = aging.load_fleischer()
    M = load_meta().reindex(m.index)
    print("unmatched", M.platform.isna().sum())
    D = m.join(M)
    D["age_bin"] = pd.cut(D.age, [17, 39, 59, 79, 100], labels=["18-39", "40-59", "60-79", "80+"])
    print(pd.crosstab(D.age_bin, D.platform)); print(pd.crosstab(D.age_bin, D.site)); print(pd.crosstab(D.age_bin, D.repository))
