"""Human orthologs of model-organism genes from the Alliance of Genome Resources combined orthology file
(pairs supported by >= 3 algorithms and best score in at least one direction)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import pandas as pd
RAW = str(_cfg.RAW)
df = pd.read_csv(RAW + r"\alliance_orthology.tsv.gz", sep="\t", comment="#")
tax = {"NCBITaxon:10090": "mmusculus", "NCBITaxon:6239": "celegans", "NCBITaxon:7227": "dmelanogaster", "NCBITaxon:559292": "scerevisiae"}
ok = (df.AlgorithmsMatch >= 3) & ((df.IsBestScore == "Yes") | (df.IsBestRevScore == "Yes"))
d = df[ok]
a = d[(d.Gene1SpeciesTaxonID == "NCBITaxon:9606") & d.Gene2SpeciesTaxonID.isin(tax)]
b = d[(d.Gene2SpeciesTaxonID == "NCBITaxon:9606") & d.Gene1SpeciesTaxonID.isin(tax)]
x = pd.concat([pd.DataFrame({"human": a.Gene1Symbol.values, "model": a.Gene2Symbol.values, "tax": a.Gene2SpeciesTaxonID.values}),
               pd.DataFrame({"human": b.Gene2Symbol.values, "model": b.Gene1Symbol.values, "tax": b.Gene1SpeciesTaxonID.values})])
for t, o in tax.items():
    s = x[x.tax == t][["human", "model"]].drop_duplicates()
    s["type"] = "best"
    s.to_csv(RAW + rf"\ortho_{o}.tsv", sep="\t", index=False)
    print(o, len(s))
