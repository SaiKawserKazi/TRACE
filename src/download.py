import requests, os, sys, time
D = r"D:\Study\RejuvBench\data\raw"
files = {
 "K562_essential_normalized_bulk_01.h5ad": "https://ndownloader.figshare.com/files/35773070",
 "rpe1_normalized_bulk_01.h5ad": "https://ndownloader.figshare.com/files/35775512",
 "K562_gwps_normalized_bulk_01.h5ad": "https://ndownloader.figshare.com/files/35773217",
 "GSE113957_fpkm.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113957/suppl/GSE113957_fpkm.txt.gz",
 "GSE113957_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113957/matrix/GSE113957_series_matrix.txt.gz",
 "gtex_fibroblasts_tpm.gct.gz": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cells_cultured_fibroblasts.gct.gz",
 "gtex_skin_ns_tpm.gct.gz": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_skin_not_sun_exposed_suprapubic.gct.gz",
 "gtex_subject.txt": "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt",
 "gtex_sample.txt": "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
 "string_links.txt.gz": "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
 "string_info.txt.gz": "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
 "genage_models.zip": "https://genomics.senescence.info/genes/models_genes.zip",
 "genage_human.zip": "https://genomics.senescence.info/genes/human_genes.zip",
 "cellage.zip": "https://genomics.senescence.info/cells/cellAge.zip",
 "hallmark.txt": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020",
 "oncokb_cancer_genes.tsv": "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt",
 "hepg2_raw_singlecell.h5ad": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad",
 "jurkat_raw_singlecell.h5ad": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad",
}
for name, url in files.items():
    p = os.path.join(D, name)
    if os.path.exists(p) and os.path.getsize(p) > 0: print("have", name, flush=True); continue
    t = time.time()
    try:
        with requests.get(url, stream=True, timeout=120, headers={"User-Agent":"Mozilla/5.0"}) as r:
            r.raise_for_status()
            with open(p + ".part", "wb") as f:
                for c in r.iter_content(1 << 22): f.write(c)
        os.replace(p + ".part", p)
        print("ok", name, os.path.getsize(p), round(time.time()-t), "s", flush=True)
    except Exception as e:
        print("FAIL", name, repr(e)[:200], flush=True)
