"""Data manifest: every external input with source URL/accession, access date, size, SHA-256 and citation.
Third-party data are not redistributed; src/download.ps1 and src/pardl.py retrieve them."""
import hashlib, os, csv

RAW = r"D:\Study\RejuvBench\data\raw"
ACCESS = "2026-10-12"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
ITEMS = [
    ("GSE113957_fpkm.txt.gz", GEO + "GSE113nnn/GSE113957/suppl/GSE113957_fpkm.txt.gz", "GEO GSE113957", "Fleischer et al. 2018 Genome Biol", "dermal fibroblast FPKM"),
    ("GSE113957-GPL16791_series_matrix.txt.gz", GEO + "GSE113nnn/GSE113957/matrix/GSE113957-GPL16791_series_matrix.txt.gz", "GEO GSE113957", "Fleischer et al. 2018", "sample metadata (HiSeq 2500); accessed 2026-10-13"),
    ("GSE113957-GPL18573_series_matrix.txt.gz", GEO + "GSE113nnn/GSE113957/matrix/GSE113957-GPL18573_series_matrix.txt.gz", "GEO GSE113957", "Fleischer et al. 2018", "sample metadata (NextSeq 500); accessed 2026-10-13"),
    ("gtex_fibroblasts_tpm.gct.gz", "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cells_cultured_fibroblasts.gct.gz", "GTEx v10", "GTEx Consortium 2020 Science", "open-access tier"),
    ("gtex_skin_ns_tpm.gct.gz", "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_skin_not_sun_exposed_suprapubic.gct.gz", "GTEx v10", "GTEx Consortium 2020", "open-access tier"),
    ("gtex_subject.txt", "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", "GTEx v10", "GTEx Consortium 2020", "open-access subject annotations (age in 10-year bins)"),
    ("gtex_sample.txt", "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", "GTEx v10", "GTEx Consortium 2020", "open-access sample attributes"),
    ("K562_essential_normalized_bulk_01.h5ad", "https://ndownloader.figshare.com/files/35780870", "Figshare article 20029387", "Replogle et al. 2022 Cell", "pseudobulk z-scores"),
    ("K562_gwps_normalized_bulk_01.h5ad", "https://ndownloader.figshare.com/files/35773217", "Figshare article 20029387", "Replogle et al. 2022", "pseudobulk z-scores"),
    ("rpe1_normalized_bulk_01.h5ad", "https://ndownloader.figshare.com/files/35775512", "Figshare article 20029387", "Replogle et al. 2022", "pseudobulk z-scores"),
    ("hepg2_raw_singlecell.h5ad", GEO + "GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad", "GEO GSE264667", "Nadig et al. 2025 Nat Genet", "raw counts; pseudobulked by src/pseudobulk_sc.py"),
    ("jurkat_raw_singlecell.h5ad", GEO + "GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad", "GEO GSE264667", "Nadig et al. 2025", "raw counts; pseudobulked by src/pseudobulk_sc.py"),
    ("string_links.txt.gz", "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz", "STRING v12.0", "Szklarczyk et al. 2023 NAR", "CC BY 4.0"),
    ("string_info.txt.gz", "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", "STRING v12.0", "Szklarczyk et al. 2023", "CC BY 4.0"),
    ("hallmark.txt", "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020", "MSigDB Hallmark (Enrichr 2020 library)", "Liberzon et al. 2015 Cell Syst", "gene sets"),
    ("cellage.zip", "https://genomics.senescence.info/cells/cellAge.zip", "CellAge (HAGR)", "Avelar et al. 2020 Genome Biol", "held-out evidence"),
    ("genage_models.zip", "https://genomics.senescence.info/genes/models_genes.zip", "GenAge model organisms (HAGR)", "de Magalhaes et al. 2024 NAR", "held-out evidence"),
    ("genage_human.zip", "https://genomics.senescence.info/genes/human_genes.zip", "GenAge human (HAGR)", "de Magalhaes et al. 2024", "not used in tests"),
    ("oncokb_cancer_genes.tsv", "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt", "OncoKB cancer gene list", "Chakravarty et al. 2017 JCO PO", "annotation only; check OncoKB terms before redistribution"),
    ("alliance_orthology.tsv.gz", "https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz", "Alliance of Genome Resources", "Alliance Consortium 2024 Genetics", "orthology"),
]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()


with open(r"D:\Study\RejuvBench\DATA_MANIFEST.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["local_file", "url", "accession_or_version", "citation", "notes", "access_date", "bytes", "sha256"])
    for name, url, acc, cit, note in ITEMS:
        p = os.path.join(RAW, name)
        w.writerow([name, url, acc, cit, note, ACCESS, os.path.getsize(p) if os.path.exists(p) else "", sha(p) if os.path.exists(p) else "missing"])
        print("ok", name, flush=True)
