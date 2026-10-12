$D = "D:\Study\RejuvBench\data\raw"
$files = [ordered]@{
 "GSE113957_fpkm.txt.gz" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113957/suppl/GSE113957_fpkm.txt.gz"
 "GSE113957_series_matrix.txt.gz" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113957/matrix/GSE113957_series_matrix.txt.gz"
 "hallmark.txt" = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
 "genage_models.zip" = "https://genomics.senescence.info/genes/models_genes.zip"
 "genage_human.zip" = "https://genomics.senescence.info/genes/human_genes.zip"
 "cellage.zip" = "https://genomics.senescence.info/cells/cellAge.zip"
 "oncokb_cancer_genes.tsv" = "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt"
 "gtex_subject.txt" = "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
 "gtex_sample.txt" = "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
 "K562_essential_normalized_bulk_01.h5ad" = "https://ndownloader.figshare.com/files/35780870"
 "rpe1_normalized_bulk_01.h5ad" = "https://ndownloader.figshare.com/files/35775512"
 "string_info.txt.gz" = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
 "string_links.txt.gz" = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
 "gtex_fibroblasts_tpm.gct.gz" = "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cells_cultured_fibroblasts.gct.gz"
 "gtex_skin_ns_tpm.gct.gz" = "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_skin_not_sun_exposed_suprapubic.gct.gz"
 "K562_gwps_normalized_bulk_01.h5ad" = "https://ndownloader.figshare.com/files/35773217"
 "hepg2_raw_singlecell.h5ad" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad"
 "jurkat_raw_singlecell.h5ad" = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad"
}
foreach ($k in $files.Keys) {
  $p = Join-Path $D $k
  if ((Test-Path $p) -and ((Get-Item $p).Length -gt 0)) { "have $k"; continue }
  $t = Get-Date
  curl.exe --doh-url https://1.1.1.1/dns-query -sSL --retry 5 --retry-delay 5 -C - -A "Mozilla/5.0" -o "$p.part" $files[$k]
  if ($LASTEXITCODE -eq 0) { Move-Item -Force "$p.part" $p; "ok $k $((Get-Item $p).Length) $([int]((Get-Date)-$t).TotalSeconds)s" } else { "FAIL $k code $LASTEXITCODE" }
}
"DONE"

