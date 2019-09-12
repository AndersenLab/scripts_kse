# updated on 8.27.18 to pull from wbgenes dataframe (using refFlat.txt) instead of gene_info dataframe
# updated on 11.20.18 to use cegwas2::query_vcf() instead of cegwas::snpeff()

library(cegwas2)
library(dplyr)

#Define function to narrow QTL for given region of interest
qtl_narrow <- function(query, sev = c("MODERATE", "HIGH")) {
  #load gene info
  # load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/wbgenes.Rda")
    goterm <- read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/GOannotations.tsv", col_names = FALSE)
    names(goterm) <- c("GO_annotation_extension", "GO_annotation_qualifier", "Data_set", "Data_set_URL", "Evidence_w_text", "Ontology_term_identifier",
                     "Ontology_term_name", "Ontology_term_description", "Subject_wormbase_ID")
    goterm <- goterm %>%
        dplyr::select(Ontology_term_identifier:Subject_wormbase_ID)
  
  # snpeff <- cegwas::snpeff(query)
  snpeff <- cegwas2::query_vcf(query, impact = sev, samples = "CB4856")
  print(paste0("Total genes in interval: ", length(unique(snpeff$gene_id))))
  
  #Only keep genes that have difference in N2/CB
  snpf <- snpeff  %>%
      dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF, "REF", "ALT"))
  
  # print(paste0("Total Coding genes with differences in N2/CB in interval: ", length(unique(snpf$gene_id))))
  
  #Add GO Terms and descriptions
  df <- snpf %>%
    select(CHROM, POS, strain = SAMPLE, REF, ALT, a1, a2, GT, query, allele, effect, impact, gene_name, gene_id, feature_type, 
           feature_id, transcript_biotype, nt_change, aa_change) %>%
    left_join(goterm, by = c("gene_id" = "Subject_wormbase_ID"))
  
  return(df)
}
