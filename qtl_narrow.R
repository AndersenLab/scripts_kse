# updated on 8.27.18 to pull from wbgenes dataframe (using refFlat.txt) instead of gene_info dataframe
# updated on 11.20.18 to use cegwas2::query_vcf() instead of cegwas::snpeff()

library(cegwas2)
library(dplyr)


load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_descriptions_WS273.Rda")

#Define function to narrow QTL for given region of interest
qtl_narrow <- function(query, sev = c("MODERATE", "HIGH", "LOW", "MODIFIER")) {
  
  # snpeff <- cegwas::snpeff(query)
  snpeff <- cegwas2::query_vcf(query, impact = sev, samples = "CB4856")
  print(paste0("Total genes in interval: ", length(unique(snpeff$gene_id))))
  
  #Only keep genes that have difference in N2/CB
  snpf <- snpeff  %>%
      dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF, "REF", "ALT"))
  
  # print(paste0("Total Coding genes with differences in N2/CB in interval: ", length(unique(snpf$gene_id))))
  
  #Add GO Terms and descriptions
  df <- snpf %>%
    dplyr::select(CHROM, POS, strain = SAMPLE, REF, ALT, a1, a2, GT, query, allele, effect, impact, gene_name, gene_id, feature_type, 
           feature_id, transcript_biotype, nt_change, aa_change) %>%
    dplyr::left_join(gene_descriptions, by = c("gene_id" = "wbgene"))
  
  return(df)
}
