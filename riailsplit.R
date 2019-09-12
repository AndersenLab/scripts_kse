library(dplyr)
library(linkagemapping)


load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/docetaxel/data/20170425_docetaxel_css_nil_gwer.RData")
q75EXT <- docetaxelRIAILsregressed %>%
  filter(trait == "q75.EXT")

#Label RIAILs as N2 or CB with phenotype listed
extract_genotype=function(cross){
  
  # Pull out the genotypes into snp x strain matrix
  genomat <- qtl::pull.geno(cross)
  class(genomat) <- "numeric"
  
  # Handle genotype encodings with heterozygous individuals
  # (encoded as 1, 2, 3) and without (encoded 1, 2)
  if(max(genomat, na.rm = TRUE) == 3) {
    genomat <- genomat - 2
  } else {
    genomat <- (genomat * 2) - 3
  }
  return(genomat)
}
getRIAILs <- function(cross, map, parent = "N2xCB4856") {
  uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
  colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
  pheno <- cross$pheno %>% dplyr::select_(map$trait[1])
  strains <- cross$pheno %>% dplyr::select(strain)
  pheno <- cbind(strains, pheno)
  geno_old <- data.frame(extract_genotype(cross)) %>% 
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
    data.frame(., pheno)
  colnames(geno)[1:(ncol(geno) - 1)] <- sapply(colnames(geno)[1:(ncol(geno) - 1)], 
                                               function(marker) {paste(unlist(peaks[peaks$marker == gsub("\\.", "-", marker),
                                                                                    c("chr", "pos")]), collapse = ":")})
  colnames(geno)[ncol(geno)] <- "pheno"
  colnames(geno)[2] <- "strain"
  colnames(geno)[1] <- "genotype"
  split <- geno
  # split <- geno %>% tidyr::gather(marker, genotype, -c(pheno, strain))
  split$genotype <- sapply(split$genotype, function(x) {
    if (x == -1) {
        "N2"
    }
    else {
      "CB4856"
    }
  })
  
  split <- split %>%
    na.omit()
  return(split)
}


# #Get RIAIL genotype data
# riails <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/RIL/RIL-2017-03-27/hmm/gt_hmm.tsv")




