library(tidyverse)

region <- "V:4227886-5800599"

# function
chrom <- stringr::str_split_fixed(region, ":", 2)[,1]
start <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,1])
end <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,2])

load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_ref_flat.Rda")
genes <- gene_ref_flat %>%
    dplyr::filter(chr == chrom,
                  txstart < end,
                  txend > start)
length(unique(genes$wbgene))
