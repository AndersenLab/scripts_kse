# make dataframe for all N2/CB NIL genotypes
library(tidyverse)
library(linkagemapping)

#load nilgenos
nilgeno1 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20170829/vcf/gt_hmm_fill.tsv") %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches)
nilgeno2 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20170128/hmm/gt_hmm_fill.tsv") %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches)
nilgeno3 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20171212/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno4 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20180409/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno5 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20180626/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno6 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20181016/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno7 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20191016/hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)
nilgeno8 <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/NIL/20200323/NIL-N2-CB4856-20200326//hmm/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)

nilgeno9 <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Loraina/DrugsandToxins/NILS/gt_hmm_fill.tsv") %>%
    dplyr::select(-rle)

N2_CB_geno <- read.csv("~/Dropbox/AndersenLab/LabFolders/Katie/projects//HTA_sorter/N2_CB_geno.csv")  %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches)
nilgeno <- rbind(nilgeno1, nilgeno2, nilgeno3, nilgeno4, nilgeno5, nilgeno6, nilgeno7, nilgeno8,nilgeno9, N2_CB_geno)

save(nilgeno, file = "~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/nilgeno.Rda")

# update RIL genotype dataframe
data("N2xCB4856cross")

set <- N2xCB4856cross$pheno
rilgeno <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/RIL/RIL-2017-03-27/vcf/gt_hmm_fill.tsv") %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches) %>%
    dplyr::left_join(set, by = c("sample" = "strain"))

save(rilgeno, file = "~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/rilgeno.Rda")
