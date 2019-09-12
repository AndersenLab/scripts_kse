# Script to run burden test
library(tidyverse)

#### step 1: prepare .ped file (phenotypes) ####
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/zinc/data/gwas/zincGWASregressed.Rda")
traits <- zincGWASregressed %>%
    dplyr::select(condition, trait, value = phenotype, strain) %>%
    dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2, trait = paste0(condition, ".", trait))%>%
    dplyr::filter(trait == "zinc.q90.EXT") %>%
    tidyr::spread(trait,value)%>%
    dplyr::select(-strain, -condition)

# take care of NA values
traits[is.na(traits)] <- -9

# save file
write.table(traits,"~/Downloads/burden_test/zinc_q90EXT.ped", quote = F, col.names = F, row.names = F)

#### step 2: run burden test on command line ####
# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno zinc_q90EXT.ped --inVcf WI.20180527.impute.vcf.gz --freqUpper .05 --freqLower 0.003 --out test_q90TOF.ped --geneFile refFlat.ws245.txt --vt price

#### step 3: analyze results ####
map <- data.table::fread("~/Downloads/burden_test/test_q90TOF.ped.VariableThresholdPrice.assoc")%>%
    dplyr::rowwise()%>%
    dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                  POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                  endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
    dplyr::mutate(size = abs(POS-endPOS)) %>%
    ungroup()%>%
    dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
    dplyr::filter(CHROM!="MtDNA",NumVar>1,size >500)

ggplot(data = map) +
    aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
    geom_point()+
    scale_color_manual(values=c("black","red"))+
    facet_grid(.~CHROM, scales = "free")+
    theme_bw()+
    theme(axis.text.x = ggplot2::element_text(size = 14),
          axis.text.y = ggplot2::element_text(size = 14),
          legend.position = "none",
          axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
          axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3))+
    labs(x = "Genomic Position (Mb)", y = "Test Statisitic")
