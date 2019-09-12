library(tidyverse)
library(linkagemapping)
data("N2xCB4856cross")

set <- N2xCB4856cross$pheno

# plot RIAIL genotypes
rilgeno <- data.table::fread("~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/RIL/RIL-2017-03-27/vcf/gt_hmm_fill.tsv") %>%
    dplyr::mutate(gt_name = ifelse(gt == 1, "N2", "CB4856")) %>%
    dplyr::select(chrom, start, end, sample, gt, gt_name, supporting_sites, sites, DP, switches) %>%
    dplyr::left_join(set, by = c("sample" = "strain"))


# Function to plot RIL genotypes based on WGS data from 20170327 - 
# includes 1200 strains, all of set 1 and set 2 and "E" strains (not set3?)
# default plots all chromosomes, you can provide vector of chromosomes to plot
# must supply either a vector of strains or a strain set (1 or 2)
# can add strain labels with strain_label = T
plot_rilgeno <- function(strains = NA, chr = c("I", "II", "III", "IV", "V", "X"), strainset = NA, theme_size = 12, strain_label = F) {
    # if set != NA, plot all RILs in the set
    if(!is.na(strainset)) {
        if(strainset %in% c(1, 2)) {
            df <- rilgeno %>%
                dplyr::filter(chrom %in% chr, set == strainset)
        } else {
            stop("Error. Must supply either vector of strains or strain set (1 or 2)")
        }
    } else {
        if(length(strains) == 1 && is.na(strains)) {
            stop("Error. Must supply either vector of strains or strain set (1 or 2)")
        } else {
            df <- rilgeno %>%
                dplyr::filter(chrom %in% chr, sample %in% strains)
        }
    }
    
    # add or remove strain labels, depending on flag
    if(strain_label == T) {
        df %>%
            ggplot(.)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
            theme_bw(theme_size) +
            theme(axis.text = element_text(color="black"),
                  axis.title = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  legend.position = "none",
                  panel.grid = element_blank()) +
            labs(x = "Genomic Position (Mb)", y = "")
    } else {
        df %>%
            ggplot(.)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
            theme_bw(theme_size) +
            theme(axis.text.x = element_text(color="black"),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  legend.position = "none",
                  panel.grid = element_blank()) +
            labs(x = "Genomic Position (Mb)", y = "")
    }
    
}
