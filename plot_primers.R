# Script to plot NIL genotype and primer locations
# you need to install googlesheets
# install.packages('gsheet')

library(tidyverse)
library(gsheet)

# primers
primers <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1LJFZZ4dZm9KnoTpZqwCkQiXAKuD4Q4IGebNu1OKSnXM/edit#gid=0') %>%
    dplyr::mutate(pos = `Physical Position` / 1e6,
                  pair = paste0(`Left primer`, "-", `Right primer`),
                  pair2 = paste0(`Right primer`, "-", `Left primer`))

# load genotypes
#load nilgenos
load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/nilgeno.Rda")

# function
# primer pairs - pair of primers in the format "oECAXXX-oECAXXX", can be multiple
# NIL - ECAxxx strain name to plot genotype (N2 is default)
# RIL = false. If true, will plot RILs instead of NILs
plot_primers <- function(primer_pairs, NIL = "N2", RIL = F, chr = c("I", "II", "III", "IV", "V", "X")) {
    
    # get primer positions
    primer_positions <- primers %>%
        dplyr::filter(pair %in% primer_pairs | pair2 %in% primer_pairs)
    
    if(RIL == T) {
        source("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/scripts_kse/RIL_genotype_plots.R")
        df <- rilgeno
    } else {
        df <- nilgeno
    }
    
    # plot NIL genotype
    df %>%
        dplyr::filter(sample %in% NIL) %>%
        dplyr::filter(chrom %in% chr) %>%
        dplyr::mutate(gt = as.character(gt)) %>%
        ggplot(.)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("1"="orange","2"="blue"))+
            theme_bw() +
            theme(axis.text.x = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_text(size=12, face="bold", color="black"),
                  axis.title.x = element_text(size=14, face="bold", color="black"),
                  axis.title.y = element_text(size=14, face="bold", color="black"),
                  strip.text = element_text(size = 12, face = "bold", color = "black"),
                  plot.title = element_text(size=24, face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            labs(x = "Genomic Position (Mb)", y = "") +
            geom_vline(data = data.frame(pos = primer_positions$pos, chrom = primer_positions$Chromosome), aes(xintercept = pos)) +
            geom_text(data = data.frame(pos = primer_positions$pos, 
                                        chrom = primer_positions$Chromosome, 
                                        strain = NIL[length(NIL)],
                                        primers = primer_positions$pair), aes(x = pos,y = strain, label = primers), 
                      angle = 90,
                      hjust = -0.1,
                      vjust = -0.5)
}


########## TEST ###########
# primer_pairs <- c("oECA1148-oECA1147", "oECA1408-oECA1409", "oECA1141-oECA1142")
# NIL <- "ECA1065"
# plot_primers(primer_pairs, NIL)

# plot_primers(primer_pairs, c("N2", "ECA1065"))
