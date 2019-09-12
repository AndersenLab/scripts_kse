# Script to find NILs in your area of interest

# The idea is to input a region (larger than your QTL) and the script will be able to tell you if we already have NILs here.
# INPUT: chromosome region ("IV:5000000-8000000")
# OUTPUT: list of [1] NIL genotypes on chrom of interest with region bounded.
#                 [2] NIL genotypes on all chrom
#                 [3] dataframe of NIL genotype

input <- "X:10000000-15000000"
chrIIInils <- paste0("ECA", c(832:877))

findNILs <- function(input) {
    
    # Load libraries
    library(tidyverse)
    library(linkagemapping)
    
    # Source nil genotype file
    source("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/scripts_kse/NIL_genotype_plots.R")
    
    # Get chromosome region
    chr <- stringr::str_split_fixed(input, ":", 2)[1]
    range <- stringr::str_split_fixed(input, ":", 2)[2]
    left_pos <- as.numeric(stringr::str_split_fixed(range, "-", 2)[1])
    right_pos <- as.numeric(stringr::str_split_fixed(range, "-", 2)[2])
    
    # Filter nil genotype data:
         # - filter to chromosome of interest
         # - filter to keep only ECA strains
         # - filter the length of a segment to be > 10,000 bp
         # - filter to keep only strains with N2 and CB on chrom
    nils <- nilgeno %>%
        dplyr::filter(!(sample %in% c("N2", "CB4856"))) %>%
        dplyr::mutate(length = end - start) %>%
        dplyr::group_by(sample, gt) %>%
        dplyr::mutate(background = sum(length)) %>%
        tidyr::spread(gt_name, background) %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(CB4856 = mean(CB4856, na.rm = T),
                      N2 = mean(N2, na.rm = T),
                      RIAIL = N2/CB4856) %>%
        dplyr::filter(RIAIL > 5 | RIAIL < 0.25) %>% #0.167 is one chromosome (CSS) but better safer than sorry probably.
        dplyr::group_by(sample) %>%
        dplyr::mutate(backgroundtype = ifelse(N2 > CB4856, "N2", ifelse(CB4856 > N2, "CB4856", NA))) %>%
        dplyr::filter(chrom == chr, grepl("ECA", sample)) %>%
        dplyr::filter(length > 100000) %>% # what should I use as the threshold? smallest NIL in chrIII nils is 241,230 bp
        dplyr::group_by(sample) %>%
        dplyr::mutate(num_genos = sum(plyr:::nunique(gt))) %>%
        dplyr::filter(num_genos == 2, sample != "ECA257") # this is the kammenga lab N2 

    # Separate into N2 and CB nils, then make sure the NIL region is within the boundary
    cbnils <- nils %>%
        dplyr::filter(backgroundtype == "N2", gt == 2, start < right_pos, end > left_pos) #1 = N2, 2 = CB
    
    n2nils <- nils %>%
        dplyr::filter(backgroundtype == "CB4856", gt == 1, start < right_pos, end > left_pos)
    
    allnils <- rbind(cbnils, n2nils)
    
    # outputs
    plot1 <- nil_plot(unique(allnils$sample), chr = chr, all.chr = F, left.bound = left_pos, right.bound = right_pos, left.cb = left_pos, left.n2 = left_pos)[[1]] +
        geom_vline(xintercept = c(left_pos/1e6, right_pos/1e6))
    
    plot2 <- nil_plot(unique(allnils$sample), chr = chr, all.chr = T, left.bound = left_pos, right.bound = right_pos, left.cb = left_pos, left.n2 = left_pos)[[1]]
    
    nildf <- nilgeno %>%
        dplyr::filter(sample %in% allnils$sample)
    
    return(list(plot1, plot2, nildf))
}

findNILs("V:7e6-12e6")[[1]]
