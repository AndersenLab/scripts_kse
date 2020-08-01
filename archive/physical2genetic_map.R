# Script to find the genetic position of QTL in RIAILs
# Can be used to find the expected recombination frequency between N2/CB

library(tidyverse)
library(linkagemapping)

# function to estimate the genetic position from a physical genomic position
# input: chr = chromosome ("I", "II", "III", "IV", "V", or "X")
# input: pos = basepair position
p2gmap <- function(chr, pos) {
    # read gff
    gff <- readr::read_tsv("~/Downloads/c_elegans.current.genetic_limits.gff2", col_names = F, col_types = cols())
    
    # add column names
    colnames(gff) <- c("chrom", "pmap", "gmap_span", "pmap1", "pmap2", "test1", "test2", "test3", "info")
    
    # tidy data
    test <- gff %>%
        dplyr::select(chrom, pmap1, pmap2, info) %>%
        tidyr::separate(info, c("gene", "gmap", "note"), sep = ";") %>%
        dplyr::mutate(gene = stringr::str_split_fixed(gene, "GMap ", 2)[,2]) %>%
        dplyr::mutate(gmap = stringr::str_split_fixed(gmap, '"', 3)[,2]) %>%
        dplyr::select(-note) %>%
        dplyr::mutate(gmap = stringr::str_split_fixed(gmap, " cM", 2)[,1]) %>%
        dplyr::mutate(gmap = as.numeric(gmap))
    
    # find genetic position with physical position closest to marker
    test2 <- test %>%
        dplyr::filter(chrom == chr, pmap1 >= pos - 100000, pmap1 <= pos + 100000) %>%
        tidyr::gather(gene2, pos2, pmap1:pmap2) %>%
        dplyr::mutate(dist = abs(pos2 - pos)) %>%
        dplyr::arrange(dist)
    
    return(test2$gmap[1])
}

# function to estimate the centimorgan distance between two physical positions
# input: region in format of `chr:pos1-pos2`
genetic_distance <- function(region) {
    # parse out chromosome, pos1, and pos2
    chr <- stringr::str_split_fixed(region, ":", 2)[,1]
    positions <- stringr::str_split_fixed(region, ":", 2)[,2]
    pos1 <- as.numeric(stringr::str_split_fixed(positions, "-", 2)[,1])
    pos2 <- as.numeric(stringr::str_split_fixed(positions, "-", 2)[,2])
    
    return(glue::glue("{abs(p2gmap(chr, pos1) - p2gmap(chr, pos2))} cM"))

}

genetic_distance("IV:3841160-5110734")

