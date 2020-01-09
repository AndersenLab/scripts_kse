# fosmid update
library(tidyverse)
library(linkagemapping)
data("AllCBfosmids")
data("AllN2fosmids")

# load data
load(CBfosmids_pacbio, file = "~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/CBfosmids_pacbio.Rda")

# function to plot all fosmids within region in CB (mapped to CB pacbio genome)
# chrom - chromosome
# startpos - start position
# endpos - end position
plot_fosmids <- function(chrom, startpos, endpos) {
    # which fosmids in your region?
    fosmid_data <- CBfosmids_pacbio %>%
        dplyr::filter(chr == chrom, start > startpos, end < endpos) %>%
        dplyr::arrange(start, end)
    
    # plot
    plot <- CBfosmids_pacbio %>%
        dplyr::filter(chr == chrom, start > startpos, end < endpos) %>%
        dplyr::arrange(start, end) %>%
        tibble::rownames_to_column(var = "row") %>%
        dplyr::mutate(row = as.numeric(row)) %>%
        ggplot(.) +
        geom_segment(aes(x = start/1e6, xend = end/1e6, y = row, yend = row)) +
        geom_text(aes(x = end/1e6, y = row + 0.3, label = id), size = 2) +
        theme_bw() +
        labs(x = "Genomic Position", y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank()) +
        facet_grid(~chr, scales = "free")
    
    return(list(fosmid_data, plot))
}

# example plotting
plot_fosmids("V", 10e6, 11e6)[[2]]

# plot only necessary (minimal) fosmids in region (overlapping ends)
# fosmid_df - dataframe comtaining all fosmids (pacbio or old from linkage mapping)
# chrom - chromosome
# startpos - start position
# endpos - end position
plot_minimal_fosmids <- function(fosmid_df, chrom, startpos, endpos) {
    # if the old dataframe, change "clone" to "id"
    if("clone" %in% names(fosmid_df)) {
        fosmid_df <- fosmid_df %>%
            dplyr::mutate(id = clone)
    }
    
    # which fosmids in your region?
    all_fosmids <- fosmid_df %>%
        dplyr::filter(chr == chrom, start > startpos, end < endpos) %>%
        dplyr::arrange(start, end) %>%
        dplyr::mutate(keep = NA)
    
    # set up for while loop
    remaining_df <- all_fosmids
    finaldf <- all_fosmids[1,] %>%
        dplyr::mutate(keep = T)
    newrow <- 1
    
    # while there is still fosmids remaining, look for overlaps
    while(nrow(remaining_df) > 1) {
        # find all fosmids within bounds and keep first and last (overlapping ends)
        df <- remaining_df %>%
            dplyr::filter(start < all_fosmids$end[newrow]) %>%
            ungroup() %>%
            dplyr::mutate(keep = ifelse(start %in% c(dplyr::first(start), dplyr::last(start)), T, F))
        
        # if there is only 1 row, there are no overlaps, so move on to the next row
        if(nrow(df) == 1) {
            # assign to final dataframe
            finaldf <- dplyr::bind_rows(finaldf, df)
            
            # remove from remaining
            remaining_df <- remaining_df %>%
                dplyr::filter(!(id %in% df$id))
            
            # move to the last fosmid kept and make it the first
            newrow <- which(all_fosmids$id == dplyr::last(df$id)) + 1
        } else {
            # assign to final dataframe
            finaldf <- dplyr::bind_rows(finaldf, data.frame(df[2:nrow(df),]))
            
            # move to the last fosmid kept and make it the first
            newrow <- which(all_fosmids$id == dplyr::last(df$id))
            
            # remove non-important fosmids from remaining_df
            df <- df %>%
                dplyr::filter(id != dplyr::last(id))
            remaining_df <- remaining_df %>%
                dplyr::filter(!(id %in% df$id))
        }
        
    }
    
    # plot
    plot <- finaldf %>%
        dplyr::filter(chr == chrom, start > startpos, end < endpos, keep == T) %>%
        dplyr::arrange(start, end) %>%
        tibble::rownames_to_column(var = "row") %>%
        dplyr::mutate(row = as.numeric(row)) %>%
        ggplot(.) +
        geom_segment(aes(x = start/1e6, xend = end/1e6, y = row, yend = row)) +
        geom_text(aes(x = end/1e6, y = row + 0.3, label = id), size = 2) +
        theme_bw() +
        labs(x = "Genomic Position", y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank()) +
        facet_grid(~chr, scales = "free")
    
    return(list(finaldf, plot))
}

# example plotting
plot_minimal_fosmids(CBfosmids_pacbio, "V", 10e6, 11e6)[[2]]
plot_minimal_fosmids(AllCBfosmids, "V", 10e6, 11e6)[[2]]
plot_minimal_fosmids(CBfosmids_pacbio, "III", 1, 800000)[[2]]

