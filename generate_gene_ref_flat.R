library(tidyverse)
library(data.table)

# get elegans_transcripts_WS280.bed.gz from wormbase:
# curl ftp://ftp.wormbase.org/pub/wormbase/releases/current-production-release/MULTI_SPECIES/hub/elegans/elegans_genes_WS280.bb > elegans_genes_WS280.bb
# # install bigbedtobed
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
# chmod a+x bigBedToBed
# module load bedtools
# # convert to bed
# ./bigBedToBed elegans_genes_WS280.bb tmp.bed
# sortBed -i tmp.bed > elegans_transcripts_WS280.bed
# bgzip -f elegans_transcripts_WS280.bed

gene_ref_flat <- data.table::fread("~/Downloads/elegans_transcripts_WS280.bed.gz", col.names = c("chr", "txstart", "txend", "gene", "score", "strand",
                                                                                         "codingstart", "codingend", "unknown", "numexons", 
                                                                                         "exonlengths", "exonstarts", "gene2", "test", "test2", 
                                                                                         "test3", "biotype", "wbgene", "gene_name", "biotype2")) %>%
    dplyr::mutate(exonstarts = gsub(',$', '', exonstarts),
                  exonlengths = gsub(',$', '', exonlengths)) %>%
    tidyr::separate_rows(c(exonstarts, exonlengths)) %>%
    dplyr::mutate(exonstarts = as.numeric(txstart) + as.numeric(exonstarts),
                  exonends = as.numeric(exonstarts) + as.numeric(exonlengths)) %>%
    dplyr::select(gene, chr, strand, txstart, txend, codingstart, codingend, numexons, exonstarts, exonends, wbgene, gene_name, biotype) %>%
    dplyr::group_by(gene, chr, strand, txstart, txend, codingstart, codingend, numexons, wbgene, gene_name, biotype) %>%
    dplyr::summarize(exonstarts = paste(exonstarts, collapse = ","),
                     exonends = paste(exonends, collapse = ",")) %>%
    dplyr::mutate(type = "Transcript")

save(gene_ref_flat, file = "~/Downloads/gene_ref_flat.Rda")
