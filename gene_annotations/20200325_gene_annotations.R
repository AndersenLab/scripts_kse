# Created on 20200325 using WS275 annotations from wormbase FTP
# not sure where the GOannotations.tsv came from (from wormbase GO somewhere...) but only used for go descriptions, not gene annotations

library(tidyverse)

# go annotations
go <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_annotations/c_elegans.PRJNA13758.current.go_annotations.gaf", 
                      skip = 3, col_names = F) %>%
    dplyr::select(X2, X3, X5:X9, X11:X15)

# for now probably only need the gene and go term, then need to figure out what these go terms are...
new_go <- go %>%
    dplyr::select(wbgene = X2, gene_name = X3, go_term = X5) %>%
    dplyr::distinct()

# old go annotations from somewhere else on wormbase... do they match?
oldgo <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_annotations/GOannotations.tsv", 
                         col_names = F) %>%
    dplyr::select(go_term = X6, go_name = X7, go_description = X8) %>%
    dplyr::distinct()

# use oldgo for descriptions (has 6503 of the 6733 in the new go... pretty good)
go_terms <- new_go %>%
    dplyr::left_join(oldgo, by = "go_term")

# functional annotations
functional <- readr::read_delim("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_annotations/c_elegans.PRJNA13758.current.functional_descriptions.txt", 
                                skip = 3, delim = "\n", col_names = F)

# make a list of the rows for each gene starting
gene_starts <- grep("WBGene", functional$X1)

# paste all the rows for each gene into one row
new_functional <- data.frame(gene = rep(NA, length(gene_starts)))
for(i in 1:length(gene_starts)) {
    # get where the gene starts and ends
    k <- gene_starts[i]
    j <- gene_starts[i+1]-2
    
    # how to handle the last gene
    if(is.na(j)) {
        j <- nrow(functional)
    }
    
    # paste all rows together
    new_functional$gene[i] <- paste(functional$X1[k:j], collapse = " ")
}

# tidy!
tidy_functional <- new_functional %>%
    tidyr::separate(gene, into = c("wbgene", "gene_name", "gene_id"), sep = "\t") %>%
    tidyr::separate(gene_id, into = c("gene_id", "concise_description"), sep = " Concise description: ") %>%
    tidyr::separate(concise_description, into = c("concise_description", "automated_description"), sep = " Automated description: ") %>%
    tidyr::separate(automated_description, into = c("automated_description", "gene_class_description"), sep = " Gene class description: ")

gene_annotations <- tidy_functional %>%
    dplyr::full_join(go_terms, by = c("wbgene", "gene_name"))

save(gene_annotations, file = "~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_annotations.Rda")
