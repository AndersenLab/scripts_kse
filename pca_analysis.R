library(tidyverse)
library(linkagemapping)

# Calculate principal components from phenotype dataframe with replicates (i.e. NILs)
# pheno - dataframe of phenotypes
# returns list of (1) pca object and (2) pca pheno
calc_pc_reps <- function(pheno) {
    # calculate PC for linkage
    pc_traits <- pheno %>%
        dplyr::mutate(drugtrait = paste0(condition, ".", trait)) %>%
        dplyr::ungroup()%>%
        dplyr::mutate(well = paste(assay, round, plate, row, col, strain, sep = "-")) %>%
        dplyr::select(well, drugtrait, phenotype)%>%
        unique() %>%
        tidyr::spread(drugtrait, phenotype) %>%
        na.omit()
    
    # keep strains as rownames and remove strain
    row.names(pc_traits) <- pc_traits$well
    pc_traits <- pc_traits %>%
        dplyr::select(-well)
    
    # scale traits
    scales_pc_traits <- as.data.frame(scale(pc_traits))
    
    # create PCs
    pca_obj <- princomp(scales_pc_traits)
    
    # figure out how many PCs to keep that will explain > 90% of the variance
    # pull the total variance explained with each PC and call it "drug.cumsum"
    cumsums <- t(as.matrix(cumsum(pca_obj$sdev^2/sum(pca_obj$sdev^2))))
    cumsums <- as.data.frame(cumsums) %>%
        tidyr::gather(comp, var) %>%
        dplyr::filter(var > 0.9) 
    
    # the first component in the dataframe is the first one that goes over 90%, so we want to keep it and everything else below it
    keep <- as.numeric(stringr::str_split_fixed(cumsums$comp[1], "Comp.", 2)[,2])
    
    # keep only those PCs phenotypes for all RIAIL strains
    colnames(pca_obj$scores) <- paste(pheno$condition[1], colnames(pca_obj$scores), sep = "_")
    
    PCpheno <- data.frame(pca_obj$scores[,1:keep]) %>%
        dplyr::mutate(well = rownames(.)) %>%
        tidyr::separate(well, into = c("assay", "round", "plate", "row", "col", "strain"), sep = "-") %>%
        tidyr::gather(trait, phenotype, -c(assay:strain)) %>%
        dplyr::mutate(trait = gsub("Comp.", "PC", trait),
                      phenotype = phenotype) %>%
        dplyr::mutate(condition = stringr::str_split_fixed(trait, "_", 2)[,1],
                      trait = stringr::str_split_fixed(trait, "_", 2)[,2]) %>%
        dplyr::mutate(phenotype = as.numeric(phenotype),
                      round = as.numeric(round),
                      plate = as.numeric(plate))
    
    return(list(pca_obj, PCpheno))
}

# Calculate principal components from phenotype dataframe with no replicates (i.e. RIAILs)
# pheno - dataframe of phenotypes
# returns list of (1) pca object and (2) pca pheno
calc_pc_noreps <- function(pheno) {
    # calculate PC for linkage
    pc_traits <- pheno %>%
        dplyr::mutate(drugtrait = paste0(condition, ".", trait)) %>%
        dplyr::ungroup()%>%
        dplyr::select(strain, drugtrait, phenotype)%>%
        unique() %>%
        tidyr::spread(drugtrait, phenotype) %>%
        na.omit()
    
    # keep strains as rownames and remove strain
    row.names(pc_traits) <- pc_traits$strain
    pc_traits <- pc_traits %>%
        dplyr::select(-strain)
    
    # scale traits
    scales_pc_traits <- as.data.frame(scale(pc_traits))
    
    # create PCs
    pca_obj <- princomp(scales_pc_traits)
    
    # figure out how many PCs to keep that will explain > 90% of the variance
    # pull the total variance explained with each PC and call it "drug.cumsum"
    cumsums <- t(as.matrix(cumsum(pca_obj$sdev^2/sum(pca_obj$sdev^2))))
    cumsums <- as.data.frame(cumsums) %>%
        tidyr::gather(comp, var) %>%
        dplyr::filter(var > 0.9) 
    
    # the first component in the dataframe is the first one that goes over 90%, so we want to keep it and everything else below it
    keep <- as.numeric(stringr::str_split_fixed(cumsums$comp[1], "Comp.", 2)[,2])
    
    # keep only those PCs phenotypes for all RIAIL strains
    colnames(pca_obj$scores) <- paste(pheno$condition[1], colnames(pca_obj$scores), sep = "_")
    
    PCpheno <- data.frame(pca_obj$scores[,1:keep]) %>%
        dplyr::mutate(strain = rownames(.)) %>%
        tidyr::gather(trait, phenotype, -strain) %>%
        dplyr::mutate(trait = gsub("Comp.", "PC", trait),
                      phenotype = phenotype) %>%
        tidyr::separate(trait, into = c("condition", "trait"), sep = "_") %>%
        dplyr::mutate(phenotype = as.numeric(phenotype))
    
    return(list(pca_obj, PCpheno))
}

# Predicts PCA phenotypes given loadings from another pca object (use for applying pca from linkage to NILs)
# pheno - dataframe of phenotypes
# pca_oject - pca object, output [[1]] from calc_pc_reps or calc_pc_noreps
# keep - how many PCs to keep? default is number of traits in pheno
predict_pc <- function(pheno, pca_obj, keep = length(unique(pheno$trait))) {
    # clean dataframe
    pc_traits <- pheno %>%
        dplyr::mutate(drugtrait = paste0(condition, ".", trait)) %>%
        dplyr::ungroup()%>%
        dplyr::mutate(well = paste(round, assay, plate, row, col, strain, sep = "-")) %>%
        dplyr::select(well, drugtrait, phenotype)%>%
        unique() %>%
        tidyr::spread(drugtrait, phenotype) %>%
        na.omit()
    
    # keep strains as rownames and remove strain
    row.names(pc_traits) <- pc_traits$well
    pc_traits <- pc_traits %>%
        dplyr::select(-well)
    
    # scale traits
    scales_pc_traits <- as.data.frame(scale(pc_traits))
    
    # predict PCA based on pca_obj
    pcapredict <- data.frame(predict(pca_obj, scales_pc_traits))
    colnames(pcapredict) <- paste(pheno$condition[1], colnames(pcapredict), sep = "_")
    
    pcaout <- pcapredict[,1:keep] %>%
        dplyr::mutate(well = rownames(.)) %>%
        tidyr::separate(well, into = c("round", "assay", "plate", "row", "col", "strain"), by = "-") %>%
        tidyr::gather(trait, phenotype, -c(round:strain)) %>%
        dplyr::mutate(trait = gsub("Comp.", "PC", trait),
                      phenotype = phenotype) %>%
        tidyr::separate(trait, into = c("condition", "trait"), sep = "_") %>%
        dplyr::mutate(phenotype = as.numeric(phenotype))
    
    return(pcaout)
}
