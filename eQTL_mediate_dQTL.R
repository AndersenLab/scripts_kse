# write a function to test a gene and qtl given 1) peak marker, 2) probe id, 3) phenodf
# peak marker in format chr:pos
# probeID from expression data from rockman
# pheno df will be in the form of strain, trait, phenotype with the strain matching RIL from set1 (where we have expression data)
library(linkagemapping)
library(mediation)
# linkagemapping::load_cross_obj("N2xCB4856cross_full")
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/eQTL_mediation/data/raw/N2xCB4856cross_full2.Rda")
expression_pheno <- read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/eQTL_mediation/data/raw/expression_phenos.tsv")

eQTL_mediate_dQTL <- function(peak, probe, phenodf, scaled = T, lm = TRUE) {
    
    # get the genotype at the peak marker
    newpeak <- gsub(":", "_", peak)
    chrom <- stringr::str_split_fixed(newpeak, "_", 2)[,1]
    geno <- data.frame(N2xCB4856cross_full2$geno[[chrom]]$data)[,newpeak]
    strains <- N2xCB4856cross_full2$pheno
    geno <- cbind(strains, geno)
    
    # get the expression pheno for that probe
    probepheno <- expression_pheno %>%
        dplyr::filter(trait == probe) %>%
        dplyr::select(strain, expression = phenotype)
    
    # if scaled = T, scale the phenotype (mean = 0, var = 1)
    if(scaled == T) {
        probepheno <- probepheno %>%
            dplyr::mutate(newpheno = (expression - mean(expression, na.rm = T)) / sd(expression, na.rm = T)) %>%
            dplyr::select(strain, expression = newpheno)
    }
    
    # merge drug phenotype and expression phenotype
    pheno <- phenodf %>%
        dplyr::left_join(probepheno, by = "strain") %>%
        dplyr::left_join(geno, by = "strain") %>%
        na.omit()
    
    if(lm) {
        ### LINEAR MODELS ###
        
        # total effect = geno estimate
        model.g <- lm(phenotype ~ geno, data = pheno)
        total <- data.frame(var = "total",
                            estimate = summary(model.g)$coef[2,1],
                            pval = summary(model.g)$coef[2,4])
        
        # direct effect = geno estimate
        model.y <- lm(phenotype ~ expression + geno, data = pheno)
        direct <- data.frame(var = "direct",
                             estimate = summary(model.y)$coef[3,1],
                             pval = summary(model.y)$coef[3,4])
        
        # mediation effect = expression estimate
        med <- data.frame(var = "med",
                          estimate = summary(model.y)$coef[2,1],
                          pval = summary(model.y)$coef[2,4])
        
        # mediation proportion = total - direct / total
        out <- rbind(total, direct, med)
    } else {
        model.m <- lm(expression ~ geno, data = pheno)
        model.y <- lm(phenotype ~ expression + geno, data = pheno)
        out <- mediation::mediate(model.m, model.y, sims = 1000, boot = T, treat = "geno", mediator = "expression")
    }
    
    return(out)
}

# function to get the summary statistics from model
summarize_model <- function(model) {
    # causal mediation effect
    acme <- data.frame(var = "ACME", 
               estimate = model$d0, 
               ci_lower = model$d0.ci[[1]], 
               ci_upper = model$d0.ci[[2]],
               prob = model$d0.p)
    
    # direct effect
    ade <- data.frame(var = "ADE", 
                       estimate = model$z0, 
                       ci_lower = model$z0.ci[[1]], 
                       ci_upper = model$z0.ci[[2]],
                       prob = model$z0.p)
    
    # total effect
    total <- data.frame(var = "total", 
                       estimate = model$tau.coef, 
                       ci_lower = model$tau.ci[[1]], 
                       ci_upper = model$tau.ci[[2]],
                       prob = model$tau.p)
    
    # prop. mediated
    med <- data.frame(var = "MED", 
                       estimate = model$n0, 
                       ci_lower = model$n0.ci[[1]], 
                       ci_upper = model$n0.ci[[2]],
                       prob = model$n0.p)
    
    # make a dataframe
    df <- rbind(acme, ade, total, med)
    return(df)
}

