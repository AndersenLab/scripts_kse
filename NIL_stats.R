# Set of functions for assessing statistical significance of NILs
library(tidyverse)

# Calculate TukeyHSD statistics between pairs of conditions (such strains in NIL phenotyping assay)
# df - dataframe of phenotypes
# trt - trait of interest
# cond - condition of interest
# plot - boolean whether to plot output, default is FALSE
quick_stats <- function(df, trt, cond, plot = FALSE){
    stat_df <- df %>%
        dplyr::filter(trait == trt, condition==cond)%>%
        dplyr::select(strain, phenotype)
    
    aov_res <- aov(stat_df$phenotype ~ stat_df$strain)
    
    if(plot == TRUE) {
        summary(aov_res)
        tuk <- TukeyHSD(aov_res)
        psig=as.numeric(apply(tuk$`stat_df$strain`[,2:3],1,prod)>=0)+1
        op=par(mar=c(4.2,9,3.8,2))
        plot(tuk,col=psig,yaxt="n")
        for (j in 1:length(psig)){
            axis(2,at=j,labels=rownames(tuk$`stat_df$strain`)[length(psig)-j+1],
                 las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
        }
        par(op)
    }
    
    pwtuk <- TukeyHSD(aov_res)
    
    return(pwtuk)
}

# Wrapper function for `quick_stats` to do statitistics for each strain pair-wise comparison for each condition-trait sample.
# dfregressed - dataframe of phenotypes
# pval - what pvalue is considered significant? Defaults to 0.05
get_stats <- function(dfregressed, pval = 0.05) {
    statsdf <- dfregressed %>%
        ungroup() %>%
        dplyr::select(condition, trait) %>%
        dplyr::distinct(condition, trait, .keep_all = T)
    
    #Add statistical significance to each pair of strains
    for(i in 1:nrow(statsdf)) {
        stats <- broom::tidy(quick_stats(df = dfregressed, trt = statsdf$trait[i], cond = statsdf$condition[i], plot = FALSE))
        for(j in 1:nrow(stats)){
            if(stats$adj.p.value[j] < pval) {
                comp <- stats$comparison[j]
                if(!(comp %in% names(statsdf))) {
                    statsdf[comp] <- c(NA)
                }
                statsdf[i, comp] <- TRUE
            }
        }
    }
    return(statsdf)
}
