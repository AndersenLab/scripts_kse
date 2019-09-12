library(easysorter)
library(dplyr)
library(broom)

# Define a vector of your experiement directories
dirs <- "~/Dropbox/HTA/Results/20180410_abachrV/"

# Read in the data
raw <- read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam <- remove_contamination(raw)

# Summarize the data
# directories = FALSE because supplied data comes from one directory
summedraw <- sumplate(raw_nocontam, directories = FALSE, quantiles = TRUE)

#Prune based on biological impossibilities
biopruned <- bioprune(summedraw) %>%
    tidyr::gather(trait, phenotype, -(date:col)) %>%
    dplyr::filter(!grepl("red|green|yellow|f.|iqr", trait))

# Prune based on bins
# don't use bamfprune for assays with high replication
#bamfpruned <- bamf_prune(biopruned, drop = TRUE)

# prune based on outliers 2 standard deviations away from the mean
pruned <- NULL
for(drug in unique(biopruned$condition)) {
    for(i in unique(biopruned$trait)) {
        test <- biopruned %>%
            dplyr::filter(trait == i, condition == drug) %>%
            dplyr::group_by(strain) %>%
            dplyr::mutate(med = median(phenotype), med1 = med + 2*IQR(phenotype),
                          mph = mean(phenotype), sph = sd(phenotype)) %>%
            dplyr::filter(phenotype <= 2*sph+mph,
                          phenotype >= mph-2*sph)
        pruned <- rbind(pruned, test)
    }
}

regressed <- regress(pruned)
