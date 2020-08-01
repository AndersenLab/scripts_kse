# Install and load the package
library(dplyr)
library(linkagemapping)

# Get the cross object
load_cross_obj("N2xCB4856cross_full")
cross <- N2xCB4856cross_full

# Get the phenotype data
fullpheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/zinc/data/allRIAILsregressed.rds") 

pheno <- fullpheno %>%
  dplyr::filter(condition == "methotrexate3125", !grepl("red|green|yellow|iqr|f.", trait))

methoRIAILsregressed <- pheno

# Merge the cross object and the phenotype data
methocross <- mergepheno(cross, pheno, set = 2)

# Perform a mapping with only 10 iterations of the phenotype data for FDR calc
methoGWER <- fsearch(methocross, permutations = 1000, markerset = NA, threshold = "GWER")

# Annotate the LOD scores
methoannotatedGWER <- annotate_lods(methoGWER, methocross)

methoannotatedGWER <- methoannotatedGWER %>%
  dplyr::mutate(pos = as.numeric(stringr::str_split_fixed(marker, "\\_", 2)[,2]),
                ci_l_pos = as.numeric(stringr::str_split_fixed(ci_l_marker, "\\_", 2)[,2]),
                ci_r_pos = as.numeric(stringr::str_split_fixed(ci_r_marker, "\\_", 2)[,2]))

# Save the data
save(methoRIAILsregressed, methocross, methoGWER, methoannotatedGWER, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/methotrexate_GWER_Data_full.Rdata")
save(methoannotatedGWER, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/methotrexate/data/methoannotatedGWER.RData")
save(methoRIAILsregressed, methocross, methoannotatedGWER, methoregressed, methostats, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/methotrexate/data/metho_nil_css_GWER_229.RData")
save(methoregressed, methostats, methoCSSregressed, methoCSSstats, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/methotrexate/data/metho_nil_css_229.RData")

