#Remapping RIAILs for drugs in which we are testing ECA231 and ECA233 with the GWER threshold
#Restarted R session before running this script
setwd("~/Dropbox/AndersenLab/LabFolders/katie/projects/chemos/chrIV_NILs")

library(linkagemapping)
library(dplyr)
library(ggplot2)

RIAILs1regressed <- read.csv("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/regressedRIAILs1.csv")
RIAILs2regressed <- read.csv("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/regressedRIAILs2.csv")
allRIAILs <- rbind(RIAILs2regressed, RIAILs1regressed)

data("N2xCB4856cross")
cross <- N2xCB4856cross

#### Methotrexate
metho <- allRIAILs %>%
  dplyr::filter(condition == "methotrexate-625")

metho$condition <- "methotrexate.625"

methocross <- mergepheno(cross, metho, set = 2)

#### Docetaxel
docetaxel <- allRIAILs %>%
  dplyr::filter(condition == "docetaxel")

docetaxelcross <- mergepheno(cross, docetaxel, set = 2)

#### Tunicamycin
tunic <- allRIAILs %>%
  dplyr::filter(condition == "tunicamycin")

tuniccross <- mergepheno(cross, tunic, set = 2)

#### Deiquat
deiquat <- allRIAILs %>%
  dplyr::filter(condition == "deiquat")

deiquatcross <- mergepheno(cross, deiquat, set = 2)

#### Fluoxotene
fluoxetine <- allRIAILs %>%
  dplyr::filter(condition == "fluoxetine-250")
fluoxetine$condition <- "fluoxetine.250"

fluoxcross <- mergepheno(cross, fluoxetine, set = 2)

#### Paraquat
paraquat <- allRIAILs %>%
  dplyr::filter(condition == "paraquat")

paracross <- mergepheno(cross, paraquat, set = 2)

#### Silver
silver <- allRIAILs %>%
  dplyr::filter(condition == "silver")

silvercross <- mergepheno(cross, silver, set = 2)

#### Cisplatin
cisplat <- allRIAILs %>%
  dplyr::filter(condition == "cisplatin-250")

cisplatcross <- mergepheno(cross, cisplat, set = 2)

save(methocross, file = "methocross.Rda")
save(docetaxelcross, file = "docetaxelcross.Rda")
save(tuniccross, file = "tuniccross.Rda")
save(deiquatcross, file = "deiquatcross.Rda")
save(fluoxcross, file = "fluoxcross.Rda")
save(paracross, file = "paracross.Rda")
save(silvercross, file = "silvercross.Rda")
save(cisplatcross, file = "cisplatcross.Rda")

### Do all the mappings
methoGWER <- fsearch(methocross, permutations = 1000, threshold = "GWER")
docetaxelGWER <- fsearch(docetaxelcross, permutations = 1000, threshold = "GWER")
tunicGWER <- fsearch(tuniccross, permutations = 1000, threshold = "GWER")
deiquatGWER <- fsearch(deiquatcross, permutations = 1000, threshold = "GWER")
fluoxGWER <- fsearch(fluoxcross, permutations = 1000, threshold = "GWER")
paraGWER <- fsearch(paracross, permutations = 1000, threshold = "GWER")
silverGWER <- fsearch(silvercross, permutations = 1000, threshold = "GWER")
cisplatGWER <- fsearch(cisplatcross, permutations = 1000, threshold = "GWER")

### Save all the mappings
save(methoGWER, file = "methoGWER.Rda")
save(docetaxelGWER, file = "docetaxelGWER.Rda")
save(tunicGWER, file = "tunicGWER.Rda")
save(deiquatGWER, file = "deiquatGWER.Rda")
save(fluoxGWER, file = "fluoxGWER.Rda")
save(paraGWER, file = "paraGWER.Rda")
save(silverGWER, file = "silverGWER.Rda")
save(cisplatGWER, file = "cisplatGWER.Rda")


### Annotate all lods
methoannotatedGWER <- annotate_lods(methoGWER, methocross)
docetaxelannotatedGWER <- annotate_lods(docetaxelGWER, docetaxelcross)
tunicannotatedGWER <- annotate_lods(tunicGWER, tuniccross)
deiquatannotatedGWER <- annotate_lods(deiquatGWER, deiquatcross)
fluoxannotatedGWER <- annotate_lods(fluoxGWER, fluoxcross)
paraannotatedGWER <- annotate_lods(paraGWER, paracross)
silverannotatedGWER <- annotate_lods(silverGWER, silvercross)
cisplatannotatedGWER <- annotate_lods(cisplatGWER, cisplatcross)

### Save all annotated lods
save(methoannotatedGWER, file = "methoannotatedGWER.Rda")
save(docetaxelannotatedGWER, file = "docetaxelannotatedGWER.Rda")
save(tunicannotatedGWER, file = "tunicannotatedGWER.Rda")
save(deiquatannotatedGWER, file = "deiquatannotatedGWER.Rda")
save(fluoxannotatedGWER, file = "fluoxannotatedGWER.Rda")
save(paraannotatedGWER, file = "paraannotatedGWER.Rda")
save(silverannotatedGWER, file = "silverannotatedGWER.Rda")
save(cisplatannotatedGWER, file = "cisplatannotatedGWER.Rda")

