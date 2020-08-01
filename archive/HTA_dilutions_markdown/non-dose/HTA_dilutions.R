library(plyr)
library(dplyr)
library(knitr)
library(DT)
library(rmarkdown)
library(kableExtra)

##---------------------------##
##    EDIT THE FOLLOWING:    ##
##---------------------------##

# Edit the input file:
concs <- read.csv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/zinc/data/nils/chrV/20190929_zincV/input.csv")


# EDIT name/directory for output file:
outputdir <- "~/Dropbox/AndersenLab/LabFolders/Katie/projects/zinc/data/nils/chrV/20190929_zincV/"
outputname <- "zincV_dilutions"

# VERSION 2 or VERSION 3 HTA?
version <- 3

##---------------------------##
##     OPTIONAL EDITS        ##
##---------------------------##

# Will require a dilution by this number of each drug prior to adding to the lysate (10 recommended to reduce pipetting error for doses)
# If you don't want to dilute drug, use 1
dilution <- 10

# If your drug is not listed here you can add it to this file:
aliquotAmt <- read.csv("~/Dropbox/AndersenLab/RCode/HTA_dilutions/drug_aliquot_amount.csv", strip.white = TRUE)

# Always 1%
dilConc <- 0.01

# Minimum amount to pipette: 0.5 uL
min <- 0.5

# Can adjust if you want more room for error, I find 5 works well usually
leftover <- 5

##---------------------------##
##                           ##
##---------------------------##

#How many plates?
concs <- concs %>%
    dplyr::mutate(totalWells = numPlates*numWells)
plates <- paste0("Need ", floor(sum(concs$totalWells/96)), " plates.")

# set volume and concentration based on version
if(version == 3) {
    vol <- 25
    concs <- concs %>%
        dplyr::mutate(intConc = startConc / dilution * 3)
} else if(version == 2) {
    vol <- 50
    concs <- concs %>%
        dplyr::mutate(intConc = startConc / dilution)
} else {
    stop("Error. Must choose version 2 or version 3")
}

#Drug Plates
drugconc <- concs %>%
  dplyr::mutate(numWells = numWells + leftover, #Use 100 wells/plate instead of 96 wells/plate to allow for error
                total_mix = numPlates * numWells * vol,
                drugAmt = (total_mix * finConc) / (intConc * 1000)) %>%
  dplyr::group_by(drug) %>%
  dplyr::mutate(totaldrug = sum(drugAmt),
                drugDilute = ceiling((totaldrug * intConc) / startConc),
                dilDilute = (drugDilute * dilution) - drugDilute,
                dilAmt = (total_mix * dilConc) - drugAmt,
                lysate = total_mix - (drugAmt + dilAmt)) %>%
  dplyr::left_join(aliquotAmt, by = "drug") %>%
  dplyr::mutate(drugTubes = ifelse((drugDilute %% aliquotAmt) == 0, ceiling(drugDilute / aliquotAmt) + 1, 
                                    ceiling(drugDilute / aliquotAmt)))

# Check to make sure dilution is enough, should not pipette less than the minimum amount
tooSmall <- drugconc %>%
  dplyr::filter(drugAmt < min, drugAmt > 0)
drugconc2 <- drugconc

while(nrow(tooSmall) > 0) {
  drugconc2 <- drugconc2 %>%
    dplyr::select(names(concs))
  oldconcs <- drugconc2 %>%
    dplyr::filter(!(drug %in% unique(tooSmall$drug)))
  concs2 <- drugconc2 %>%
    dplyr::filter(drug %in% unique(tooSmall$drug)) %>%
    dplyr::mutate(intConc = intConc / 10)
  newconcs <- rbind(oldconcs, concs2)
  
  drugconc2 <- newconcs %>%
    dplyr::mutate(total_mix = numPlates * numWells * vol,
                  drugAmt = (total_mix * finConc) / (intConc * 1000)) %>%
    dplyr::group_by(drug) %>%
    dplyr::mutate(totaldrug = sum(drugAmt),
                  drugDilute = ceiling((totaldrug * intConc) / startConc),
                  dilDilute = (drugDilute * startConc/intConc) - drugDilute,
                  dilAmt = (total_mix * dilConc) - drugAmt,
                  lysate = total_mix - (drugAmt + dilAmt)) %>%
    dplyr::left_join(aliquotAmt, by = "drug") %>%
    dplyr::mutate(drugTubes = ifelse((drugDilute %% aliquotAmt) == 0, ceiling(drugDilute / aliquotAmt) + 1, 
                                     ceiling(drugDilute / aliquotAmt)))
  
  tooSmall <- drugconc2 %>%
    dplyr::filter(drugAmt < min, drugAmt > 0)
}

concs <- drugconc2

#Special conditions
for(i in 1:length(concs$drug)){
  if(grepl("DMSO|water", concs$drug[i])) {
    concs$dilAmt[i] <- concs$total_mix[i]*dilConc
    concs$drugAmt[i] <- 0
    concs$totaldrug[i] <- 0
    concs$drugDilute[i] <- 0
    concs$dilDilute[i] <- 0
    concs$aliquotAmt[i] <- 0
    concs$lysate[i] <-  concs$total_mix[i] - concs$dilAmt[i]
    concs$drugTubes[i] <- 1
  } else
    if(concs$drug[i] == "cisplatin") {
      concs$dilAmt[i] <- 0
    }
}

#Truncate decimals
concs$dilAmt <- round(concs$dilAmt, digits = 2)
concs$drugAmt <- round(concs$drugAmt, digits = 2)
concs <- concs %>%
  dplyr::select(drugTubes, drug, diluent, drugDilute, dilDilute, drugAmt, dilAmt, lysate, total_mix, startConc, intConc, finConc, numPlates)

# change concentration of lysate depending on version 2 or 3
if(version == 2) {
    # make 10 mg/mL for 96 hours
    #If lysate is exact, add another tube
    lysate <- sum(concs$total_mix)
    if((lysate %% 10000) == 0) {
        lysate <- lysate + 10000
    }
    
    lysate_mix <- (paste0("Total lysate mix needed: ", lysate / 1000, " mL (10 mg/mL), actually make: ", ceiling(lysate / 10000)*10, " mL"))
    lysate_tubes <- (paste0("Number of lysate tubes needed: ", ceiling(lysate / 10000)))
    lysate <- ceiling(lysate / 10000) *10000
    lysate_K <- (paste0("Add ",  lysate / 10000 * 9, " mL of K medium to lysate."))
    kan <- (paste0("Add ", 50*lysate/80000, " uL of Kanamycin to lysate mix."))
} else {
    # make 15 mg/mL which will be diluted 1:3 in the well
    #If lysate is exact, add another tube
    lysate <- sum(concs$total_mix)
    if((lysate %% 6666) == 0) {
        lysate <- lysate + 10000
    }
    
    lysate_mix <- (paste0("Total lysate mix needed: ", lysate / 1000, " mL (15 mg/mL), actually make: ", ceiling(lysate / 6666)*6.666, " mL"))
    lysate_tubes <- (paste0("Number of lysate tubes needed: ", ceiling(lysate / 6666)))
    lysate <- ceiling(lysate / 6666) *6666
    lysate_K <- (paste0("Add ",  lysate/1000 - ceiling(lysate / 6666), " mL of K medium to lysate."))
    kan <- (paste0("Add ", 3*50*lysate/80000, " uL of Kanamycin to lysate mix."))
}

drug_dilution <- concs %>%
  dplyr::mutate(drugDilute = ifelse(dilDilute == 0, NA, drugDilute), dilDilute = ifelse(dilDilute == 0, NA, dilDilute)) %>%
  dplyr::select(tubes = drugTubes, drug, diluent, `drugDilute (uL)` = drugDilute, `dilDilute (uL)` = dilDilute, `startConc (mM)` = startConc, 
                `intConc (mM)` = intConc) %>%
  dplyr::filter(!(drug %in% c("DMSO", "water"))) %>%
  dplyr::distinct() 


plate_dilution <- concs %>%
  dplyr::mutate(lysate = lysate / 1000,
                total_mix = total_mix / 1000) %>%
  dplyr::select(drug, diluent, `drug (uL)` = drugAmt, `dil (uL)` = dilAmt, `lysate (mL)` = lysate, `total (mL)` = total_mix, 
                `startConc (mM)` = startConc, `finConc (uM)` = finConc, plates = numPlates)

special <- NULL
if("cisplatin" %in% unique(concs$drug)){
  special <- ("For cisplatin: heat at 55°C with shaking for 5 minutes and make dilution right after.")
}


#Render markdown
rmarkdown::render("~/Dropbox/AndersenLab/RCode/HTA_dilutions/non-dose/sorter_dilution.Rmd", output_file = paste0(outputdir, outputname, ".pdf")) 

