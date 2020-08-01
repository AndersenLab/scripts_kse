# Setup script for HTA_dose_dilutions

##############
#You need to edit these lines to reflect your directory setup and file name
##############

#Get the compound names, molecular weights, solvents, and dilution curves from a csv file entitled "DrugDoseResponses.csv"
curves <- read.csv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/HTA_sorter/20190305_mtx_zincV/day1_input copy.csv", stringsAsFactors = F)

# EDIT name/directory for output file:
outputdir <- "~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/HTA_sorter/20190305_mtx_zincV/"
outputname <- "zinc_v4"

# If your drug is not listed here you can add it to this file:
aliquotAmt <- read.csv("~/Dropbox/AndersenLab/RCode/HTA_dilutions/drug_aliquot_amount.csv", strip.white = TRUE)

options(knitr.duplicate.label = "allow")
#Render markdown
rmarkdown::render("~/Dropbox/AndersenLab/RCode/HTA_dilutions/doseresponse/HTA_dose_dilutions.Rmd", output_file = paste0(outputdir, outputname, ".pdf")) 
