# README.md

*updated by Katie Evans 20180824*
	- changed doseresponse HTA dilution script to save output file in directory of user choice with user defined name
*written by Katie Evans 20180510*

Scripts in this folder will generate dilutions for HTA phenotyping using the sorter for general drug assays. Different scripts are required for dose assay (located in the 'doseresponse' folder) versus a normal drug assay (located in the 'non-dose' folder).

## Dose Response

For a dose response, use template called 'HTA_dose_dilutions_setup.R' and follow steps below:
1. Edit the location of your input file 'curves'
	* curves is a .csv file in the following format: Drug, Diluent, Working, Dose_concentrations. An example can be found in the file labeled 'curves.csv'
		* Drug is the name of your compound. Can be hyphenated with '-low' or '-high' to designated two different doses of the same compound.
		* Diluent is the name of the diluent (DMSO or water usually)
		* Working is the stock concentration of the drug, in mM
		* Dose_concentrations is a vector of (usually 5) concentrations (in uM) to use in the assay, separated by commas (e.g. 0, 1, 2, 3, 4)
2. Edit the output directory and name where you want the final file to be saved (likely your own folder)
4. If necessary, you can edit the 'drug_aliquot_amount.csv' file to include your drug and the amount (in uL) that is aliquoted into each tube to figure out how many drug tubes you need to thaw. If your drug is not in the file, it will be designated NA but not return an error.
5. Run the code and knit the markdown.

The output file will be located in the directory you designated and the file and have the output name you designated (as a pdf). The first page gives you an overview of your assay: how many plates you need, how much lysate you should make, how much kanamycin to add to the lysate, etc. It also shows you all your drug dilutions you will be making as well as how many tubes of drug you will need to thaw, if applicable.

The remaining pages will each be for a plate/compound. At the top of the page it shows you again how you will dilute your drug. Then it shows how much of the diluted drug you should add to the lysate to make all your doses. Sometimes you need two different stock concentrations to keep the diluent at 1 percent, if this is the case the 'use higher stock (if necessary) will have another table showing your higher dose at the higher concentration.

**make sure your diluent + drug always is 1 percent of the total volume (lysate + diluent + drug). For doses, this is usually 12 uL diluent + drug and 1188 uL lysate (1200 uL total).

## Traditional drug HTA

For a non-dose response assay (phenotyping parents and NILs at one drug concentration, or using multiple drugs at one drug concentration each) you can use the template 'HTA_dilutions.R' and follow the steps below:
1. Edit the location of your input file 'concs'
	* concs is a .csv file in the following format: drug, diluent, startConc, finConc, numPlates, numWells
		* drug is the name of your compound
		* diluent is the name of the diluent (DMSO or water, usually)
		* startConc is the stock concentration (in mM)
		* finConc is the final concentration (in uM)
		* numPlates is the number of plates you are running with this drug
		* numWells is the number of wells on each plate you have (usually 100 is fine, unless you are running half a plate for some reason)
	* NOTE: you will want to have both your compound and diluent listed as separate drugs (e.g. to run docetaxel I need 2 lines in 'concs.csv: one for docetaxel (in DMSO) and one for DMSO (with no diluent)
	* an example can be found in the file labeled 'concs.csv'
2. Edit the output directory and name where you want the final file to be saved (likely your own folder)
3. If you want, you can choose to make a dilution of the drug (especially if you are adding ~1 uL of drug to 20 mL of lysate for example). To do this, you can change the 'dilution' variable from 1 (no dilution) to 10 (or any number, to be read as a 1:10 dilution of the drug). 
4. Another optional edit is the minimum number of uL you want to pipette for your drug (set at 0.5). Probably easier to manually set the dilution variable than mess with the minimum pipette. But either is an option.
	* Note, this script makes just the right amount of lysate needed for these assays. If you want more wiggle room, you can increase the 'leftover' variable. Default is 5.
5. Run the entire script, look for your output file in the set output directory

The top of the page gives you an overview of your assay: how many plates you need, how much lysate you should make, how much kanamycin to add to the lysate, etc. It also shows you all your drug dilutions (if you are making any) you will be making as well as how many tubes of drug you will need to thaw, if applicable. Under the 'plate dilutions' section, you will find how much of your drug and diluent you should add to your lysate to make enough for all the plates you have listed.


