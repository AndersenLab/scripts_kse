# scripts_kse

## NIL_genotype_plots.R
To plot genotypes of N2/CB NILs using WGS data from the AndersenLab.

### nil_plot(strains, chr, ...)
*options*
- strains - vector of strain names 
- chr - NIL chromosome
- left.cb - if looking at NIL breakups, this number corresponds to the left flank for CB4856 NILs
- left.n2 - if looking at NIL breakups, this number corresponds to the left flank for N2 NILs
- left.bound - left boundary for region of interest
- right.bound - right boundary for region of interest
- scan.range - cutoff for small genomic regions when identifying left nils
- all.chr - plot all chromosomes or just the one with the NIL
- section - "all" returns all NILs, "N2-NILs" returns only NILs N2 > CB. "CB-NILs" returns only NILs CB > N2
- background - FALSE returns the genotype of just the chromosome or all chromosomes. TRUE returns just the chrom of interest and the "genome" genotype
- ci - default is NA (no lines drawn) otherwise input vector of positions for confidence intervals of QTL

### plot_genopheno(pheno, cond, trt, chrom)
Use with `NIL_phenotype_plots.R` to create side-by-side plots of NIL genotypes and phenotypes
- pheno - dataframe of phenotypes
- cond - drug condition
- trt - trait
- chrom - NIL chromosome
- back - show the NIL background (genotype of parental background)? Default FALSE
- conf- default is NA (no lines drawn) otherwise input vector of positions for confidence intervals of QTL

## NIL_phenotype_plots.R
To plot a variety of phenotypes plots for NILs and RIAILs.

### pxgplot_kt(cross, map)

### pxgplot_par_kt(cross, map, parpheno)

### maxlodplot_kt(map)

### quick_plot_breakup_flip(df, cond, pltrt)

### quick_plot(df, cond, pltrt)

### all_lod_plots(annotatedmap)

## NIL_stats.R

### quick_stats(df, trt, cond)

### get_stats(dfregressed)
Wrapper function for `quick_stats` to do statitistics for each strain pair-wise comparison for each condition-trait sample.

## pca_analysis.R

### calc_pc2(pheno)
Returns a list of (1) pca object and (2) pca phenotypes

### calc_pc(pheno)
*how is this different from calc_pc2?*

### predict_pc(pheno, pca_obj)
Predicts pca phenotypes given a pca object (use for applying pca from linkage to NILs)

## RIL_genotype_plots.R

### plot_rilgeno(strains, chr, strainset) 
Function to plot RIL genotypes based on WGS data from 20170327 - includes 1200 strains, all of set 1 and set 2 and "E" strains (not set3?). default plots all chromosomes, you can provide vector of chromosomes to plot. must supply either a vector of strains or a strain set (1 or 2). can add strain labels with strain_label = T

## qtl_narrow.R

### qtl_narrow(query)
Wrapper function for `cegwas::query_vcf()` to only apply for N2/CB variants and add GO terms and gene descriptions
*are these the most up to date go and descriptions???*

## plot_primers.R

### plot_primers(primer_pairs)
Script to plot NIL genotype and primer locations
- primer pairs - pair of primers in the format "oECAXXX-oECAXXX", can be multiple (vector)
- NIL - ECAxxx strain name to plot genotype (N2 is default)
- RIL = false. If true, will plot RILs instead of NILs
- chr = default show all chr but you can just choose one (or multiple)

## gene_model_plot.R

### gene_model(genename, WBID)
Provide genename or Wormbase ID to plot gene model (exons/introns)
- can provide color with `gene_color = "purple"` and `intron_color = "black"` and `utr3_color = "gray60"` and `utr5_color = "gray60"`

## findNILs.R

### findNILs(input)
Script to find NILs in your area of interest. Input a region larger than your QTL and the script will be able to tell you if we already have NILs here.
- INPUT: chromosome region ("IV:5000000-8000000")
- OUTPUT: list of [1] NIL genotypes on chrom of interest with region bounded.
				  [2] NIL genotypes on all chrom
				  [3] dataframe of NIL genotype

### eQTL_mediate_dQTL.R

## eQTL_mediate_dQTL(peak, probe, phenodf)

### burdentest.R

### easysorter_script.R

### GWERmapping.R

### heritability.R

### HTA_sorter_template.Rmd & HTA_sorter_template2.Rmd

### linkagemapping.R

### physical2genetic_map.R
Script can be used to find the genetic position of QTL in RIAILs. 
Can also be used to find the expected recombination frequency between N2/CB

## p2gmap(chr, pos)
function to estimate genetic position from a physical genomic positon 
- chr = chromosome
- pos = basepair position

## genetic_distance(region)
estimate centimorgan distance between two physical positions
- input: region in format of `chr:pos1-pos2`

### powercalc.R

### RIAIL_to_NIL.R

### rialsplit.R