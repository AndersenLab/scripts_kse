# KSE Scripts for Andersen Lab

## NIL_genotype_plots.R
*To plot genotypes of N2/CB NILs using WGS data from the Andersen Lab.*

### nil_plot(strains, chr, ...)

**Inputs**
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

**Output**  
Output is a list of plots/dataframes.  
[[1]]: NIL genotype plot  
[[2]]: Simplififed dataframe of NIL genotypes to plot, ordered as they are plotted  
[[3]]: Complete NIL genotype dataframe at chromosome of interest or all chromosomes  

**Example**  
`nil_plot(strains = c("N2", "CB4856", "ECA411", "ECA481"), chr = "V", ci = 10, background = T)[[1]]`


### plot_genopheno(pheno, cond, trt, chrom)
Use with `NIL_phenotype_plots.R` to create side-by-side plots of NIL genotypes and phenotypes

**Inputs**
- pheno - dataframe of phenotypes
- cond - drug condition
- trt - trait
- chrom - NIL chromosome
- back - show the NIL background (genotype of parental background)? Default FALSE
- conf- default is NA (no lines drawn) otherwise input vector of positions for confidence intervals of QTL

**Output**  
Returns a plot of NIL phenotype/genotype.

**Example**  
`plot_genopheno(pheno = regressed_df, cond = "zinc", trt = "mean.EXT", chrom = "V", back = T, conf = 10)`

### plot_rilgeno(strains, chr, strainset) 
Function to plot RIL genotypes based on WGS data from 20170327 - includes 1200 strains, all of set 1 and set 2 and "E" strains (not set3?). Default plots all chromosomes, you can provide vector of chromosomes to plot. Must supply *either* a vector of strains *or* a strain set (1 or 2).

**Inputs**
- strains - vector of strains to plot
- chr - vector of chromosomes to include, default is all
- strainset - can use ***instead of*** `strains` to plot all srains in set 1 or set 2 (not set 3)
- theme_size - size of plot theme (text etc)
- strain_label - boolean to include strain names on plot, default to FALSE

**Output**
Returns plot of RIAIL genotypes

**Example**
`plot_rilgeno(chr = "III", strainset = 2, strainl_label = FALSE)`

---

## NIL_phenotype_plots.R
*To plot a variety of phenotypes plots for NILs and RIAILs.*

### quick_plot(df, cond, pltrt)
Quickly plot NIL phenotype. N2 will be colored orange and CB blue (NILs are grey)

**Inputs**
- df - dataframe of NIL phenotype (from sorter)
- cond - condition to plot
- pltrt - trait to plot
- ylab - y label, default is condition.trait
- textsize - size of text in plot, default is 12
- titlesize - size of titles in plot, defalt is 14
- pointsize - size of points in plot, default is 0.5

**Output**  
Returns plot of NIL phenotypes (strain on x axis and phenotype on y axis)

**Example**  
`quick_plot(df = regressed_df, cond = "zinc", pltrt = "mean.TOF")`

### quick_plot_breakup_flip(df, cond, pltrt)
Quickly plot NIL phenotype. N2 will be colored orange and CB blue (NILs are grey). Orientation will be flipped to show nil genotypes in combination with phenotypes (doesn't plot genotypes though)

**Inputs**
- df - dataframe of NIL phenotype (from sorter)
- cond - condition to plot
- pltrt - trait to plot
- ylab - y label, default is condition.trait
- textsize - size of text in plot, default is 12
- titlesize - size of titles in plot, defalt is 14
- pointsize - size of points in plot, default is 0.5

**Output**  
Returns plot of NIL phenotypes (strain on y axis and phenotype on x axis)

**Example**  
`quick_plot_breakup_flip(df = regressed_df, cond = "zinc", pltrt = "mean.TOF")`

### all_lod_plots(annotatedmap)
Plot all QTL from linkagemapping for several traits/conditions (like Figure 1 from QTL hotspot paper Evans and Brady et al., 2018)

**Inputs**
- annotatedmap - annotated mapping (result from `linkagemapping::annotate_lods()`)
- nils - buggy, might not work. Supply a dataframe of nil genotype information as `ci_l_pos` and `ci_r_pos` define the region of the NIL. Will be plotted as a red rectangle on the plot.

**Output**  
Plot with positions of multiple QTL as a dot (peak marker) and line (confidence interval)

**Example**  
`all_lod_plots(annotatedmap = zincGWER)`

### pxgplot_kt(cross, map)
Plot phenotype x genotype splits for RIAILs

**Inputs**
- cross - cross object containing genotype and phenotype of RIAILs (output of `linkagemapping::mergepheno()`)
- map - annotated mapping (result from `linkagemapping::annotate_lods()`)
- parent - usually "N2xCB4856"
- tit - title for plot. Default is None.
- ylab - add y label for plot, default is None.
- textsize - size of text in plot, default is 8
- titlesize - size of titles in plot, defalt is 16
- pointsize - size of points in plot, default is 0.5

**Output**  
Plot phenotype (y axis) by genotype (x axis) for RIAILs

**Example**  
`pxgplot_kt(cross = zinccross, map = zincGWER_meanTOF)`

### pxgplot_par_kt(cross, map, parpheno)
Plot phenotype x genotype splits for RIAILs, including N2/CB parents!

**Inputs**
- cross - cross object containing genotype and phenotype of RIAILs (output of `linkagemapping::mergepheno()`)
- map - annotated mapping (result from `linkagemapping::annotate_lods()`)
- parpheno - dataframe for parent phenotype (this is not found in the crossobject so must be supplemented)
- tit - title for plot. Default is None.
- ylab - add y label for plot, default is None.
- textsize - size of text in plot, default is 8
- titlesize - size of titles in plot, defalt is 16
- pointsize - size of points in plot, default is 0.5

**Output**  
Same as `pxgplot_kt` but includes facet of parents AND RIAILs

**Example**  
`pxgplot_par_kt(cross = zinccross, map = zincGWER_meanTOF, parpheno = zincparents)`

### maxlodplot_kt(map)
Plot linkagemapping results

**Inputs**
- map - annotated mapping (result from `linkagemapping::annotate_lods()`)
- textsize - size of text in plot, default is 12
- titlesize - size of titles in plot, defalt is 16
- linesize - changes the size of the line of the linkagemap (important mostly for generating large figures for posters), default is 1
- col - color of confidence interval, default is blue

**Output**  
LOD plot for linkage mapping trait

**Example**  
`maxlodplot_kt(zincGER_meanTOF, col = "pink")`

---

## NIL_stats.R
*Set of functions for assessing statistical significance of NILs*

### quick_stats(df, trt, cond)
Calculate TukeyHSD statistics between pairs of conditions (such strains in NIL phenotyping assay)

**Inputs**
- df - dataframe of phenotypes
- trt - trait of interest
- cond - condition of interest
- plot - boolean whether to plot output, default is FALSE

**Output**  
Returns matrix of statistical significance for all pairwise comparisons in dataframe

**Example**  
`quick_stats(df = regressed_df, trt = "mean.TOF", cond = "zinc")`

### get_stats(dfregressed)
Wrapper function for `quick_stats` to do statitistics for each strain pair-wise comparison for each condition-trait sample.

**Inputs**
- dfregressed - dataframe of phenotypes
- pval - what pvalue is considered significant? Defaults to 0.05

**Output**  
Returns dataframe of statistical significance for all pairwise comparisons, tidied with `broom::tidy` and TRUE/FALSE for significance, based on the given p-value

**Example**  
`get_stats(dfregressed = regressed_df, pval = 0.01)`

---

## pca_analysis.R
*Set of functions to calculate principal components of phenotypes*

### calc_pc_reps(pheno)
Calculate principal components from phenotype dataframe with replicates (i.e. NILs)

**Inputs**
- pheno - dataframe of phenotypes

**Output**  
Returns a list of:  
[[1]] PCA object (including loadings)  
[[2]] PCA phenotypes (only including PCs that account for 90% of variation in phenotype)

**Example**  
`calc_pc_reps(pheno = nil_pheno_regressed)`

### calc_pc_noreps(pheno)
Calculate principal components from phenotype dataframe with no replicates (i.e. RIAILs)

**Inputs**
- pheno - dataframe of phenotypes

**Output**  
Returns a list of:  
[[1]] PCA object (including loadings)  
[[2]] PCA phenotypes (only including PCs that account for 90% of variation in phenotype)

**Example**  
`calc_pc_reps(pheno = riail_pheno_regressed)`

### predict_pc(pheno, pca_obj)
Predicts PCA phenotypes given loadings from another pca object (use for applying pca from linkage to NILs)

**Inputs**
- pheno - dataframe of phenotypes
- pca_oject - pca object, output [[1]] from calc_pc_reps or calc_pc_noreps
- keep - how many PCs to keep? default is number of traits in pheno

**Output**  
Returns dataframe of PC phenotypes

**Example**  
`predict_pc(pheno = nil_pheno_regressed, pca_obj = riail_pca_object, keep = 6)`

---

## NIL_narrowing.R

### qtl_narrow(query)
Wrapper function for `cegwas::query_vcf()` to show variants between N2/CB within a QTL interval and add gene functions/GO terms

**Inputs**
- query - region of interest (III:1000-700000)
- sev - vector of severity to include (MODIFIER, LOW, MODERATE, HIGH), default is all

**Output**  
Dataframe with all variants in region, annotated with gene description, GO term, and whether CB is alt/ref

**Example**  
`qtl_narrow("III:1000-700000", sev = c("MODERATE", "HIGH"))`

### findNILs(input)
Script to find NILs in your area of interest. Input a region larger than your QTL and the script will be able to tell you if we already have NILs here.

**Inputs**
- input - chromosome region ("IV:5000000-8000000")

**Output**  
Returns a list of  
[1] NIL genotypes on chrom of interest with region bounded.  
[2] NIL genotypes on all chrom  
[3] dataframe of NIL genotype  

**Example**  
`findNILs(input = "IV:5000000-8000000")`

### plot_primers(primer_pairs)
Function to plot NIL genotype and primer locations

**Inputs**
- primer_pairs - pair of primers in the format "oECAXXX-oECAXXX", can be multiple (vector)
- NIL - ECAxxx strain name to plot genotype (N2 is default)
- RIL = false. If true, will plot RILs instead of NILs
- chr = default show all chr but you can just choose one (or multiple)

**Output**  
Returns plot with locations of primers and genotypes of NIL/RIAILs

**Example**  
`plot_primers(primer_pairs = c("oECA1148-oECA1147", "oECA1408-oECA1409", "oECA1141-oECA1142"), NIL = "ECA1065")`

### genetic_distance(region)
Estimate centimorgan distance between two physical positions

**Inputs**
- region - region in format of `chr:pos1-pos2`

**Output**
Returns

### p2gmap(chr, pos)
function to estimate genetic position from a physical genomic positon 
- chr = chromosome
- pos = basepair position


---



## gene_model_plot.R

### gene_model(genename, WBID)
Provide genename or Wormbase ID to plot gene model (exons/introns)
- can provide color with `gene_color = "purple"` and `intron_color = "black"` and `utr3_color = "gray60"` and `utr5_color = "gray60"`

---

## eQTL_mediate_dQTL.R

### eQTL_mediate_dQTL(peak, probe, phenodf)

---

## burdentest.R

---

## easysorter_script.R

---

## GWERmapping.R

---

## heritability.R

---

## HTA_sorter_template.Rmd & HTA_sorter_template2.Rmd

---

## linkagemapping.R


---

## powercalc.R

---

## RIAIL_to_NIL.R

---

## rialsplit.R

---