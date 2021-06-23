library(cegwas2)
library(tidyverse)
library(linkagemapping)
library(gsheet)

# Source nil genotype file
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_genotype_plots.R")

# load gene descriptions
load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_annotations.Rda")
data("eQTLpeaks")
data("probe_info")
load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/nilgeno.Rda")

# load primersprimers
primers <- gsheet::gsheet2tbl('https://docs.google.com/spreadsheets/d/1LJFZZ4dZm9KnoTpZqwCkQiXAKuD4Q4IGebNu1OKSnXM/edit#gid=0') %>%
    dplyr::mutate(pos = `Physical Position` / 1e6,
                  pair = paste0(`Left primer`, "-", `Right primer`),
                  pair2 = paste0(`Right primer`, "-", `Left primer`))


# Look for genes in interval
# variables
# region <- "V:5260997-5906132"

# update with "better" eQTL info?
query_genes <- function(region, GO = NULL, strain = "CB4856", v = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/WI.20200815.hard-filter.vcf.gz") {
    
    # filter eqtl to > 5% VE
    eqtlmap2 <- eQTLpeaks %>%
        dplyr::filter(var_exp >= 0.05)
    
    # how many genes are in the interval?
    all_genes <- cegwas2::query_vcf(region, impact = "ALL", samples = strain, vcf = v)
    print(glue::glue("There are {length(unique(all_genes$gene_id))} genes in the interval {region}"))
    
    # how many eQTL map to this region?
    chrom <- stringr::str_split_fixed(region, ":", 2)[,1]
    left_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,1])
    right_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,2])
    
    all_eQTL <- eqtlmap2 %>%
        dplyr::filter(chr == chrom,
                      ci_l_pos < right_pos,
                      ci_r_pos > left_pos)
    print(glue::glue("There are {nrow(all_eQTL)} eQTL ({length(unique(all_eQTL$trait))} traits) with VE > 5% that map to {region}"))
    
    # all eQTL probes
    all_eQTL_probes <- probe_info %>%
        dplyr::filter(probe %in% all_eQTL$trait)
    
    ##############################
    # if wbgene is NA - try to fix
    ##############################
    
    # filter na
    na1 <- all_eQTL_probes %>%
        dplyr::group_by(probe) %>%
        dplyr::mutate(num_na = sum(is.na(wbgene))/length(wbgene)) %>%
        dplyr::filter(num_na == 1)
    
    unique_probes <- paste(unique(na1$probe), collapse = ",")
    print(glue::glue("There are {nrow(na1)} genes with an eQTL that need to be hand curated: {unique_probes}"))
    
    ##################################
    
    # which of the eQTL are overlapping with genes in interval?
    eQTL_outside_CI <- all_eQTL_probes %>%
        dplyr::filter(!wbgene %in% all_genes$gene_id)
    print(glue::glue("There are {nrow(all_eQTL)-length(unique(eQTL_outside_CI$wbgene))-nrow(na1)} genes in the region with an eQTL, {length(unique(eQTL_outside_CI$wbgene))} genes outside the region with an eQTL, and {nrow(na1)} unknown"))
    
    # Total genes of interest:
    print(glue::glue("There are at least {length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} total genes of interest."))
    
    # how many of the genes in interval have variation?
    vars <- all_genes %>%
        dplyr::mutate(GT = ifelse(a1 == REF, "ref", "alt")) %>%
        dplyr::filter(GT == "alt")
    
    # genes with protein coding vars
    proteincode <- vars %>%
        dplyr::filter(impact %in% c("MODERATE", "HIGH"))
    print(glue::glue("There are {length(unique(vars$gene_id))}/{length(unique(all_genes$gene_id))} genes in interval with genetic variation, {length(unique(proteincode$gene_id))}/{length(unique(vars$gene_id))} have protein-coding variation"))
    
    
    # should I look at GO annotations?
    if(!is.null(GO)) {
        # total genes with GO annotations
        go_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_genes$wbgene))}/{length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} genes with {GO} annotation"))
        
        # genes with GO annotations and variation
        go_var <- gene_annotations %>%
            dplyr::filter(wbgene %in% vars$gene_id) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_var$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND genetic variation"))
        
        # genes with GO annotation and protein-coding variation
        go_pcvar <- gene_annotations %>%
            dplyr::filter(wbgene %in% proteincode$gene_id) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_pcvar$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND protein-coding genetic variation"))
        
        # genes with GO annotation and eQTL
        go_eqtl <- gene_annotations %>%
            dplyr::filter(wbgene %in% all_eQTL_probes$wbgene) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_eqtl$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND eQTL"))
        
        # return final dataframe with all info (might be off, only has 133 instead of 134?)
        total_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                          eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                          vars = ifelse(wbgene %in% vars$gene_id, T, F),
                          pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                          go_annotation = ifelse(wbgene %in% go_genes$wbgene, T, F))
        
        distinct <- total_genes %>%
            dplyr::distinct(wbgene, .keep_all = T)
        
        # pc alone
        pc <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == F))
        
        # pc + eQTL
        pc_eqtl <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == T))
        
        # eQTL alone
        e <- nrow(distinct %>% dplyr::filter(pc_vars == F, eqtl == T))
        
        print(glue::glue("There are {pc + pc_eqtl + e} genes with protein-coding variation and/or an eQTL (top priority)"))
    } else {
        
        # return final dataframe with all info (might be off, only has 133 instead of 134?)
        total_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                          eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                          vars = ifelse(wbgene %in% vars$gene_id, T, F),
                          pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                          go_annotation = NA)
        
        distinct <- total_genes %>%
            dplyr::distinct(wbgene, .keep_all = T)
        
        # pc alone
        pc <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == F))
        
        # pc + eQTL
        pc_eqtl <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == T))
        
        # eQTL alone
        e <- nrow(distinct %>% dplyr::filter(pc_vars == F, eqtl == T))
        
        print(glue::glue("There are {pc + pc_eqtl + e} genes with protein-coding variation and/or an eQTL (top priority)"))
    }
    
    return(total_genes)
}    

# test <- query_genes("III:1-500000")

# update all_probe_info with new info
# example: probe A_12_P101610 corresponds to gene WBGene00010280 and probe A_12_P111099 = WBGene00010281
# probe_ids = c("A_12_P101610", "A_12_P111099")
# gene_ids = c("WBGene00010280", "WBGene00010281")
# you can find gene ids for probes by searching the probe id on wormbase. Not all probes have a gene...

## THIS IS NOT VALID RIGHT NOW, BECAUSE I STARTED USING DATA DIRECTLY FROM LM PACKAGE
# update_probes <- function(probe_ids, gene_ids) {
#     
#     update <- all_probe_info %>%
#         dplyr::filter(probe %in% probe_ids) %>%
#         dplyr::mutate(wbgene = gene_ids,
#                       gene = NA) %>%
#         dplyr::select(probe:wbgene) %>%
#         dplyr::left_join(gene_annotations)
#     
#     new <- all_probe_info %>%
#         dplyr::filter(!probe %in% probe_ids) %>%
#         dplyr::bind_rows(update)
#     
#     # save old copy just in case of errors
#     save(all_probe_info, file = glue::glue("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/archive/{Sys.Date()}_all_probe_info.Rda"))
#     
#     
#     all_probe_info <- new
#     
#     save(all_probe_info, file = "~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/all_probe_info.Rda")
#     
#     # reload new data
#     load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/all_probe_info.Rda")
#     
#     
# }

# update sqst-5
# update_probes("A_12_P104472", "WBGene00269434")



# Script to find NILs in your area of interest

# The idea is to input a region (larger than your QTL) and the script will be able to tell you if we already have NILs here.
# INPUT: chromosome region ("IV:5000000-8000000")
# OUTPUT: list of [1] NIL genotypes on chrom of interest with region bounded.
#                 [2] NIL genotypes on all chrom
#                 [3] dataframe of NIL genotype
findNILs <- function(input) {
    
    # Get chromosome region
    chr <- stringr::str_split_fixed(input, ":", 2)[1]
    range <- stringr::str_split_fixed(input, ":", 2)[2]
    left_pos <- as.numeric(stringr::str_split_fixed(range, "-", 2)[1])
    right_pos <- as.numeric(stringr::str_split_fixed(range, "-", 2)[2])
    
    # Filter nil genotype data:
    # - filter to chromosome of interest
    # - filter to keep only ECA strains
    # - filter the length of a segment to be > 10,000 bp
    # - filter to keep only strains with N2 and CB on chrom
    nils <- nilgeno %>%
        dplyr::filter(!(sample %in% c("N2", "CB4856"))) %>%
        dplyr::mutate(length = end - start) %>%
        dplyr::group_by(sample, gt) %>%
        dplyr::mutate(background = sum(length)) %>%
        tidyr::spread(gt_name, background) %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(CB4856 = mean(CB4856, na.rm = T),
                      N2 = mean(N2, na.rm = T),
                      RIAIL = N2/CB4856) %>%
        dplyr::filter(RIAIL > 5 | RIAIL < 0.25) %>% #0.167 is one chromosome (CSS) but better safer than sorry probably.
        dplyr::group_by(sample) %>%
        dplyr::mutate(backgroundtype = ifelse(N2 > CB4856, "N2", ifelse(CB4856 > N2, "CB4856", NA))) %>%
        dplyr::filter(chrom == chr, grepl("ECA", sample)) %>%
        dplyr::filter(length > 100000) %>% # what should I use as the threshold? smallest NIL in chrIII nils is 241,230 bp
        dplyr::group_by(sample) %>%
        dplyr::mutate(num_genos = sum(plyr:::nunique(gt))) %>%
        dplyr::filter(num_genos == 2, sample != "ECA257") # this is the kammenga lab N2 
    
    # Separate into N2 and CB nils, then make sure the NIL region is within the boundary
    cbnils <- nils %>%
        dplyr::filter(backgroundtype == "N2", gt == 2, start < right_pos, end > left_pos) #1 = N2, 2 = CB
    
    n2nils <- nils %>%
        dplyr::filter(backgroundtype == "CB4856", gt == 1, start < right_pos, end > left_pos)
    
    allnils <- rbind(cbnils, n2nils)
    
    # outputs
    plot1 <- nil_plot(unique(allnils$sample), chr = chr, all.chr = F, left.bound = left_pos, right.bound = right_pos, left.cb = left_pos, left.n2 = left_pos)[[1]] +
        geom_vline(xintercept = c(left_pos/1e6, right_pos/1e6))
    
    plot2 <- nil_plot(unique(allnils$sample), chr = chr, all.chr = T, left.bound = left_pos, right.bound = right_pos, left.cb = left_pos, left.n2 = left_pos)[[1]]
    
    nildf <- nilgeno %>%
        dplyr::filter(sample %in% allnils$sample)
    
    return(list(plot1, plot2, nildf))
}

# Function to plot NIL genotype and primer locations
# primer pairs - pair of primers in the format "oECAXXX-oECAXXX", can be multiple
# NIL - ECAxxx strain name to plot genotype (N2 is default)
# RIL = false. If true, will plot RILs instead of NILs
plot_primers <- function(primer_pairs, NIL = "N2", RIL = F, chr = c("I", "II", "III", "IV", "V", "X")) {
    
    # get primer positions
    primer_positions <- primers %>%
        dplyr::filter(pair %in% primer_pairs | pair2 %in% primer_pairs)
    
    if(RIL == T) {
        df <- rilgeno
    } else {
        df <- nilgeno
    }
    
    # plot NIL genotype
    df %>%
        dplyr::filter(sample %in% NIL) %>%
        dplyr::filter(chrom %in% chr) %>%
        dplyr::mutate(gt = as.character(gt)) %>%
        ggplot2::ggplot(.)+
        ggplot2::geom_segment(ggplot2::aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
        ggplot2::facet_grid(~chrom, scales = "free",  space = "free")+
        ggplot2::scale_color_manual(values=c("1"="orange","2"="blue"))+
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=12, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=14, face="bold", color="black"),
                       axis.title.y = ggplot2::element_text(size=14, face="bold", color="black"),
                       strip.text = ggplot2::element_text(size = 12, face = "bold", color = "black"),
                       plot.title = ggplot2::element_text(size=24, face="bold"),
                       legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank())+
        ggplot2::labs(x = "Genomic Position (Mb)", y = "") +
        ggplot2::geom_vline(data = data.frame(pos = primer_positions$pos, chrom = primer_positions$Chromosome), 
                            ggplot2::aes(xintercept = pos)) +
        ggplot2::geom_text(data = data.frame(pos = primer_positions$pos, 
                                    chrom = primer_positions$Chromosome, 
                                    strain = NIL[length(NIL)],
                                    primers = primer_positions$pair), aes(x = pos,y = strain, label = primers), 
                  angle = 90,
                  hjust = -0.1,
                  vjust = -0.5)
}

# plot_primers(c("oECA799-oECA800",
#                "oECA745-oECA746"), NIL = "ECA232", chr = "V")


# function to estimate the genetic position from a physical genomic position
# input: region ("III:1000000-5000000")
# uses Ce_Genetic_Map_WS276.bed.gz downloaded and reformatted from wormbase
calc_cM <- function(region) {
    
    # parse out chromosome, pos1, and pos2
    chr <- stringr::str_split_fixed(region, ":", 2)[,1]
    positions <- stringr::str_split_fixed(region, ":", 2)[,2]
    pos1 <- as.numeric(stringr::str_split_fixed(positions, "-", 2)[,1])
    pos2 <- as.numeric(stringr::str_split_fixed(positions, "-", 2)[,2])
    
    # read gff
    gff <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/Ce_Genetic_Map_WS276.bed.gz")
    colnames(gff) <- c("chrom", "start", "stop", "cM")
    
    # convert physical to genetic 
    df <- gff %>%
        dplyr::filter(chrom == chr,
                      start >= pos1,
                      stop <= pos2)
    
    # subtract first from last cM - distance!
    glue::glue("{df$cM[1]} - {df$cM[nrow(df)]} cM ({df$cM[nrow(df)] - df$cM[1]} cM)")
}

