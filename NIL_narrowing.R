library(cegwas2)
library(tidyverse)
library(linkagemapping)
library(gsheet)

# Source nil genotype file
source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_genotype_plots.R")

# load gene descriptions
load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/gene_descriptions_WS273.Rda")

# load primersprimers
primers <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1LJFZZ4dZm9KnoTpZqwCkQiXAKuD4Q4IGebNu1OKSnXM/edit#gid=0') %>%
    dplyr::mutate(pos = `Physical Position` / 1e6,
                  pair = paste0(`Left primer`, "-", `Right primer`),
                  pair2 = paste0(`Right primer`, "-", `Left primer`))

# Define function to show variants between N2/CB within a QTL interval and add gene functions/GO terms
# query - region of interest (III:1000-700000)
# sev - vector of severity to include (MODIFIER, LOW, MODERATE, HIGH), default is all
qtl_narrow <- function(query, sev = c("MODIFIER", "LOW", "MODERATE", "HIGH")) {
    
    # snpeff <- cegwas::snpeff(query)
    snpeff <- cegwas2::query_vcf(query, impact = sev, samples = "CB4856")
    print(paste0("Total genes in interval: ", length(unique(snpeff$gene_id))))
    
    #Only keep genes that have difference in N2/CB
    snpf <- snpeff  %>%
        dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF, "REF", "ALT"))
    
    #Add GO Terms and descriptions
    df <- snpf %>%
        dplyr::select(CHROM, POS, strain = SAMPLE, REF, ALT, a1, a2, GT, query, allele, effect, impact, gene_name, gene_id, feature_type, 
                      feature_id, transcript_biotype, nt_change, aa_change) %>%
        dplyr::left_join(gene_descriptions, by = c("gene_id" = "wbgene"))
    
    return(df)
}

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
        ggplot(.)+
        geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
        facet_grid(~chrom, scales = "free",  space = "free")+
        scale_color_manual(values=c("1"="orange","2"="blue"))+
        theme_bw() +
        theme(axis.text.x = element_text(size=12, face="bold", color="black"),
              axis.text.y = element_text(size=12, face="bold", color="black"),
              axis.title.x = element_text(size=14, face="bold", color="black"),
              axis.title.y = element_text(size=14, face="bold", color="black"),
              strip.text = element_text(size = 12, face = "bold", color = "black"),
              plot.title = element_text(size=24, face="bold"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        labs(x = "Genomic Position (Mb)", y = "") +
        geom_vline(data = data.frame(pos = primer_positions$pos, chrom = primer_positions$Chromosome), aes(xintercept = pos)) +
        geom_text(data = data.frame(pos = primer_positions$pos, 
                                    chrom = primer_positions$Chromosome, 
                                    strain = NIL[length(NIL)],
                                    primers = primer_positions$pair), aes(x = pos,y = strain, label = primers), 
                  angle = 90,
                  hjust = -0.1,
                  vjust = -0.5)
}

# function to estimate the genetic position from a physical genomic position
# input: chr = chromosome ("I", "II", "III", "IV", "V", or "X")
# input: pos = basepair position
p2gmap <- function(chr, pos) {
    # read gff
    gff <- readr::read_tsv("~/Downloads/c_elegans.current.genetic_limits.gff2", col_names = F, col_types = cols())
    
    # add column names
    colnames(gff) <- c("chrom", "pmap", "gmap_span", "pmap1", "pmap2", "test1", "test2", "test3", "info")
    
    # tidy data
    test <- gff %>%
        dplyr::select(chrom, pmap1, pmap2, info) %>%
        tidyr::separate(info, c("gene", "gmap", "note"), sep = ";") %>%
        dplyr::mutate(gene = stringr::str_split_fixed(gene, "GMap ", 2)[,2]) %>%
        dplyr::mutate(gmap = stringr::str_split_fixed(gmap, '"', 3)[,2]) %>%
        dplyr::select(-note) %>%
        dplyr::mutate(gmap = stringr::str_split_fixed(gmap, " cM", 2)[,1]) %>%
        dplyr::mutate(gmap = as.numeric(gmap))
    
    # find genetic position with physical position closest to marker
    test2 <- test %>%
        dplyr::filter(chrom == chr, pmap1 >= pos - 100000, pmap1 <= pos + 100000) %>%
        tidyr::gather(gene2, pos2, pmap1:pmap2) %>%
        dplyr::mutate(dist = abs(pos2 - pos)) %>%
        dplyr::arrange(dist)
    
    return(test2$gmap[1])
}

# function to estimate the centimorgan distance between two physical positions
# input: region in format of `chr:pos1-pos2`
genetic_distance <- function(region) {
    # parse out chromosome, pos1, and pos2
    chr <- stringr::str_split_fixed(region, ":", 2)[,1]
    positions <- stringr::str_split_fixed(region, ":", 2)[,2]
    pos1 <- as.numeric(stringr::str_split_fixed(positions, "-", 2)[,1])
    pos2 <- as.numeric(stringr::str_split_fixed(positions, "-", 2)[,2])
    
    return(glue::glue("{abs(p2gmap(chr, pos1) - p2gmap(chr, pos2))} cM"))
    
}

