library(tidyverse)
library(ggplot2)

#load nilgenos
load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/nilgeno.Rda")

# load RIAIL genotypes
load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/rilgeno.Rda")


# df - nil genotype df in segment format
# chr - NIL chromosome
# left.cb - if looking at NIL breakups, this number corresponds to the left flank for CB4856 NILs
# left.n2 - if looking at NIL breakups, this number corresponds to the left flank for N2 NILs
# left.bound - left boundary for region of interest
# right.bound - right boundary for region of interest
# scan.range - cutoff for small genomic regions when identifying left nils
# all.chr - plot all chromosomes or just the one with the NIL
# section - "all" returns all NILs, "N2-NILs" returns only NILs N2 > CB. "CB-NILs" returns only NILs CB > N2
# background - FALSE returns the genotype of just the chromosome or all chromosomes. TRUE returns just the chrom of interest and the "genome" genotype
# ci - default is NA (no lines drawn) otherwise input vector of positions for confidence intervals of QTL - plots only on chromosome of interest (chr)
# order - should the function guess the best order for your NILs or should you define this? Default is TRUE (function decides). If FALSE, the order of your NILs
#         is determined by the order you write them in "strains"
# elements - should the function return each individual plot and dataset (T) or should it just plot it all (F)? For ease, choose F,
#         if you want to manipulate the plots later, maybe choose T
# n2_color and cb_color - can change to a different "orange" and "blue" if preferred. Pretty option: 
# tsize - theme size. Default is 12

nil_plot <- function(strains, chr, left.cb = 0, left.n2 = 0, left.bound = 0, right.bound = 19e6, scan.range = 2e4, 
                     all.chr=F, section = "all", background = F, ci = NA, elements = F, order = T, n2_color = "orange", cb_color = "blue",
                     tsize = 12){
    
    # copy nilgeno
    nilgeno2 <- nilgeno
    
    # check to see if all strains are represented in the nilgeno2 file
    df <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::distinct(sample)
    
    if(nrow(df) != length(strains)) {
        unknown_strains <- setdiff(strains, df$sample)
       
        # make fake genotype for unknown strains
        for(s in unknown_strains) {
            fake_geno <- nilgeno2 %>%
                dplyr::filter(sample == "N2")%>%
                dplyr::distinct(chrom, .keep_all = T) %>%
                dplyr::mutate(sample = s,
                              gt = 3,
                              gt_name = "unknown")
            
            # add fake geno to nil geno
            nilgeno2 <- rbind(nilgeno2, fake_geno)
        }
    }
    
    # # # determine if NILs are CB or N2
    nilsII_sort_type <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::group_by(sample)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::group_by(sample, gt)%>%
        dplyr::mutate(major_gt = sum(gt_ct))%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(major_gt))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::mutate(nil_type = ifelse(gt == 1, "CB", "N2"))%>%
        dplyr::select(sample, nil_type)

    # # # keep NILs that lost NIL genotype on right side
    nilsII_left <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::filter((start > left.cb - .3e6 & start < left.cb + .3e6 & nil_type == "CB") |
                          (start > left.n2 - .3e6 & start < left.n2 + .3e6 & nil_type == "N2") )%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::mutate(side = "LEFT")%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(size))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::filter(size > scan.range)%>% # # # remove small (likely wrong calls) around interval site
        dplyr::select(sample, side)

    # # # keep NILs that lost NIL genotype on left side
    nilsII_right <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::filter(!sample %in% nilsII_left$sample)%>%
        dplyr::mutate(side = "RIGHT")%>%
        dplyr::distinct(sample,.keep_all = T)%>%
        dplyr::select(sample, side)

    nil_sides <- bind_rows(nilsII_left,nilsII_right)


    nilsII_sort_left <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::left_join(.,nil_sides, by = "sample")%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::group_by(sample)%>%
        dplyr::filter(start > left.bound & start < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::filter(side == "LEFT")%>%
        dplyr::arrange( desc(gt_ct))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::arrange(nil_type, desc(gt_ct))

    #Here
    nilsII_sort_right <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::left_join(.,nil_sides, by = "sample")%>%
        dplyr::filter(side == "RIGHT")%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::group_by(sample)%>%
        dplyr::filter(start > left.bound & start < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::arrange( desc(gt_ct))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::ungroup()%>%
        dplyr::arrange( nil_type, gt_ct)

    N2CB <- nilgeno2 %>%
        dplyr::filter(sample %in% c("N2", "CB4856"), chrom == chr) %>%
        dplyr::mutate(nil_type = "parent", side = NA, size = NA, gt_ct = NA)

    nilsII_sort <- bind_rows(nilsII_sort_right, nilsII_sort_left) %>%
        arrange(desc(nil_type), desc(side), size)

    nilsII_sort <- rbind(nilsII_sort, N2CB)
    
    if(nrow(df) != length(strains)) {
        nilsII_sort <- rbind(nilsII_sort, nilgeno2 %>%
                             dplyr::filter(sample %in% unknown_strains,
                                           chrom == chr) %>%
                                 dplyr::mutate(nil_type = "unkown", side = NA, size = NA, gt_ct = NA))

    }
    
    if (all.chr == T){

        nilsII <- nilgeno2 %>%
            dplyr::filter(sample %in% strains) %>%
            dplyr::filter(chrom != "MtDNA")%>%
            dplyr::left_join(.,nilsII_sort_type, by = "sample")
        
        if(nrow(df) != length(strains)) {
            nilsII <- nilsII %>%
                dplyr::mutate(nil_type = ifelse(gt_name == "unknown", "unknown", nil_type))
        }

        if(order == T) {
            nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
        } else {
            nilsII$sample <- factor(nilsII$sample, levels = unique(strains), ordered = T)
        }
        nilsII$gt <- as.character(nilsII$gt)
    } else {
        nilsII <- nilgeno2 %>%
            dplyr::filter(sample %in% strains) %>%
            dplyr::filter(chrom == chr, chrom != "MtDNA")%>%
            dplyr::left_join(.,nilsII_sort_type, by = "sample")
        
        if(nrow(df) != length(strains)) {
            nilsII <- nilsII %>%
                dplyr::mutate(nil_type = ifelse(gt_name == "unknown", "unknown", nil_type))
        }

        if(order == T) {
            nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
        } else {
            nilsII$sample <- factor(nilsII$sample, levels = unique(strains), ordered = T)
        }
        nilsII$gt <- as.character(nilsII$gt)
    }
    
    # make fake 'background' genome
    if(background == T) {
        # cannot show all chromosomes and the background
        if(all.chr == T) {
            bgplot <- NULL
        } else {
            bg <- data.frame(chrom = NA, start = NA, end = NA, sample = NA, gt = NA, gt_name = NA, supporting_sites = NA, sites = NA,
                             DP = NA, switches = NA, nil_type = NA)
            for(i in unique(nilsII$sample)) {
                # get the NIL type
                type <- nilsII %>%
                    dplyr::filter(sample == i) %>%
                    dplyr::distinct(nil_type) %>%
                    dplyr::pull(nil_type)
                
                # add a new row with the background genotype in a "new" chromosome
                bg <- rbind(bg, c("Genome", 1, 3e6, i, 
                                  ifelse(type == "N2", 2,
                                         ifelse(type == "CB", 1,
                                                ifelse(type == "unknown", 3, NA))),
                                  ifelse(type == "N2", "CB4856",
                                         ifelse(type == "CB", "N2",
                                                ifelse(type == "unknown", "unknown", NA))),
                                  NA, NA, NA, NA, type))
            }
            bg <- bg %>%
                tidyr::drop_na(sample) %>%
                dplyr::mutate(start = as.numeric(start), end = as.numeric(end))
            
            # add "chrom" to chromosome for plotting
            nils2 <- nilsII %>%
                dplyr::mutate(chrom = paste0("chr", chrom))
            
            # background plot
            bgplot <- ggplot(bg)+
                geom_segment(aes(x = start/1e6, y = factor(sample, levels = levels(nils2$sample)), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
                scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
                # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
                facet_grid(~chrom, scales = "free",  space = "free")+
                theme_bw(tsize) +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      strip.text = element_text(face = "bold", color = "black"),
                      plot.title = element_text(face="bold"),
                      legend.position = "none",
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      axis.ticks = element_blank())
        }
    } else { bgplot <- NULL }
    
    # only plot the N2-NILs
    if(section == "N2-NILs") {
        nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
            theme_bw(tsize) +
            theme(axis.text.x = element_text(face="bold", color="black"),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.title.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  plot.title = element_text(face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(data = filter(nilsII, chrom == chr) %>%
                           dplyr::distinct(sample, chrom) %>%
                           dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                           tidyr::separate_rows(xint), 
                       aes(xintercept = as.numeric(xint)), color = "red") +
            labs(x = "Genomic position (Mb)", y = "")
    } else if(section == "CB-NILs") {
        # only plot the CB-NILs
        nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
            theme_bw(tsize) +
            theme(axis.text.x = element_text(face="bold", color="black"),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.title.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  plot.title = element_text(face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(data = filter(nilsII, chrom == chr) %>%
                           dplyr::distinct(sample, chrom) %>%
                           dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                           tidyr::separate_rows(xint), 
                       aes(xintercept = as.numeric(xint)), color = "red") +
            labs(x = "Genomic position (Mb)", y = "")
    } else {
        
        # plot all NILs
        nl.pl <- ggplot(nilsII)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
            theme_bw(tsize) +
            theme(axis.text.x = element_text(face="bold", color="black"),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.title.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  plot.title = element_text(face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(data = filter(nilsII, chrom == chr) %>%
                           dplyr::distinct(sample, chrom) %>%
                           dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                           tidyr::separate_rows(xint), 
                       aes(xintercept = as.numeric(xint)), color = "red") +
            labs(x = "Genomic position (Mb)", y = "")
    }
    
    # return plots and dataframe
    if(is.null(bgplot)) {
        return(list(nl.pl, nilsII_sort, nilsII))
    } else {
        if(elements == T) {
            return(list(nl.pl, bgplot, nilsII_sort))
        } else {
            return(list(cowplot::plot_grid(nl.pl, bgplot, nrow = 1, ncol = 2, align = "h", axis = "b", rel_widths = c(1, 0.3)),
                        nilsII_sort,
                        nilsII))
            # return(list(nl.pl,
            #             nilsII_sort,
            #             nilsII,
            #             bgplot))
        }
    }
}


# plots genotype on left and phenotype on right for dataframe, condition, trait
# dataframe needs condition and trait columns
plot_genopheno <- function(pheno, cond, trt, chrom, back = F, conf = 1) {
    # source phenotype plot code
    source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_phenotype_plots.R")
    
    # plot genotype
    nilgeno_plot <- nil_plot(unique(pheno$strain), chr = chrom, background = back, ci = conf, left.bound = 0)
    
    # order phenotype in same level as genotype plot
    pheno$strain <- factor(pheno$strain, levels = unique(nilgeno_plot[[2]]$sample), 
                           labels = unique(nilgeno_plot[[2]]$sample), ordered = T)
    
    # Plot NIL phenotypes in condition
    plot <- cowplot::plot_grid(nilgeno_plot[[1]], 
                               quick_plot_breakup_flip(pheno, cond, trt) + facet_grid(~trait), nrow = 1, ncol = 2)
    
    # return plot
    return(plot)
}

# Function to plot RIL genotypes based on WGS data from 20170327 - 
# includes 1200 strains, all of set 1 and set 2 and "E" strains (not set3?)
# default plots all chromosomes, you can provide vector of chromosomes to plot
# must supply either a vector of strains or a strain set (1 or 2)
# can add strain labels with strain_label = T
plot_rilgeno <- function(strains = NA, chr = c("I", "II", "III", "IV", "V", "X"), strainset = NA, theme_size = 12, strain_label = F) {
    # if set != NA, plot all RILs in the set
    if(!is.na(strainset)) {
        if(strainset %in% c(1, 2)) {
            df <- rilgeno %>%
                dplyr::filter(chrom %in% chr, set == strainset)
        } else {
            stop("Error. Must supply either vector of strains or strain set (1 or 2)")
        }
    } else {
        if(length(strains) == 1 && is.na(strains)) {
            stop("Error. Must supply either vector of strains or strain set (1 or 2)")
        } else {
            df <- rilgeno %>%
                dplyr::filter(chrom %in% chr, sample %in% strains)
        }
    }
    
    # add or remove strain labels, depending on flag
    if(strain_label == T) {
        df %>%
            ggplot(.)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
            theme_bw(theme_size) +
            theme(axis.text = element_text(color="black"),
                  axis.title = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  legend.position = "none",
                  panel.grid = element_blank()) +
            labs(x = "Genomic Position (Mb)", y = "")
    } else {
        df %>%
            ggplot(.)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
            theme_bw(theme_size) +
            theme(axis.text.x = element_text(color="black"),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  legend.position = "none",
                  panel.grid = element_blank()) +
            labs(x = "Genomic Position (Mb)", y = "")
    }
    
}


