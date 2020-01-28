library(tidyverse)
library(ggplot2)
library(cowplot)

# quickly plot NIL phenotype. N2 will be colored orange and CB blue (NILs are grey)
# df - dataframe of NIL phenotype (from sorter)
# cond - condition to plot
# pltrt - trait to plot
# ylab - y label, default is condition.trait
# textsize - size of text in plot, default is 12
# titlesize - size of titles in plot, defalt is 14
# pointsize - size of points in plot, default is 0.5
quick_plot <- function(df, cond, pltrt, ylab = paste0(cond, ".", pltrt),
                       textsize = 12, titlesize = 14, pointsize = 0.5) {
    phen_gen <- df %>%
        dplyr::filter(trait == pltrt, condition == cond) %>%
        dplyr::mutate(type = ifelse(as.character(strain) == "N2", "N2_parent", ifelse(as.character(strain) == "CB4856", "CB_parent", "NIL")))
    cbdf <- phen_gen[phen_gen$strain=="CB4856",]
    n2df <- phen_gen[phen_gen$strain=="N2",]
    phen_gen %>%
        ggplot(.) +
        aes(x = factor(strain),
            y = phenotype, 
            fill=factor(type)) +
        geom_jitter(size = pointsize, width = 0.1)+
        geom_boxplot(outlier.colour = NA, alpha = 0.7)+
        scale_fill_manual(values = c("N2_parent" = "orange", "CB_parent" = "blue", "NIL" = "gray"))+
        theme_bw()+
        theme(axis.text.x = element_text(size=textsize, face="bold", color="black", angle = 90),
              axis.text.y = element_text(size=textsize, face="bold", color="black"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=titlesize, face="bold", color="black"),
              strip.text.x = element_text(size=textsize, face="bold", color="black"),
              strip.text.y = element_text(size=textsize, face="bold", color="black"),
              plot.title = element_text(size=titlesize, face="bold", vjust = 1),
              legend.position="none",
              panel.background = element_rect( color="black",size=1.2),
              strip.background = element_rect(color = "black", size = 1.2),
              panel.border = element_rect( colour = "black"))+
        labs(title = ylab)
}

# quickly plot NIL phenotype. N2 will be colored orange and CB blue (NILs are grey). Orientation will be flipped to show nil genotypes in combination with phenotypes (doesn't plot genotypes though)
# df - dataframe of NIL phenotype (from sorter)
# cond - condition to plot
# pltrt - trait to plot
# ylab - y label, default is condition.trait
# textsize - size of text in plot, default is 12
# titlesize - size of titles in plot, defalt is 14
# pointsize - size of points in plot, default is 0.5
quick_plot_breakup_flip <- function(df, cond, pltrt, ylab = paste0(cond, ".", pltrt),
                                    textsize = 12, titlesize = 14, pointsize = 0.5) {
    phen_gen <- df %>%
        dplyr::filter(trait == pltrt, condition == cond) %>%
        dplyr::mutate(type = ifelse(as.character(strain) == "N2", "N2_parent", ifelse(as.character(strain) == "CB4856", "CB_parent", "NIL")))
    cbdf <- phen_gen[phen_gen$strain=="CB4856",]
    n2df <- phen_gen[phen_gen$strain=="N2",]
    phen_gen %>%
        ggplot(.) +
        aes(x = factor(strain),
            y = phenotype, 
            fill=factor(type)) +
        geom_jitter(size = pointsize, width = 0.1)+
        geom_boxplot(outlier.colour = NA, alpha = 0.7)+
        # scale_fill_manual(values = c("N2_parent" = "orange", "CB_parent" = "blue", "NIL" = "gray"))+
        scale_fill_manual(values = c("N2_parent" = "#F0BA51", "CB_parent" = "#484DA0", "NIL" = "gray"))+
        theme_bw()+
        coord_flip()+
        theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.title.x = element_text(size=titlesize, face="bold", color="black", vjust=-.3),
              axis.title.y = element_blank(),
              strip.text.x = element_text(size=textsize, face="bold", color="black"),
              strip.text.y = element_text(size=textsize, face="bold", color="black"),
              plot.title = element_text(size=titlesize, face="bold", vjust = 1),
              legend.position="none",
              panel.background = element_rect( color="black",size=1.2),
              strip.background = element_rect(color = "black", size = 1.2),
              panel.border = element_rect( colour = "black"))+
        labs(y = ylab)
}

# plot all QTL from linkagemapping for several traits/conditions
# annotatedmap - annotated mapping (result from `linkagemapping::annotate_lods()`)
# nils - buggy, might not work. Supply a dataframe of nil genotype information as ci_l_pos and ci_r_pos define the region of the NIL. Will be plotted as a red rectangle on the plot.
all_lod_plots <- function(annotatedmap, nils = NULL) {
    newmap <- annotatedmap %>%
        arrange(chr, ci_l_pos, ci_r_pos) %>%
        na.omit()
    
    newmap <- newmap %>%
        dplyr::mutate(condition = stringr::str_split_fixed(.$trait, "\\.", 2)[,1],
                      trait = stringr::str_split_fixed(.$trait, "\\.", 2)[,2],
                      n2res = ifelse(eff_size < 0, "yes", "no"))
    
    faketrait <- newmap$trait[1]
    fakecond <- newmap$condition[1]
    #Set chromosome boundaries
    newrows <- newmap[1,] 
    newrows[1,] = c(NA,"I",5000000,faketrait,0,NA,NA,NA,NA,NA,1,NA,14972282, fakecond, "yes")
    newrows[2,] = c(NA,"II",5000000,faketrait,0,NA,NA,NA,NA,NA,1,NA,15173999, fakecond, "yes")
    newrows[3,] = c(NA,"III",5000000,faketrait,0,NA,NA,NA,NA,NA,1,NA,13829314, fakecond, "yes")
    newrows[4,] = c(NA,"IV",5000000,faketrait,0,NA,NA,NA,NA,NA,1,NA,17450860, fakecond, "yes")
    newrows[5,] = c(NA,"V",5000000,faketrait,0,NA,NA,NA,NA,NA,1,NA,20914693, fakecond, "yes")
    newrows[6,] = c(NA,"X",5000000,faketrait,0,NA,NA,NA,NA,NA,1,NA,17748731, fakecond, "yes")
    newrows$ci_l_pos <- as.numeric(newrows$ci_l_pos)
    newrows$ci_r_pos <- as.numeric(newrows$ci_r_pos)
    newrows$pos <- as.numeric(newrows$pos)
    newrows$lod <- as.numeric(newrows$lod)
    
    if(is.null(nils)) {
        #Plot
        ggplot(newmap)+
            aes(x=pos/1E6, y=trait)+
            theme_bw() +
            viridis::scale_fill_viridis(name = "LOD") + 
            viridis::scale_color_viridis(name = "LOD") +
            geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
            geom_segment(data=newrows,aes(x = 0, y = trait, xend = ci_r_pos/1e6, yend = trait), size = 2.5, alpha = 0) +
            geom_point(aes(fill = lod, shape = n2res), color = "black",size = 3, alpha = 1)+
            scale_shape_manual(values = c("yes" = 24, "no" = 25)) +
            xlab("Genomic position (Mb)") + ylab("") +
            guides(shape = FALSE) +
            theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                  axis.ticks.y = element_blank(),
                  legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
                  legend.key.size = unit(.75, "cm"),
                  panel.grid.major.x = element_line(),
                  panel.grid.major.y = element_line(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.y = element_text(size = 10, face = "bold", color = "black"),
                  axis.title.x = element_text(size=12, face="bold", color= "black"),
                  axis.title.y = element_blank(),
                  strip.text.x = element_text(size=12, face="bold", color="black"),
                  strip.text.y = element_text(size=12, face="bold", color="black", angle = 0),
                  strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                  plot.title = element_text(size=12, face="bold")) +
            facet_grid(condition ~ chr, scales = "free_x", space = "free")
    } else {
        #Plot
        ggplot(newmap)+
            aes(x=pos/1E6, y=trait)+
            theme_bw() +
            viridis::scale_fill_viridis(name = "LOD") + viridis::scale_color_viridis(name = "LOD") +
            geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
            geom_segment(data=newrows,aes(x = 0, y = trait, xend = ci_r_pos/1e6, yend = trait), size = 2.5, alpha = 0) +
            geom_rect(data=nils, aes(xmin = ci_l_pos/1e6, ymin = "cv.EXT", xmax = ci_r_pos/1e6, ymax = "var.TOF"), size = 2, alpha = 0.2, fill = "red")+
            geom_point(aes(fill = lod, shape = n2res), color = "black",size = 3, alpha = 1)+
            scale_shape_manual(values = c("yes" = 24, "no" = 25)) +
            xlab("Genomic position (Mb)") + ylab("") +
            guides(shape = F) +
            theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                  axis.ticks.y = element_blank(),
                  legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
                  legend.key.size = unit(.75, "cm"),
                  panel.grid.major.x = element_line(),
                  panel.grid.major.y = element_line(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.y = element_text(size = 10, face = "bold", color = "black"),
                  axis.title.x = element_text(size=12, face="bold", color= "black"),
                  axis.title.y = element_blank(),
                  strip.text.x = element_text(size=12, face="bold", color="black"),
                  strip.text.y = element_text(size=12, face="bold", color="black", angle = 0),
                  strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                  plot.title = element_text(size=12, face="bold")) +
            facet_grid(condition ~ chr, scales = "free_x", space = "free")
        
        # need to fix the NIL segment to be dynamic based on how many traits there are
    }
    
}

# plot phenotype x genotype splits for RIAILs
# cross - cross object containing genotype and phenotype of RIAILs (output of `linkagemapping::mergepheno()`)
# map - annotated mapping (result from `linkagemapping::annotate_lods()`)
# parent - usually "N2xCB4856"
# tit - title for plot. Default is None.
# ylab - add y label for plot, default is None.
# textsize - size of text in plot, default is 8
# titlesize - size of titles in plot, defalt is 16
# pointsize - size of points in plot, default is 0.5
pxgplot_kt <- function (cross, map, parent = "N2xCB4856", tit = "", ylab = "",
                        textsize = 8, titlesize = 16, pointsize = 0.5) {
    peaks <- map %>% dplyr::group_by(iteration) %>% dplyr::filter(!is.na(var_exp)) %>% 
        dplyr::do(head(., n = 1))
    if (nrow(peaks) == 0) {
        stop("No QTL identified")
    }
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    pheno <- cross$pheno %>% 
        dplyr::select_(map$trait[1])
    geno <- data.frame(linkagemapping:::extract_genotype(cross)) %>% 
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
        data.frame(., pheno)
    colnames(geno)[1:(ncol(geno) - 1)] <- sapply(colnames(geno)[1:(ncol(geno) - 1)], 
                                                 function(marker) {paste(unlist(peaks[peaks$marker == gsub("\\.", "-", marker),
                                                                                      c("chr", "pos")]), collapse = ":")})
    colnames(geno)[ncol(geno)] <- "pheno"
    split <- tidyr::gather(geno, marker, genotype, -pheno)
    split$genotype <- sapply(split$genotype, function(x) {
        if (is.na(x)) {
            return(NA)
        }
        if (parent == "N2xCB4856") {
            if (x == -1) {
                "N2"
            }
            else {
                "CB4856"
            }
        }
        else if (parent == "N2xLSJ2") {
            if (x == -1) {
                "N2"
            }
            else {
                "LSJ2"
            }
        }
        else if (parent == "AF16xHK104") {
            if (x == -1) {
                "AF16"
            }
            else {
                "HK104"
            }
        }
    })
    # split$genotype <- factor(split$genotype, levels = c("N2", "CB4856", "LSJ2", "AF16", "HK104"),
    #                          labels = c("N2-RIAILs", "CB-RIAILs", "LSJ2-RIAILs", "AF16-RIAILs", "HK104-RIAILs"))
    split$genotype <- factor(split$genotype, levels = c("N2", "CB4856", "LSJ2", "AF16", "HK104"),
                             labels = c("N2 Allele", "CB4856 Allele", "LSJ2-RIAILs", "AF16-RIAILs", "HK104-RIAILs"))
    split <- split %>%
        tidyr::drop_na(genotype) %>%
        dplyr::mutate(chr = (as.character(stringr::str_split_fixed(marker, ":", 2)[,1])),
                      pos = as.numeric(as.character(stringr::str_split_fixed(marker, ":", 2)[,2]))) %>%
        dplyr::arrange(chr, pos)
    split$marker <- factor(split$marker, levels = unique(split$marker))
    ggplot2::ggplot(split) + 
        ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = 1, size = pointsize, width = 0.1) + 
        ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype, alpha = 0.8), outlier.shape = NA) + 
        ggplot2::scale_fill_manual(values = c(`N2 Allele` = "orange", `CB4856 Allele` = "blue", LSJ2 = "green", AF16 = "indianred", HK104 = "gold")) + 
        ggplot2::facet_wrap(~marker, ncol = 5) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = textsize,face = "bold", color = "black"), 
                       axis.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"), 
                       axis.title.x = ggplot2::element_text(size = titlesize, face = "bold", color = "black", vjust = -0.3), 
                       axis.title.y = ggplot2::element_text(size = titlesize, face = "bold", color = "black"), 
                       strip.text.x = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
                       strip.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
                       plot.title = ggplot2::element_blank(), 
                       legend.position = "none", 
                       panel.background = ggplot2::element_rect(color = "black",size = 1.2)) + 
        ggplot2::labs(x = "", y = ylab)
}




# plot phenotype x genotype splits for RIAILs, including N2/CB parents!
# cross - cross object containing genotype and phenotype of RIAILs (output of `linkagemapping::mergepheno()`)
# map - annotated mapping (result from `linkagemapping::annotate_lods()`)
# parpheno - dataframe for parent phenotype (this is not found in the crossobject so must be supplemented)
# tit - title for plot. Default is None.
# ylab - add y label for plot, default is None.
# textsize - size of text in plot, default is 8
# titlesize - size of titles in plot, defalt is 16
# pointsize - size of points in plot, default is 0.5
pxgplot_par_kt <- function (cross, map, parpheno, tit = "", ylab = "",
                        textsize = 8, titlesize = 16, pointsize = 0.5) {
    peaks <- map %>% 
        dplyr::group_by(iteration) %>% 
        dplyr::filter(!is.na(var_exp)) %>% 
        dplyr::do(head(., n = 1))
    
    if (nrow(peaks) == 0) {
        stop("No QTL identified")
    }
    
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    
    pheno <- cross$pheno %>% 
        dplyr::select_(map$trait[1])
    
    geno <- data.frame(linkagemapping:::extract_genotype(cross)) %>% 
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
        data.frame(., pheno)
    
    colnames(geno)[1:(ncol(geno) - 1)] <- sapply(colnames(geno)[1:(ncol(geno) - 1)], 
                                                 function(marker) {paste(unlist(peaks[peaks$marker == gsub("\\.", "-", marker),
                                                                                      c("chr", "pos")]), collapse = ":")})
    colnames(geno)[ncol(geno)] <- "pheno"
    
    split <- tidyr::gather(geno, marker, genotype, -pheno) %>%
        tidyr::drop_na(genotype)
    
    split$genotype <- sapply(split$genotype, function(x) {
        if (x == -1) {
            "N2-RIAILs"
        }
        else {
            "CB-RIAILs"
        }
    })
    
    # add parent phenotype
    parpheno <- parpheno %>% 
        dplyr::mutate(marker = " Parental") %>%
        dplyr::mutate(genotype = strain) %>%
        dplyr::select(pheno = phenotype, marker, genotype)
    
    split <- split %>%
        dplyr::bind_rows(parpheno)
    
    split$genotype <- factor(split$genotype, levels = c("N2", "CB4856", "N2-RIAILs", "CB-RIAILs"))
    
    split <- split %>%
        tidyr::drop_na(genotype) %>%
        dplyr::mutate(chr = (as.character(stringr::str_split_fixed(marker, ":", 2)[,1])),
                      pos = as.numeric(as.character(stringr::str_split_fixed(marker, ":", 2)[,2]))) %>%
        dplyr::arrange(chr, pos)
    split$marker <- factor(split$marker, levels = unique(split$marker))

    ggplot2::ggplot(split) + 
        ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = 1, size = pointsize, width = 0.1) + 
        ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype, alpha = 0.8), outlier.shape = NA) + 
        ggplot2::scale_fill_manual(values = c(`N2-RIAILs` = "orange", `CB-RIAILs` = "blue", "N2" = "orange", "CB4856" = "blue")) + 
        ggplot2::facet_wrap(~marker, ncol = 5, scales = "free_x") + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = textsize,face = "bold", color = "black"), 
                       axis.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"), 
                       axis.title.x = ggplot2::element_text(size = titlesize, face = "bold", color = "black", vjust = -0.3), 
                       axis.title.y = ggplot2::element_text(size = titlesize, face = "bold", color = "black"), 
                       strip.text.x = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
                       strip.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
                       plot.title = ggplot2::element_blank(), 
                       legend.position = "none", 
                       panel.background = ggplot2::element_rect(color = "black",size = 1.2)) + 
        ggplot2::labs(x = "", y = ylab)
}

# plot linkagemapping results
# map - annotated mapping (result from `linkagemapping::annotate_lods()`)
# textsize - size of text in plot, default is 12
# titlesize - size of titles in plot, defalt is 16
# linesize - changes the size of the line of the linkagemap (important mostly for generating large figures for posters), default is 1
# col - color of confidence interval, default is blue
maxlodplot_kt <- function (map, textsize = 12, titlesize = 16, linesize = 1, col = "blue") {
    map1 <- map %>% dplyr::group_by(marker) %>% dplyr::filter(lod ==max(lod))
    cis <- map %>% dplyr::group_by(marker) %>% dplyr::mutate(maxlod = max(lod)) %>%
        dplyr::group_by(iteration) %>% dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n = 1))
    if (nrow(cis) == 0) {
        plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype,y = pheno)) + ggplot2::geom_blank()
        return(plot)
    }
    map1 <- linkagemapping:::cidefiner(cis, map1)
    plot <- ggplot2::ggplot(map1) + 
        ggplot2::aes(x = pos/1e+06,y = lod)
    if (nrow(cis) != 0) {
        plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = col, alpha = 0.5) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                                size = 3.2, show.legend = FALSE) + 
            ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2 * maxlod), label = paste0(100 *round(var_exp, digits = 4), "%")), 
                               colour = "black",size = textsize / 4, hjust = "inward")
    }
    plot <- plot + ggplot2::geom_line(size = linesize, alpha = 0.85) +
        ggplot2::facet_grid(. ~ chr, scales = "free", space = "free") +
        ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
        ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
        ggplot2::ggtitle(map1$trait[1]) + 
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = textsize, color = "black", face = "bold"),
            axis.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
            axis.title.x = ggplot2::element_text(size = titlesize, face="bold", color="black"),
            axis.title.y = ggplot2::element_text(size = titlesize, face="bold", color="black"),
            strip.text.x = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
            strip.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
            plot.title = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(color = "black", size = 1.2))
    return(plot)
}

maxlodplot_kt_noVE <- function (map, textsize = 12, titlesize = 16, linesize = 1, col = "blue") {
    map1 <- map %>% dplyr::group_by(marker) %>% dplyr::filter(lod ==max(lod))
    cis <- map %>% dplyr::group_by(marker) %>% dplyr::mutate(maxlod = max(lod)) %>%
        dplyr::group_by(iteration) %>% dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n = 1))
    if (nrow(cis) == 0) {
        plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype,y = pheno)) + ggplot2::geom_blank()
        return(plot)
    }
    map1 <- linkagemapping:::cidefiner(cis, map1)
    plot <- ggplot2::ggplot(map1) + 
        ggplot2::aes(x = pos/1e+06,y = lod)
    if (nrow(cis) != 0) {
        plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = col, alpha = 0.5) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                                size = textsize/4, show.legend = FALSE)
            # ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2 * maxlod), label = paste0(100 *round(var_exp, digits = 4), "%")), 
                               # colour = "black",size = textsize / 4, hjust = "inward")
    }
    plot <- plot + ggplot2::geom_line(size = linesize, alpha = 0.85) +
        ggplot2::facet_grid(. ~ chr, scales = "free", space = "free") +
        ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
        ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
        ggplot2::ggtitle(map1$trait[1]) + 
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = textsize, color = "black", face = "bold"),
            axis.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
            axis.title.x = ggplot2::element_text(size = titlesize, face="bold", color="black"),
            axis.title.y = ggplot2::element_text(size = titlesize, face="bold", color="black"),
            strip.text.x = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
            strip.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
            plot.title = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(color = "black", size = 1.2))
    return(plot)
}
