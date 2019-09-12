library(tidyverse)

load("~/Dropbox/AndersenLab/Reagents/WormReagents/Variation/N2_GFF_WS245/gene_ref_flat.Rda")

gene_model <- function(df, genename = NA, WBID = NA, gene_color = "purple", intron_color = "black",
                       utr3_color = "gray60", utr5_color = "gray60", gene_alpha = 0.5) {
    if(is.na(genename)){
        if(is.na(WBID)) {
            stop("No gene provided. Give either gene name or WB ID")
        } else {
            selection <- df %>%
                dplyr::filter(wbgene == WBID)
        }
    } else {
        selection <- df %>%
            dplyr::filter(gene == genename)
    }
    
    if(nrow(selection) > 1){
        stop("multiple entries with given gene name.")
    }
    
    if(selection$strand == "+"){
        
        selection$txend <- as.numeric(selection$txend)
        selection$codingend <- as.numeric(selection$codingend)
        exonstarts <- matrix(as.numeric(unlist(strsplit(selection$exonstarts, ",")), ncol = selection$numexons))
        
        exonends <- matrix(as.numeric(unlist(strsplit(selection$exonends, ",")), ncol = selection$numexons))
        
        allexons <- data.frame("starts" = exonstarts, "ends" = exonends)
        
        intronstarts <- exonends[1:(nrow(exonends)-1)]
        
        intronends <- exonstarts[2:(nrow(exonends))]
        
        allintrons <- data.frame("starts" = intronstarts, "ends" = intronends) %>%
            dplyr::mutate(midpoint = (starts+ends)/2)
        
        endUTR <- data.frame("x" = c(selection$codingend, selection$codingend, selection$txend), "y" = c(1, -1, 0))
        
        plot <- ggplot(allexons)+
            geom_segment(aes(x = txstart, xend = txend, y = 0, yend = 0), data = selection, color = "black")+
            geom_rect( aes(xmin =  starts, xmax = ends, ymin = -1 , ymax = 1), fill = gene_color, color = "black", alpha = gene_alpha)+
            geom_segment(aes(x = starts, y = 1, xend = midpoint, yend = 2), data = allintrons, color = intron_color)+
            geom_segment(aes(x = midpoint, y = 2, xend = ends, yend = 1), data = allintrons, color = intron_color)+
            geom_rect(aes(xmin = txstart, xmax = codingstart, ymin = -1, ymax = 1), data = selection, fill = utr3_color)+
            geom_rect(aes(xmin = codingend, xmax = txend, ymin = -1, ymax = 1), data = selection, fill= "white", color = "white", lwd = 1.2)+
            geom_polygon(aes(x = x, y = y), data = endUTR, fill = utr5_color, color = "black")+
            theme_void()
        
        return(plot)
        
    } else {
        
        selection$txend <- as.numeric(selection$txend)
        selection$codingend <- as.numeric(selection$codingend)
        exonstarts <- matrix(as.numeric(unlist(strsplit(selection$exonstarts, ",")), ncol = selection$numexons))
        
        exonends <- matrix(as.numeric(unlist(strsplit(selection$exonends, ",")), ncol = selection$numexons))
        
        allexons <- data.frame("starts" = exonstarts, "ends" = exonends)
        
        intronstarts <- exonends[1:(nrow(exonends)-1)]
        
        intronends <- exonstarts[2:(nrow(exonends))]
        
        allintrons <- data.frame("starts" = intronstarts, "ends" = intronends) %>%
            dplyr::mutate(midpoint = (starts+ends)/2)
        
        endUTR <- data.frame("x" = c(selection$codingstart, selection$codingstart, selection$txstart), "y" = c(1, -1, 0))
        
        plot <- ggplot(allexons)+
            geom_segment(aes(x = txend, xend = txstart, y = 0, yend = 0), data = selection, color = "black")+
            geom_rect( aes(xmin =  starts, xmax = ends, ymin = -1 , ymax = 1), fill = gene_color, color = "black", alpha = gene_alpha)+
            geom_segment(aes(x = starts, y = 1, xend = midpoint, yend = 2), data = allintrons, color = intron_color)+
            geom_segment(aes(x = midpoint, y = 2, xend = ends, yend = 1), data = allintrons, color = intron_color)+
            geom_rect(aes(xmin = txend, xmax = codingend, ymin = -1, ymax = 1), data = selection, fill = utr3_color)+
            geom_rect(aes(xmin = codingstart, xmax = txstart, ymin = -1, ymax = 1), data = selection, fill= "white", color = "white", lwd = 1.2)+
            geom_polygon(aes(x = x, y = y), data = endUTR, fill = utr5_color, color = "black")+
            theme_void()
        
        return(plot)
    }
}

# # gene model
# gene_model(gene_ref_flat, WBID = "WBGene00010471")
# 
# # deletions
# gene_model(gene_ref_flat, WBID = "WBGene00010471")+
#     geom_segment(aes(x = 12413036, xend = 12413576, y = -1.25, yend = -1.25), color = "red", lwd =3) +
#     geom_segment(aes(x = 12413023, xend = 12413753, y = -1.5, yend = -1.5), color = "blue", lwd =3) +
#     geom_segment(aes(x = 12413047, xend = 12413719, y = -1.75, yend = -1.75), color = "green", lwd =3) +
#     geom_segment(aes(x = 12412981, xend = 12413792, y = -2, yend = -2), color = "orange", lwd =3)
# 
# ggsave("~/Dropbox/AndersenLab/LabFolders/Katie/projects/chemos/zinc/figures/20190513_N2_cdr6_deletions.pdf", height = 4, width = 10)

