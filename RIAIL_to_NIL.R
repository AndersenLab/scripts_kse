#### Script takes genotype data for RIAILs and Identifies RIAILs that can be used for generating NILs

#### Input
#CI.L - left marker name for confidence interval
#CI.R - right marker name for confidence interval
#chromosome - chromosome that interval is on
#dStrain - desired genotype you would like to generate NIL block with ("N2", or "CB4856")
#dir - directory where genotype data is located
#     NOTE** on formatting of genotype data
#           column structure - 
#                    id - marker name
#                    V2 - chromosome number
#                    V3 - centimorgan distance
#                    followed by one column for each strain
#                      Genotypes should be in "AA" or "AB" format
#nstrains - number of desired strains you want output to have breaks

#### Output 

#test[[1]] - Complete genotypes for RIAILS (1 = N2, 0 = CB4856)
#test[[2]] - data for NIL block
#test[[3]] - data for NIL block plot
#test[[4]] - NIL block
#test[[5]] - NIL breaks
#test[[6]] - data for NIL break plot

kammengaNILS <- function(CI.L, CI.R, dir, chromosome){
  nils <- fread(dir)
  
  nils1 <- nils %>%
    gather(strain,geno,-V1)%>%
    select(-N2,-wn009,  -wn037,  -wn055,  -wn106)%>%
    separate(V1, into = c("chr","pos"),sep="bin")%>%
    select(-chr)%>%
    mutate(chr = str_extract(pos,"[[:upper:]]+"))%>%
    mutate(pos1 = str_extract(pos,"[[:digit:]]+"))%>%
    mutate(pos2 = as.numeric(pos1))%>%
    mutate(geno1 = ifelse(geno == 1, geno, 
                          ifelse(geno == -1, 0, NA)))%>%
    select(chr,pos=pos2,strain,geno=geno1)%>%
    filter(!is.na(geno))%>%
    filter(!grepl("wn",strain))%>%
    as.data.frame()%>%
    group_by(strain,chr)%>%
    mutate(index = seq(1:n()))
  
  
  Lpos <- CI.L
  Rpos <- CI.R
  
  nils2 <- nils1%>%
    filter(chr == chromosome)%>%
    filter(pos<Rpos)%>%
    filter(pos>Lpos)%>%
    group_by(strain)%>%
    mutate(CBid = ifelse(sum(geno)<n(),1,0))%>%
    filter(CBid == 1)
  
  if(nrow(nils2)==0){
    nils2 <- nils1%>%
      filter(chr == chromosome)%>%
      filter(pos<Rpos+4e5 & pos>Lpos-4e5)%>%
      group_by(strain)%>%
      mutate(CBid = ifelse(sum(geno)<n(),1,0))%>%
      filter(CBid == 1)
  }
  
  if(nrow(nils2)==0){
    print("NO USEFUL KAMMENGA NILS")
  }else{
    strains <- as.character(unique(nils2$strain))
    
    output <-  nils1 %>%
      filter(strain %in% strains)%>%
      filter(chr == chromosome)%>%
      mutate(intL = which.min(abs(pos-Lpos)),
             intR = which.min(abs(pos-Rpos)))%>%
      ggplot(.)+
      aes(x=index,y = strain,fill=ifelse(geno==1,"N2","CB4856"))+
      geom_tile()+
      geom_vline(aes(xintercept=min(intL)))+
      geom_vline(aes(xintercept=max(intR)))+
      scale_fill_manual(values=c("N2"="orange","CB4856"="blue"),name = "Genotype")+
      theme(axis.text.x = element_text(size=16, face="bold", color="black"),
            axis.text.y = element_text(size=16, face="bold", color="black"),
            axis.title.x = element_text(size=20, face="bold", color="black"),
            axis.title.y = element_text(size=20, face="bold", color="black"),
            strip.text.x = element_text(size=20, face="bold", color="black"),
            strip.text.y = element_text(size=20, face="bold", color="black"),
            plot.title = element_text(size=24, face="bold"))+
      labs(x="Marker Index", y= "Strain", title="Available Kammenga NILs")
  }
  return(output)
}
  
checkNILS <- function(CI.L, CI.R, chromosome, dir, wiggle){
  nils <- fread(dir)
  
  nils1 <- nils %>%
    filter(chr==chromosome)%>%
    rename(end=stop)
  
  if(nrow(nils1)>0){
    nils2 <- nils1%>%
      mutate(identify = ifelse((CI.L<end &CI.L>start),1,
                               ifelse((CI.R<end&CI.R>start), 1,
                                      ifelse((CI.L<(end+wiggle) &CI.L>(start-wiggle)),1,
                                             ifelse((CI.R<(end+wiggle) &CI.R>(start-wiggle)),1,0)))))%>%
      filter(identify ==1)
    if(nrow(nils2)>0){
      ggplot(nils2)+
        aes(y=strain,color = ifelse(nil_genotype=="N2","N2",
                                    ifelse(nil_genotype=="CB4856","CB",NA)))+
        scale_color_manual(values=c("N2"="orange","CB"="blue"),name="Genotype")+
        geom_segment(mapping= aes(x = start/1e6, xend = end/1e6, yend=strain), alpha=1, size = 2)+
        xlim(0,20)+
        geom_vline(xintercept = CI.L/1e6,linetype="dashed")+
        geom_vline(xintercept = CI.R/1e6,linetype="dashed")+
        theme(axis.text.x = element_text(size=16, face="bold", color="black"),
              axis.text.y = element_text(size=16, face="bold", color="black"),
              axis.title.x = element_text(size=20, face="bold", color="black"),
              axis.title.y = element_text(size=20, face="bold", color="black"),
              strip.text.x = element_text(size=20, face="bold", color="black"),
              strip.text.y = element_text(size=20, face="bold", color="black"),
              plot.title = element_text(size=24, face="bold"))+
        labs(x = "Genomic Position", y= "Strain", title = "Available NILs")
    }else{
      print("NO NILS WITHIN 500Kb")
    }
  }else{
#     print("NO NILS")
  }
}

# function to find RIAILs for generating NILs

# input
# CI.L <- Left marker e.g. "UCE1-526"
#CI.R <-  Right marker e.g. "CE1-108"
# chromosome <- chromosome # roman numeral e.g. "I"

# outputs a list
#[[1]] - CompleteInterval1 - dataframe of strains with a continous genotype block across user input interval
# strains also have a breakpoint in the surrounding region
# used to generate plot "complete"
#[[2]] - BreakInterval - dataframe of strains with breaks in the interval and outside the interval
# used to generate plot "complete"
#[[3]] - breaks - ggplot object of all RIAILs identified to have break points in and around the interval
#[[4]] - complete - ggplot object of all RIAILs identified to have a continous genotype stretch across the interval

FindRIAILSforNILs <- function(CI.L, CI.R, chromosome){
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  genos <- data.frame(fread("~/Dropbox/AndersenLab/RCode/NILs/Ancillary/N2xCB4856_598RIAILs_gen.csv",
                            header = T))
  
  load("~/Dropbox/AndersenLab/RCode/NILs/Ancillary/marker_pos_conversion.Rda")
  
  test <- do.call(cbind, lapply(genos[,4:ncol(genos)], function(x){
    temp <- data.frame(gen = x)
    temp <- mutate(temp, num = ifelse(gen =="AA",1,
                                      ifelse(gen=="AB",0,NA)))
    temp <- select(temp, -gen)
  }))
  
  test1 <- data.frame(genos[,1:3], test)
  colnames(test1) <- c("id","chr","CM",seq(1,ncol(test1)-3))
  conv$id <- as.character(conv$id)
  
  test1 <- left_join(test1, conv, by = "id")

  interval <- c(CI.L,CI.R)
  
  # rown <- which(test1$id%in%interval)
  rown <- which(test1$pos %in% interval)
  
  
  
  startDist <- 10
  CompleteInterval1 <- data.frame()
  
  while(length(unique(CompleteInterval1$riail)) < 20){
    
    if(rown[1] > startDist & rown[2] < 1444)
    {
      
      markers <- test1 %>%
        slice(rown[1]:rown[2])%>%
        select(id)
      
      plumin <- c(rown[1]-startDist,rown[2]+startDist)
      
      leftmarks <- test1 %>%
        slice(plumin[1]:(rown[1]-1))%>%
        select(id)
      
      rightmarks <- test1 %>%
        slice((rown[2]+1):plumin[2])%>%
        select(id)
      
      intervals <- test1 %>%
        slice(rown[1]:rown[2])%>%
        select(pos)
      
      ints <- c(min(intervals$pos), max(intervals$pos))
      
      CompleteInterval <- test1 %>%
        slice(plumin[1]:plumin[2])%>%
        mutate(interval = ifelse(id%in%markers$id, "interval",
                                 ifelse(id%in%leftmarks$id, "left",
                                        ifelse(id%in%rightmarks$id, "right", NA))))%>%
        gather(riail, geno, -id ,-chr,-pos,  -CM,-interval)%>%
        group_by(riail,interval)%>%
        mutate(block = sum(geno, na.rm=T))%>%
        ungroup()%>%
        group_by(riail)%>%
        mutate(inte = ifelse((interval == "interval") & (block == 0 | block == nrow(markers)), "complete", 
                             ifelse((interval == "interval") & (block != 0 | block != nrow(markers)), "break",NA)),
               left = ifelse((interval == "left") & (block == 0 | block == nrow(leftmarks)), "completeL", 
                             ifelse((interval == "left") & (block != 0 | block != nrow(leftmarks)), "breakL", NA)),
               right = ifelse( (interval == "right") & (block == 0 | block == nrow(rightmarks)), "completeR",
                               ifelse((interval == "right") & (block != 0 | block != nrow(rightmarks)), "breakR",NA)))%>%
        # ungroup()%>%
        mutate(inte1 = rep(unique(grep("[:alpha:]",inte, value = T))),
               left1 = rep(unique(grep("[:alpha:]",left, value = T))),
               right1 = rep(unique(grep("[:alpha:]",right, value = T))))%>%
        select(-inte,-left,-right)%>%
        rename(inte = inte1, left = left1, right = right1)%>%
        ungroup()%>%
        mutate(bothBreak = ifelse(left == "breakL" & right == "breakR", "breakB",NA))%>%
        mutate(cbgen = ifelse(geno==1,0,
                              ifelse(geno==0,1,NA)),
               Lint = ints[1],
               Rint = ints[2])  
      
      
      
      CompleteInterval1 <- CompleteInterval %>%
        filter(bothBreak == "breakB" & inte == "complete")
      
      startDist <- startDist + 5
      if(startDist > 40){
        break
      }
    }
    else if(rown[1] <= startDist)
    {
      print("You are at the start of Chromosome I")
      markers <- test1 %>%
        slice(rown[1]:rown[2])%>%
        select(id)
      
      plumin <- c(1,rown[2]+startDist)
      
      leftmarks <- test1 %>%
        slice(plumin[1]:(rown[1]-1))%>%
        select(id)
      
      rightmarks <- test1 %>%
        slice((rown[2]+1):plumin[2])%>%
        select(id)
      
      intervals <- test1 %>%
        slice(rown[1]:rown[2])%>%
        select(pos)
      
      ints <- c(min(intervals$pos), max(intervals$pos))
      
      CompleteInterval <- test1 %>%
        slice(plumin[1]:plumin[2])%>%
        mutate(interval = ifelse(id%in%markers$id, "interval",
                                 ifelse(id%in%leftmarks$id, "left",
                                        ifelse(id%in%rightmarks$id, "right", NA))))%>%
        gather(riail, geno, -id ,-chr,-pos,  -CM,-interval)%>%
        group_by(riail,interval)%>%
        mutate(block = sum(geno, na.rm=T))%>%
        ungroup()%>%
        group_by(riail)%>%
        mutate(inte = ifelse((interval == "interval") & (block == 0 | block == nrow(markers)), "complete", 
                             ifelse((interval == "interval") & (block != 0 | block != nrow(markers)), "break",NA)),
               left = ifelse((interval == "left") & (block == 0 | block == nrow(leftmarks)), "completeL", 
                             ifelse((interval == "left") & (block != 0 | block != nrow(leftmarks)), "breakL", NA)),
               right = ifelse( (interval == "right") & (block == 0 | block == nrow(rightmarks)), "completeR",
                               ifelse((interval == "right") & (block != 0 | block != nrow(rightmarks)), "breakR",NA)))%>%
        # ungroup()%>%
        mutate(inte1 = rep(unique(grep("[:alpha:]",inte, value = T))),
               left1 = rep(unique(grep("[:alpha:]",left, value = T))),
               right1 = rep(unique(grep("[:alpha:]",right, value = T))))%>%
        select(-inte,-left,-right)%>%
        rename(inte = inte1, left = left1, right = right1)%>%
        ungroup()%>%
        mutate(bothBreak = ifelse(left == "breakL" & right == "breakR", "breakB",NA))%>%
        mutate(cbgen = ifelse(geno==1,0,
                              ifelse(geno==0,1,NA)),
               Lint = ints[1],
               Rint = ints[2])  
      
      
      
      CompleteInterval1 <- CompleteInterval %>%
        filter(right == "breakR" & inte == "complete")
      
      startDist <- startDist + 5
      if(startDist > 40){
        break
      }
    }
    else if(rown[2] > (1454 - startDist) )
    {
      print("You are at the end of Chromosome X")
      markers <- test1 %>%
        slice(rown[1]:rown[2])%>%
        select(id)
      
      plumin <- c(rown[1]-startDist,1454)
      
      leftmarks <- test1 %>%
        slice(plumin[1]:(rown[1]-1))%>%
        select(id)
      
      rightmarks <- test1 %>%
        slice((rown[2]+1):plumin[2])%>%
        select(id)
      
      intervals <- test1 %>%
        slice(rown[1]:rown[2])%>%
        select(pos)
      
      ints <- c(min(intervals$pos), max(intervals$pos))
      
      CompleteInterval <- test1 %>%
        slice(plumin[1]:plumin[2])%>%
        mutate(interval = ifelse(id%in%markers$id, "interval",
                                 ifelse(id%in%leftmarks$id, "left",
                                        ifelse(id%in%rightmarks$id, "right", NA))))%>%
        gather(riail, geno, -id ,-chr,-pos,  -CM,-interval)%>%
        group_by(riail,interval)%>%
        mutate(block = sum(geno, na.rm=T))%>%
        ungroup()%>%
        group_by(riail)%>%
        mutate(inte = ifelse((interval == "interval") & (block == 0 | block == nrow(markers)), "complete", 
                             ifelse((interval == "interval") & (block != 0 | block != nrow(markers)), "break",NA)),
               left = ifelse((interval == "left") & (block == 0 | block == nrow(leftmarks)), "completeL", 
                             ifelse((interval == "left") & (block != 0 | block != nrow(leftmarks)), "breakL", NA)),
               right = ifelse( (interval == "right") & (block == 0 | block == nrow(rightmarks)), "completeR",
                               ifelse((interval == "right") & (block != 0 | block != nrow(rightmarks)), "breakR",NA)))%>%
        # ungroup()%>%
        mutate(inte1 = rep(unique(grep("[:alpha:]",inte, value = T))),
               left1 = rep(unique(grep("[:alpha:]",left, value = T))),
               right1 = rep(unique(grep("[:alpha:]",right, value = T))))%>%
        select(-inte,-left,-right)%>%
        rename(inte = inte1, left = left1, right = right1)%>%
        ungroup()%>%
        mutate(bothBreak = ifelse(left == "breakL" & right == "breakR", "breakB",NA))%>%
        mutate(cbgen = ifelse(geno==1,0,
                              ifelse(geno==0,1,NA)),
               Lint = ints[1],
               Rint = ints[2])  
      
      
      CompleteInterval1 <- CompleteInterval %>%
        filter(left == "breakL" & inte == "complete")
      
      startDist <- startDist + 5
      if(startDist > 40){
        break
      }
    }
  }
  
  if(length(unique(CompleteInterval1$chr)) > 1){
    CompleteInterval1 <- CompleteInterval1 %>%
      filter(chr == chromosome)
    print("You are at the end of a chromosome")
  }
  
  BreakInterval <- CompleteInterval %>%
    filter((left == "breakL" | right == "breakR") & inte == "break")
  
  if(length(unique(BreakInterval$chr)) > 1){
    BreakInterval <- BreakInterval %>%
      filter(chr == chromosome)
    print("You are at the end of a chromosome")
  }
  
  breaks <- ggplot(BreakInterval)+
    aes(x = pos/1e6, y = geno)+
    facet_grid(riail~chr, scales = "free_y")+
    geom_ribbon(aes(ymin = 0, ymax = geno, alpha=.5), fill = "orange")+
    geom_ribbon(aes(ymin = 0, ymax = cbgen, alpha=.5), fill = "blue")+
    geom_vline(aes(xintercept = Lint/1e6))+
    geom_vline(aes(xintercept = Rint/1e6))+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=0, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black"),
          axis.title.y = element_text(size=20, face="bold", color="black"),
          strip.text.x = element_text(size=20,face="bold", color="black"),
          strip.text.y = element_text(size=20, angle =0, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold"),
          legend.position = "none")+
    labs(x = "Genomic Position (Mb)", y = "RIAIL")
  
  complete <- ggplot(CompleteInterval1)+
    aes(x = pos/1e6, y = geno)+
    facet_grid(riail~chr, scales = "free_y")+
    geom_ribbon(aes(ymin = 0, ymax = geno, alpha=.5), fill ="orange")+
    geom_ribbon(aes(ymin = 0, ymax = cbgen, alpha=.5), fill ="blue")+
    geom_vline(aes(xintercept = Lint/1e6))+
    geom_vline(aes(xintercept = Rint/1e6))+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=0, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black"),
          axis.title.y = element_text(size=20, face="bold", color="black"),
          strip.text.x = element_text(size=20,face="bold", color="black"),
          strip.text.y = element_text(size=20, angle =0, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold"),
          legend.position = "none")+
    labs(x = "Genomic Position (Mb)", y = "RIAIL")
  
  return(list(CompleteInterval1, BreakInterval, breaks, complete))
}        

FindRIAILSforNILs(9093651, 11342635, "X")
FindRIAILSforNILs(5274023, 7727638, "X")

# linkage PxG function

# PxGlinkage("UCE5-1890","bortezomib_resid.q90.TOF")

PxGlinkage <- function(peakMarker, condtrt){
  
  
  temp <- str_split_fixed(condtrt,pattern="_",n=2)
  trt <- temp[,2]
  condition <- temp[,1]
  
  
  gen <- Lgeno %>%
    filter(snp==peakMarker)
  
  posit <- gen[1,5]
  
  Lphen%>%
    filter(grepl(condition,trait))%>%
    filter(grepl(trt,trait))%>%
    left_join(.,gen,by="strain")%>%
    mutate(color = ifelse(strain=="N2", "N2",
                          ifelse(strain == "CB4856", "CB4856",
                                 ifelse(geno == 1, "RIAILs-N2",
                                        ifelse(geno == 2, "RIAILs-CB",0)))),
           label = ifelse(strain=="N2", "N2",
                          ifelse(strain == "CB4856", "CB4856",
                                 ifelse(geno == 1, "RIAILs-N2",
                                        ifelse(geno == 2, "RIAILs-CB",0)))))%>%
    filter(!is.na(color))%>%
    ggplot(.)+
    aes(x = color, y = value,fill = ifelse(strain=="N2", "N2",
                                           ifelse(strain == "CB4856", "CB4856",
                                                  ifelse(geno == 1, "RIAILs-N2",
                                                         ifelse(geno == 2, "RIAILs-CB",0)))),
        color = ifelse(strain=="N2", "N2",
                       ifelse(strain == "CB4856", "CB4856",
                              ifelse(geno == 1, "RIAILs-N2",
                                     ifelse(geno == 2, "RIAILs-CB",0)))))+
    scale_fill_manual(values = c("N2"="orange","CB4856"="blue","RIAILs-CB"="blue","RIAILs-N2"="orange"), name = "Genotype")+
    scale_color_manual(values = c("N2"="black","CB4856"="black","RIAILs-CB"="#666666","RIAILs-N2"="#666666"),name="Genotype")+
    geom_boxplot(outlier.shape=NA)+
    theme_bw()+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=16, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black"),
          axis.title.y = element_text(size=20, face="bold", color="black",vjust=1),
          strip.text.x = element_text(size=20, face="bold", color="black"),
          strip.text.y = element_text(size=20, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold",vjust=1),
          legend.title = element_text(size=14),
          panel.border = element_rect(size=2))+
    geom_jitter(alpha=.75)+
    labs(title=paste( trt,peakMarker,posit,sep="_"),y=trt, x="Genotype" )
}

# function takes linkage mapping data and plots trait of interest

# input
# condition_trait - this is the condition and trait of the mapping of interest separated by a "_"

# example input::

#linkagePlot("bleomycin_resid.q75.TOF")

linkagePlot <- function(condition_trait){
  #linkage plot data
  linkplot2 <- linkage%>%
    filter(pheno == condition_trait)%>%
    group_by(SNP)%>%
    arrange(desc(LOD))%>%
    mutate(maxlod = max(LOD))%>% #create maxlod column with max lod
    mutate(QTLleft = max(CI.L.pos,na.rm=TRUE))%>%
    mutate(QTLright = max(CI.R.pos,na.rm=TRUE))%>%
    mutate(varexp = max(var.exp,na.rm=TRUE))%>% #create varexp column with whatever iteration has it calculated
    distinct(SNP)
  #filter for list of QTL peaks
  linkpeaks <- linkplot2%>%
    filter(!is.na(var.exp))
  #add LOD2 column containing LOD values only for CI ranges and 0 for everything else
  linkplot2 <- mutate(linkplot2, LOD2 = 0) #initialize column
  for(i in 1:nrow(linkpeaks)){
    linkplot2$LOD2 <- ifelse(linkplot2$pos >= linkpeaks$QTLleft[i] & linkplot2$pos <= linkpeaks$QTLright[i] & linkplot2$chr == linkpeaks$chr[i], linkplot2$LOD, linkplot2$LOD2)
  }
  
  p1a<- ggplot(linkplot2)+
    aes(x=pos/1e6,y=LOD)+
    geom_line(size = 1, alpha = 0.85)+
    facet_grid(.~chr, scales ="free")+
    #ylim (0, 16) +
    labs(x = "Position (Mb)", y = "LOD score") +
    # theme_bw() + 
    geom_point(data=subset(linkplot2, linkplot2$varexp > 0 ), aes(x=pos/1e6,y=maxlod + 0.75), fill ="red", shape=25, size=3.2, show_guide = FALSE) +
    geom_text(data=subset(linkplot2, linkplot2$varexp > 0), aes(x=pos/1e6,y=maxlod + 4,label= paste0(100*round(varexp, digits = 4),"%")), size=5) +
    geom_hline(aes(yintercept=3), color = "red", linetype="dashed", size=0.5) +
    geom_ribbon(aes(ymin = 0, ymax = LOD2), fill = "blue", alpha = 0.5)+
    theme_bw()
  
  p1a <- p1a +   theme(axis.text.x = element_text(size=24, face="bold", color="black"),
                       axis.text.y = element_text(size=24, face="bold", color="black"),
                       axis.title.x = element_text(size=24, face="bold", color="black", vjust=-.3),
                       axis.title.y = element_text(size=24, face="bold", color="black"),
                       strip.text.x = element_text(size=24, face="bold", color="black"),
                       strip.text.y = element_text(size=24, face="bold", color="black"),
                       plot.title = element_text(size=24, face="bold", vjust = 1),
                       legend.position="none",
                       panel.background = element_rect( color="black",size=1.2))
  
  p1a
}


# plot riail genotype function
# # etoposide
# plotRIAILgenotypes(c(103,322,327,462), "II", 11492171, 12008179)

plotRIAILgenotypes <- function(riails, ch, leftPos, rightPos){
  
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  genos <- data.frame(fread("~/Dropbox/AndersenLab/RCode/NILs/Ancillary/N2xCB4856_598RIAILs_gen.csv",
                            header = T))
  
  load("~/Dropbox/AndersenLab/RCode/NILs/Ancillary/marker_pos_conversion.Rda")
  
  test <- do.call(cbind, lapply(genos[,4:ncol(genos)], function(x){
    temp <- data.frame(gen = x)
    temp <- mutate(temp, num = ifelse(gen =="AA",1,
                                      ifelse(gen=="AB",0,NA)))
    temp <- select(temp, -gen)
  }))
  
  test1 <- data.frame(genos[,1:3], test)
  colnames(test1) <- c("id","chr","CM",seq(1,ncol(test1)-3))
  
  test1 <- left_join(test1, conv, by = "id")
  
  pl <- test1 %>%
    gather(riail,geno, -id, -chr,-CM,-pos)%>%
    filter(riail %in% riails)%>%
    mutate(cbgen = ifelse(geno==1,0,
                          ifelse(geno==0,1,NA)),
           Lint = leftPos,
           Rint = rightPos)%>%
    filter(chr == ch)%>%
    ggplot(.)+
    aes(x = pos/1e6, y = geno)+
    facet_grid(riail~chr, scales = "free_y")+
    geom_ribbon(aes(ymin = 0, ymax = geno, alpha=.5), fill ="orange")+
    geom_ribbon(aes(ymin = 0, ymax = cbgen, alpha=.5), fill ="blue")+
    geom_vline(aes(xintercept = Lint/1e6))+
    geom_vline(aes(xintercept = Rint/1e6))+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=0, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black"),
          axis.title.y = element_text(size=20, face="bold", color="black"),
          strip.text.x = element_text(size=20,face="bold", color="black"),
          strip.text.y = element_text(size=20, angle =0, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold"),
          legend.position = "none")+
    labs(x = "Genomic Position (Mb)", y = "RIAIL")
  
  df <- test1 %>%
    gather(riail,geno, -id, -chr,-CM,-pos)%>%
    filter(riail %in% riails)%>%
    filter(chr == ch)%>%
    # select(riail, pos, geno)%>%
    spread(riail,geno)%>%
    arrange(pos)  %>%
    mutate( Lint = leftPos,
            Rint = rightPos)
  
  return(list(pl,df))
}


# Example inputs:

# showQTLintervals(1000,15000000,"I",NA,NA)
# showQTLintervals(1000,15000000,"II",NA,NA)
# showQTLintervals(1000,15000000,"III",NA,NA)
# showQTLintervals(1000,17000000,"IV",NA,NA)
# showQTLintervals(1000,20000000,"V",NA,NA)
# showQTLintervals(1000,20000000,"X",NA,NA)
# showQTLintervals(15400000,15600000,"IV",NA,NA)
# showQTLintervals(15400000,15600000,"I","Chemotherapeutic",NA)
# showQTLintervals(1000,14779486,"I","Chemotherapeutic",NA)
# showQTLintervals(1000,14779486,"I","Chemotherapeutic",c("topotecan",
#                                                         "vincristine",
#                                                         "irinotecan",
#                                                         "bortezomib"))
# 
# showQTLintervals(8000000,13000000,"V","Chemotherapeutic",NA)
# 
# 
# showQTLintervals(1000,14779486,"III",NA,NA)
# showQTLintervals(1000,14779486,"IV",NA,NA)
# showQTLintervals(1000,14779486,"I","Chemotherapeutic",NA)
# 
# showQTLintervals(1000,20000000,"V","Neuroactive",NA)


showQTLintervals <- function(left, 
                             right, 
                             chromosome, 
                             treatments,
                             traits){
  
  library(dplyr)
  library(ggplot2)
  
  load("~/Dropbox/AndersenLab/RCode/NILs/Data/LinkagePeaks.Rda")
  
  df <- linkPEAKS
  
  if(is.na(treatments)){
    df %>%
      filter(!is.na(var.exp))%>%
      filter(chr == chromosome)%>%
      filter(pos < right+1e6)%>%
      filter(pos > left-1e6)%>%
      mutate(lint = left,
             rint = right)%>%
      group_by(condition) %>%
      filter(LOD == max(LOD))%>%
      ggplot(.)+
      aes(x = pos/1e6, y = LOD, fill = condition)+
      geom_rect(aes(xmin = CI.L.pos/1e6, xmax=CI.R.pos/1e6, ymax = LOD+1, ymin = LOD-1), color = "black", alpha = .5)+
      geom_vline(aes(xintercept = lint/1e6), size = 1)+
      geom_vline(aes(xintercept = rint/1e6), size = 1)+
      facet_grid(treatment_type~., scales = "free")+
      labs(y = "LOD", x = "Genomic Position (Mb)")+
      theme(axis.text.x = element_text(size=16, face="bold", color="black"),
            axis.text.y = element_text(size=16, face="bold", color="black"),
            axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
            axis.title.y = element_text(size=20, face="bold", color="black"),
            strip.text.x = element_text(size=20, face="bold", color="black"),
            strip.text.y = element_text(size=20, face="bold", color="black", angle = 0),
            plot.title = element_text(size=24, face="bold", vjust = 1))
    
  }else if(is.na(traits)){
    df %>%
      filter(!is.na(var.exp))%>%
      filter(chr == chromosome)%>%
      filter(treatment_type == treatments)%>%
      filter(pos < right+1e6)%>%
      filter(pos > left-1e6)%>%
      mutate(lint = left,
             rint = right)%>%
      #       mutate(roundedpos = signif(pos,digits=1))%>%
      group_by(condition) %>%
      filter(LOD == max(LOD))%>%
      ggplot(.)+
      aes(x = pos/1e6, y = LOD, fill = condition)+
      geom_rect(aes(xmin = CI.L.pos/1e6, xmax=CI.R.pos/1e6, ymax = LOD+1, ymin = LOD-1), color = "black", alpha = .5)+
      geom_point(aes(x=pos/1e6,y=LOD))+
      geom_text(aes(x=pos/1e6,y=LOD, label=signif(var.exp*100,digits=2)),size=6,hjust=1.5)+
      geom_vline(aes(xintercept = lint/1e6), size = 1)+
      geom_vline(aes(xintercept = rint/1e6), size = 1)+
      facet_grid(condition~., scales = "free")+
      labs(y = "LOD", x = "Genomic Position (Mb)")+
      theme(axis.text.x = element_text(size=16, face="bold", color="black"),
            axis.text.y = element_text(size=10, face="bold", color="black"),
            axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
            axis.title.y = element_text(size=20, face="bold", color="black"),
            strip.text.x = element_text(size=20, face="bold", color="black"),
            strip.text.y = element_text(size=20, face="bold", color="black", angle = 0),
            plot.title = element_text(size=24, face="bold", vjust = 1),
            legend.position="none")
  }else{
    df %>%
      filter(!is.na(var.exp))%>%
      filter(chr == chromosome)%>%
      filter(treatment_type == treatments)%>%
      filter(condition %in% traits)%>%
      filter(pos < right+1e6)%>%
      filter(pos > left-1e6)%>%
      mutate(lint = left,
             rint = right)%>%
      mutate(roundedpos = signif(pos,digits=1))%>%
      group_by(condition,roundedpos) %>%
      filter(LOD==max(LOD))%>%
      ggplot(.)+
      aes(x = pos/1e6, y = LOD, fill = condition)+
      geom_rect(aes(xmin = CI.L.pos/1e6, xmax=CI.R.pos/1e6, ymax = LOD+1, ymin = LOD-1), color = "black", alpha = .5)+
      geom_point(aes(x=pos/1e6,y=LOD))+
      geom_text(aes(x=pos/1e6,y=LOD, label=signif(var.exp*100,digits=2)),size=6,hjust=1.5)+
      geom_vline(aes(xintercept = lint/1e6), size = 1)+
      geom_vline(aes(xintercept = rint/1e6), size = 1)+
      facet_grid(trait~., scales = "free")+
      labs(y = "LOD", x = "Genomic Position (Mb)")+
      theme(axis.text.x = element_text(size=16, face="bold", color="black"),
            axis.text.y = element_text(size=10, face="bold", color="black"),
            axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
            axis.title.y = element_text(size=20, face="bold", color="black"),
            strip.text.x = element_text(size=20, face="bold", color="black"),
            strip.text.y = element_text(size=20, face="bold", color="black", angle = 0),
            plot.title = element_text(size=24, face="bold", vjust = 1))
  }
  
}

# showQTLintervals(6000000,13000000,"V","Chemotherapeutic",NA)


