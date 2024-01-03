library(ggplot2)
library(ggalluvial)
library(ggpubr)
library(ggtern)
library(tidyverse)
library(gridExtra)
pops<-c("CLM","MXL","PEL","PUR","ACB","ASW")
anc1_prop<-c(82,43,32,138,881,765)
anc2_prop<-c(277,500,761,150,5,41)
pulses_false_true<-c("false","true")

##functions
plotsfunction <- function(run){
  
  #####ANCESTRIES PROP PLOT
  ancestries_plot <-  ancestries_ind_df_allruns[which( ancestries_ind_df_allruns$run==run),]
  ancestries_plot <-  ancestries_ind_df_allruns
  
  ggplot(ancestries_plot, aes(x = ancestries_plot$generation, y = ancestries_plot$child.prop_aut, fill=ancestries_plot$ancestry, colour=ancestries_plot$ancestry)) + 
    stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.25), fun.max = function(x) quantile(x, 0.75),alpha=0.25,colour = NA) + 
    stat_summary(geom="line", size=0.5,fun=median)+
    xlab("Generation")+ xlim(3,GEN)+
    ylab("Fraction of ancestry")+
    scale_color_manual("Median and quartiles of ancestry",values=c("#b84a62","#bdb246","#437c90"))+
    scale_fill_manual("Median and quartiles of ancestry",values=c("#b84a62","#bdb246","#437c90"))+
    theme(legend.position="none",plot.title = element_text(size = 7),axis.title=element_text(size=6),axis.text=element_text(size=5))+
    facet_wrap(~paste("AM:",AM1,AM2,AM3,"; SB:",SB1,",",SB2,"; POP1",POP1,"POP2",POP2,"POP3",POP3,"; run",run, sep=" "))
}



plotsmatingbarsfunction <- function(run,g){
  ancestries_ind_df_allruns <- read.table(file=paste(path,file_anc,sep=""),header = TRUE,sep=";")
  ancestries_ind_df_allruns$run<-run
  ###COUPLES PLOT
  ancestries_plot <-  ancestries_ind_df_allruns[which(ancestries_ind_df_allruns$run==run),]
  mating <- ancestries_plot[which(ancestries_plot$generation==g),]
  #colnames(mating) <- c("generation","ancestry","child.index","child.sex","child.prop_aut","child.prop_X","child.prop_Y","mother.index","mother.prop","father.index","father.prop","run","AM1","AM2","AM3","SB1","SB2","POP1","POP2","POP3","totalGEN")
  mating$ancestry <- as.factor(as.character(mating$ancestry))
  mating$mother.index <- as.factor(as.character(mating$mother.index))
  mating$father.index <- as.factor(as.character(mating$father.index))
  mating$mother.prop <- as.numeric(as.character(mating$mother.prop))
  mating$father.prop <- as.numeric(as.character(mating$father.prop))
  ##sort individualb by ancestry. This is the only purpose of matingtogether dataframe, to get sorted_mothers and sorted_fathers.
  matingtogether <- data.frame(matrix(nrow=length(mating$generation)/3,ncol=0))
  for (anc in c(1:3)){
    mat_anc <- mating[which(mating$ancestry==anc),-c(1:4,6,8)]
    names(mat_anc) <- paste(names(mat_anc),anc,sep=".")
    matingtogether<-cbind(matingtogether,mat_anc)
  }
  matingtogether$mother.index<-as.character(paste(mating$mother.index[which(mating$ancestry==1)],".mother",sep=""))
  matingtogether$father.index<-as.character(paste(mating$father.index[which(mating$ancestry==1)],".father",sep=""))
  matingtogether$mother.maxanc <- as.factor(as.character(sapply(c(1:length(matingtogether$mother.prop.1)),function(i) which.max(c(matingtogether$mother.prop.1[i],matingtogether$mother.prop.2[i],matingtogether$mother.prop.3[i])))))
  matingtogether$mother.maxprop <- sapply(c(1:length(matingtogether$mother.prop.1)),function(i) max(c(matingtogether$mother.prop.1[i],matingtogether$mother.prop.2[i],matingtogether$mother.prop.3[i])))
  matingtogether$father.maxanc <- as.factor(as.character(sapply(c(1:length(matingtogether$father.prop.1)),function(i) which.max(c(matingtogether$father.prop.1[i],matingtogether$father.prop.2[i],matingtogether$father.prop.3[i])))))
  matingtogether$father.maxprop<- sapply(c(1:length(matingtogether$father.prop.1)),function(i) max(c(matingtogether$father.prop.1[i],matingtogether$father.prop.2[i],matingtogether$father.prop.3[i])))
  sorted_mothers <- matingtogether$mother.index[order(matingtogether$mother.maxanc,-matingtogether$mother.maxprop)]
  sorted_fathers <- matingtogether$father.index[order(matingtogether$father.maxanc,-matingtogether$father.maxprop)]
  ##keep ordering dataset, need single ancestry column for both father and mother,etc.
  mating$couple <-  as.factor(as.character(sort(rep(c(1:(length(mating$generation)/3)),3))))
  mating$mother.maxanc <-  unlist(lapply(matingtogether$mother.maxanc,function(x)rep(x,3)))
  mating$father.maxanc <-  unlist(lapply(matingtogether$father.maxanc,function(x)rep(x,3)))
  mating$mother.maxprop <-  unlist(lapply(matingtogether$mother.maxprop,function(x)rep(x,3)))
  mating$father.maxprop<-  unlist(lapply(matingtogether$father.maxprop,function(x)rep(x,3)))
  mating$mother.index <- factor(paste(mating$mother.index,".mother",sep=""),levels=unique(sorted_mothers))
  mating$father.index <- factor(paste(mating$father.index,".father",sep=""),levels=unique(sorted_fathers))
  matingmother <- mating[,c(c("generation","couple","ancestry"),names(mating)[grep("mother",names(mating))])]
  matingfather <- mating[,c(c("generation","couple","ancestry"),names(mating)[grep("father",names(mating))])]
  matingmother$parent <- "mother"
  matingmother$parent <- factor(matingmother$parent,levels=c("mother","father"))
  matingfather$parent <- "father"
  matingfather$parent <- factor(matingfather$parent,levels=c("mother","father"))
  names(matingmother) <- gsub("mother.","",names(matingmother))
  names(matingfather) <- gsub("father.","",names(matingfather))
  matingmotherfather <- rbind(matingmother,matingfather)
  matingmotherfather <- matingmotherfather[order(matingmotherfather$index),]
  matingmotherfather$id <- factor(as.character(unlist(lapply(rep(1:(length(matingmotherfather$generation)/6),2),function(x) rep(x,3)))),levels=rev(unique(as.character(unlist(lapply(rep(1:(length(matingmotherfather$generation)/6),2),function(x) rep(x,3)))))))
  ##need single row per couple, then filtering only one ancestry.
  matingmotherfather$maxanc <- as.character(matingmotherfather$maxanc)
  mating_anc1 <- matingmotherfather[which(matingmotherfather$ancestry==1),]
  mating_anc1$sharedanc <- unlist(lapply(unique(mating_anc1$couple),function(x) length(unique(mating_anc1$maxanc[which(mating_anc1$couple == x)]))))
  mating_anc1$sharedanc[which(mating_anc1$sharedanc==2)] <- "Different ancestries"
  mating_anc1$sharedanc[which(mating_anc1$sharedanc==1)] <- mating_anc1$maxanc[which(mating_anc1$sharedanc==1)]
  mating_anc1 <- mating_anc1[order(mating_anc1$index),]
  matingmotherfather$maxanc <- factor(as.character(matingmotherfather$maxanc), levels=unique(sort(as.character(mating_anc1$maxanc))))
  mating_anc1$maxanc <- factor(as.character(mating_anc1$maxanc), levels=unique(sort(as.character(mating_anc1$maxanc))))
  colors<-as.data.frame(cbind(c("Different ancestries",1,2,3),c("gray20","#b84a62","#bdb246","#437c90")))
  names(colors) <- c("maxanc","maxanccolor")
  mating_anc1$sharedanc <- as.factor(mating_anc1$sharedanc)
  mating_anc1$sharedanc <- factor(mating_anc1$sharedanc, levels = levels(mating_anc1$sharedanc)[c(grep("Different ancestries",levels(mating_anc1$sharedanc)),grep("Different ancestries",levels(mating_anc1$sharedanc),invert=TRUE))])
  mu <- mating_anc1
  mating_anc1 <<- mu
  fillcolors <<- as.character(unlist(lapply(levels(mating_anc1$sharedanc),function(x) colors$maxanccolor[which(colors$maxanc==x)])))
  gg <- ggplot(mating_anc1, aes(x=mating_anc1$parent, alluvium=mating_anc1$couple,
                                stratum=mating_anc1$index, alpha=mating_anc1$maxprop )) 
  gg_nostratum <- gg + geom_flow(aes(fill=mating_anc1$sharedanc,color=mating_anc1$sharedanc),width=1/12,stat = "alluvium", lode.guidance = "rightward",alpha=0.3) + 
    scale_fill_manual("Ancestry",values=as.character(fillcolors))+
    scale_color_manual("Ancestry",values=as.character(fillcolors))+
    theme(axis.title=element_blank(),
          legend.position="none",
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin=unit(c(0,-3,0,-1), "cm")) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) 
  mating_mother <<- matingmotherfather[which(matingmotherfather$parent=="mother"),]
  ggmother <- ggplot(mating_mother, aes(x=mating_mother$id, y=mating_mother$prop, fill=mating_mother$ancestry)) +
    geom_col() + scale_fill_manual("Ancestry",values=c("#b84a62","#bdb246","#437c90","black")) + 
    theme(axis.title=element_blank(),
          legend.position="none",
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin=unit(c(0,0,0.1,0), "cm"))+
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    coord_flip()
  mating_father <<- matingmotherfather[which(matingmotherfather$parent=="father"),]
  ggfather <- ggplot(mating_father, aes(x=mating_father$id, y=mating_father$prop, fill=mating_father$ancestry)) +
    geom_col() + scale_fill_manual("Ancestry",values=c("#b84a62","#bdb246","#437c90","black")) + 
    theme(axis.title=element_blank(),
          legend.position="none",
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin=unit(c(0,0,0.1,0), "cm"))+
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    coord_flip()
  ggarrange(ggmother,gg_nostratum,ggfather,ncol=3,hjust=0,vjust=0,widths=c(1,3,1))
}

plotsmatingternfunction <- function(run,g){
  ancestries_ind_df_allruns <- read.table(file=paste(path,file_anc,sep=""),header = TRUE,sep=";")
  ancestries_ind_df_allruns$run<-run
  ###MATINGTERNARY
  ancestries_plot <-  ancestries_ind_df_allruns[which( ancestries_ind_df_allruns$run==run),]
  mating <- ancestries_plot[which(ancestries_plot$generation==g),]
  #colnames(mating) <- c("generation","ancestry","child.index","child.sex","child.prop_aut","child.prop_X","child.prop_Y","mother.index","mother.prop","father.index","father.prop","run","AM1","AM2","AM3","SB1","SB2","POP1","POP2","POP3","totalGEN")
  mating$ancestry <- as.factor(as.character(mating$ancestry))
  mating$mother.index <- as.factor(as.character(mating$mother.index))
  mating$father.index <- as.factor(as.character(mating$father.index))
  ##sort individualb by ancestry. This is the only purpose of matingtogether dataframe, to get sorted_mothers and sorted_fathers.
  matingtogether <- data.frame(matrix(nrow=length(mating$generation)/3,ncol=0))
  for (anc in c(1:3)){
    mat_anc <- mating[which(mating$ancestry==anc),-c(1:4,6,8)]
    names(mat_anc) <- paste(names(mat_anc),anc,sep=".")
    matingtogether<-cbind(matingtogether,mat_anc)
  }
  matingtogether$mother.index<-as.character(paste(mating$mother.index[which(mating$ancestry==1)],".mother",sep=""))
  matingtogether$father.index<-as.character(paste(mating$father.index[which(mating$ancestry==1)],".father",sep=""))
  matingtogether$mother.maxanc <- as.factor(as.character(sapply(c(1:length(matingtogether$mother.prop.1)),function(i) which.max(c(matingtogether$mother.prop.1[i],matingtogether$mother.prop.2[i],matingtogether$mother.prop.3[i])))))
  matingtogether$mother.maxprop <- sapply(c(1:length(matingtogether$mother.prop.1)),function(i) max(c(matingtogether$mother.prop.1[i],matingtogether$mother.prop.2[i],matingtogether$mother.prop.3[i])))
  matingtogether$father.maxanc <- as.factor(as.character(sapply(c(1:length(matingtogether$father.prop.1)),function(i) which.max(c(matingtogether$father.prop.1[i],matingtogether$father.prop.2[i],matingtogether$father.prop.3[i])))))
  matingtogether$father.maxprop<- sapply(c(1:length(matingtogether$father.prop.1)),function(i) max(c(matingtogether$father.prop.1[i],matingtogether$father.prop.2[i],matingtogether$father.prop.3[i])))
  sorted_mothers <- matingtogether$mother.index[order(matingtogether$mother.maxanc,-as.numeric(matingtogether$mother.maxprop))]
  sorted_fathers <- matingtogether$father.index[order(matingtogether$father.maxanc,-as.numeric(matingtogether$father.maxprop))]
  ##keep ordering dataset, need single ancestry column for both father and mother,etc.
  mating$couple <-  as.factor(as.character(sort(rep(c(1:(length(mating$generation)/3)),3))))
  mating$mother.maxanc <-  unlist(lapply(matingtogether$mother.maxanc,function(x)rep(x,3)))
  mating$father.maxanc <-  unlist(lapply(matingtogether$father.maxanc,function(x)rep(x,3)))
  mating$mother.maxprop <-  unlist(lapply(matingtogether$mother.maxprop,function(x)rep(x,3)))
  mating$father.maxprop<-  unlist(lapply(matingtogether$father.maxprop,function(x)rep(x,3)))
  mating$mother.index <- factor(paste(mating$mother.index,".mother",sep=""),levels=unique(sorted_mothers))
  mating$father.index <- factor(paste(mating$father.index,".father",sep=""),levels=unique(sorted_fathers))
  matingmother <- mating[,c(c("generation","couple","ancestry"),names(mating)[grep("mother",names(mating))])]
  matingfather <- mating[,c(c("generation","couple","ancestry"),names(mating)[grep("father",names(mating))])]
  matingmother$parent <- "mother"
  matingmother$parent <- factor(matingmother$parent,levels=c("mother","father"))
  matingfather$parent <- "father"
  matingfather$parent <- factor(matingfather$parent,levels=c("mother","father"))
  names(matingmother) <- gsub("mother.","",names(matingmother))
  names(matingfather) <- gsub("father.","",names(matingfather))
  mating <- rbind(matingmother,matingfather)
  mating <- mating[order(mating$index),]
  mating$id <- factor(as.character(unlist(lapply(rep(1:(length(mating$generation)/6),2),function(x) rep(x,3)))),levels=rev(unique(as.character(unlist(lapply(rep(1:(length(mating$generation)/6),2),function(x) rep(x,3)))))))
  ##need single row per couple, then filtering only one ancestry.
  mating$maxanc <- as.character(mating$maxanc)
  mating_anc1 <- mating[which(mating$ancestry==1),]
  mating_anc1$sharedanc <- unlist(lapply(unique(mating_anc1$couple),function(x) length(unique(mating_anc1$maxanc[which(mating_anc1$couple == x)]))))
  mating_anc1$sharedanc[which(mating_anc1$sharedanc==2)] <- "Different ancestries"
  mating_anc1$sharedanc[which(mating_anc1$sharedanc==1)] <- mating_anc1$maxanc[which(mating_anc1$sharedanc==1)]
  mating_anc1 <- mating_anc1[order(mating_anc1$index),]
  mating$maxanc <- factor(as.character(mating$maxanc), levels=c("Different ancestries",unique(sort(as.character(mating_anc1$maxanc)))))
  matingtern <- mating
  matingtern <- matingtern[order(matingtern$couple,matingtern$parent,matingtern$ancestry),]
  matingtern$sharedanc <- unlist(lapply(unique(matingtern$couple),function(x) rep(length(unique(matingtern$maxanc[which(matingtern$couple == x)])),6)))
  matingtern$sharedanc[which(matingtern$sharedanc==2)] <- "Different ancestries"
  matingtern$sharedanc[which(matingtern$sharedanc==1)] <- as.character(matingtern$maxanc[which(matingtern$sharedanc==1)])
  matingtern$sharedanc <- factor(matingtern$sharedanc, levels=c("Different ancestries",unique(sort(as.character(matingtern$maxanc)))))
  matingtern <- matingtern[order(matingtern$parent,matingtern$couple,matingtern$id),]
  matingtern_collapsed <-matingtern[which(matingtern$ancestry=="1"),]
  matingtern_collapsed$prop2 <- matingtern$prop[which(matingtern$ancestry=="2")]
  matingtern_collapsed$prop3 <- matingtern$prop[which(matingtern$ancestry=="3")]
  matingtern_collapsed <- matingtern_collapsed[,c(c(1:2),c(5:10),4,c(11:12))]
  matingtern_collapsed$mother1 <- rep(NA,length(matingtern_collapsed$generation))
  matingtern_collapsed$mother2 <- rep(NA,length(matingtern_collapsed$generation))
  matingtern_collapsed$mother3 <- rep(NA,length(matingtern_collapsed$generation))
  mating_segments<-as.data.frame(cbind(
    generation=as.character(matingtern_collapsed$generation[which(matingtern_collapsed$parent=="mother")]),
    couple=as.character(matingtern_collapsed$couple[which(matingtern_collapsed$parent=="mother")]),
    index=as.character(matingtern_collapsed$index[which(matingtern_collapsed$parent=="mother")]),
    maxanc=as.character(matingtern_collapsed$maxanc[which(matingtern_collapsed$parent=="mother")]),
    maxprop=as.character(matingtern_collapsed$maxprop[which(matingtern_collapsed$parent=="mother")]),
    parent=as.character(matingtern_collapsed$parent[which(matingtern_collapsed$parent=="mother")]),
    id=as.character(matingtern_collapsed$id[which(matingtern_collapsed$parent=="mother")]),
    sharedanc=as.character(matingtern_collapsed$sharedanc[which(matingtern_collapsed$parent=="mother")]),
    prop1=as.character(matingtern_collapsed$prop[which(matingtern_collapsed$parent=="mother")]),
    prop2=as.character(matingtern_collapsed$prop2[which(matingtern_collapsed$parent=="mother")]),
    prop3=as.character(matingtern_collapsed$prop3[which(matingtern_collapsed$parent=="mother")]),
    father1=as.character(matingtern_collapsed$prop[which(matingtern_collapsed$parent=="father")]),
    father2=as.character(matingtern_collapsed$prop2[which(matingtern_collapsed$parent=="father")]),
    father3=as.character(matingtern_collapsed$prop3[which(matingtern_collapsed$parent=="father")])))
  mating_segments$prop1 <-as.numeric(as.character(mating_segments$prop1))
  mating_segments$prop2 <-as.numeric(as.character(mating_segments$prop2))
  mating_segments$prop3 <-as.numeric(as.character(mating_segments$prop3))
  mating_segments$maxprop <-as.numeric(as.character(mating_segments$maxprop))
  mating_segments$father1 <-as.numeric(as.character(mating_segments$father1))
  mating_segments$father2 <-as.numeric(as.character(mating_segments$father2))
  mating_segments$father3 <-as.numeric(as.character(mating_segments$father3))
  mating_segments$sharedanc <- factor(mating_segments$sharedanc,levels=c("Different ancestries","1","2","3"))
  colors <- as.data.frame(cbind(c("Different ancestries","1","2","3"),c("gray40","#b84a62", "#bdb246", "#437c90")))
  names(colors) <- c("maxanc","maxanccolor")
  fillcolors <- as.character(unlist(lapply(as.character(unique(sort(mating_segments$sharedanc))),function(x) colors$maxanccolor[which(colors$maxanc==x)])))
  
  mat.tern.plot <- mating_segments
  mat.tern.plot$vector1 <- mat.tern.plot$prop1-mat.tern.plot$father1
  mat.tern.plot$vector1_meandif <- abs(mat.tern.plot$vector1 -mean(mat.tern.plot$vector1))
  mat.tern.plot$vector2 <- mat.tern.plot$prop2-mat.tern.plot$father2
  mat.tern.plot$vector2_meandif <- abs(mat.tern.plot$vector2 -mean(mat.tern.plot$vector2))
  mat.tern.plot$vector3 <- mat.tern.plot$prop3-mat.tern.plot$father3
  mat.tern.plot$vector3_meandif <- abs(mat.tern.plot$vector3 -mean(mat.tern.plot$vector3))
  mat.tern.plot$vectorlength <- sqrt((mat.tern.plot$vector1)^2+(mat.tern.plot$vector2)^2+(mat.tern.plot$vector3)^2)
  mat.tern.plot$vectorlength_meandif <- abs(mean(mat.tern.plot$vectorlength) - mat.tern.plot$vectorlength)
  highlighted <- sample(1:length(mat.tern.plot$sharedanc),10,replace=F)
  mat.tern.plot <- rbind(mat.tern.plot[-highlighted,],mat.tern.plot[highlighted,])
  rownames(mat.tern.plot) <- c(1:length(mat.tern.plot$index))
  #This allows to avoid coloured arrows from 1,0,0 to 0.5,0.5,0 at g=3 but could decolour also any to 0.5,0.25,0.25 and many others
  mat.tern.plot$sharedanc[((mat.tern.plot$prop1==0.5)&(mat.tern.plot$prop2%in%c(0.5,0.0)))]<-"Different ancestries"
  mat.tern.plot$sharedanc[((mat.tern.plot$prop1==0.0)&(mat.tern.plot$prop2%in%c(0.5,0.0,1)))]<-"Different ancestries"
  #add a 0,0 cocuple per ancestry
  mat.tern.plot.falsecouples<-as.data.frame(rbind(c(g,"false_couple_anc1","NA","1",1,"NA","NA","1",1,0,0,1,0,0,0,0,0,0,0,0,0,0),
  c(g,"false_couple_anc2","NA","2",1,"NA","NA","2",0,1,0,0,1,0,0,0,0,0,0,0,0,0),
  c(g,"falsee_couple_anc3","NA","3",1,"NA","NA","3",0,0,1,0,0,1,0,0,0,0,0,0,0,0),
  c(g,"falsee_couple_difanc","NA","3",1,"NA","NA","Different ancestries",0,0,1,0,0,1,0,0,0,0,0,0,0,0)))
  names(mat.tern.plot.falsecouples)<-names(mat.tern.plot)
  mat.tern.plot.falsecouples[,unlist(which(sapply(mat.tern.plot,class)=="numeric"))] <- sapply(mat.tern.plot.falsecouples[,unlist(which(sapply(mat.tern.plot,class)=="numeric"))],as.numeric)
  mat.tern.plot.falsecouples$sharedanc <- factor(mat.tern.plot.falsecouples$sharedanc,levels=c("Different ancestries","1","2","3"))   
  mat.tern.plot<-rbind(mat.tern.plot,mat.tern.plot.falsecouples)
  ##
  highlighted <- unlist(sapply(c(c(1:3),"Different ancestries"),function(anc_i)sample(which(mat.tern.plot$sharedanc==anc_i),min(4,sum(mat.tern.plot$sharedanc==anc_i)))))
  alpha_values <- rep(0.3,length(mat.tern.plot$sharedanc))
  alpha_values[highlighted] <- 1
  alpha_values[c((length(alpha_values)-3):length(alpha_values))]<-0
  size_values <- rep(0.8,length(mat.tern.plot$sharedanc))
  size_values[highlighted] <- 1.2
  size_values[c((length(size_values)-3):length(size_values))]<-0
  
  mating.ternary = ggtern(data=mat.tern.plot, 
                          aes(x=prop1,
                              y=prop2,
                              z=prop3,                     
                              xend=father1,
                              yend=father2,
                              zend=father3,
                              color=sharedanc,
                              fill=sharedanc)) +
    geom_segment(linetype=1,
                 alpha=alpha_values,
                 size=size_values,
                 arrow = arrow(length = unit(0.3,"cm"))
    ) +
    scale_color_manual("Ancestry",values=colors$maxanccolor)+
    scale_fill_manual("Ancestry",values=colors$maxanccolor)+
    Larrowlab("More ancestry 1")+
    Tarrowlab("More ancestry 2")+
    Rarrowlab("More ancestry 3")+
    tern_limit(T = 1.1, L = 1.1, R = 1.1)+
    theme_showarrows() +
    theme(legend.position = "none",
          #plot.margin = margin(10, 10, 10, 10, "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          #tern.panel.grid.major = element_line(linetype="solid",color="white"),
          tern.axis.title.show=F,
          tern.axis.arrow.text.T= element_text(color="#bdb246"),
          tern.axis.arrow.T=element_line(color="#bdb246"),
          tern.axis.line.T=element_line(color="#bdb246"),
          tern.axis.ticks.major.T=element_line(color="#bdb246"),
          tern.axis.text.T= element_text(color="#bdb246"),
          tern.axis.arrow.text.L= element_text(color="#b84a62"),
          tern.axis.arrow.L=element_line(color="#b84a62"),
          tern.axis.line.L=element_line(color="#b84a62"),
          tern.axis.ticks.major.L=element_line(color="#b84a62"),
          tern.axis.text.L=element_text(color="#b84a62"),
          tern.axis.arrow.text.R= element_text(color="#437c90"),
          tern.axis.arrow.R=element_line(color="#437c90"),
          tern.axis.line.R=element_line(color="#437c90"),
          tern.axis.ticks.major.R=element_line(color="#437c90"),
          tern.axis.text.R= element_text(color="#437c90"),
          plot.title = element_text(hjust = 0.2,vjust=-10))
    #ggtitle(paste("Mating at generation",gen))
  plottern <- mating.ternary
  return(plottern)
}

file_anc<-paste0("Ancestries_AutX_log_pulses_false_newgenmap_predicted_PUR_AM1_18.33351_AM2_0.34484_AM3_18.29093_SB1_0.09201_SB2_-0.28452_GEN_19_POP1_138_POP2_150_POP3_712_scale_1000_run_",run,".txt")

for (gen in c(2:19)){
  ggsave(paste(path,"training_tern_plot_PUR_gen_",sprintf("%03d", gen),"_run_",sprintf("%03d", run),".jpg",sep=""),device="jpg",width = 7,height=7,units="in",plot=plotsmatingternfunction(run,gen))
  ggsave(paste(path,"training_tern_plot_PUR_gen_",sprintf("%03d", gen),"_run_",sprintf("%03d", run),".pdf",sep=""),device="pdf",width = 7,height=7,units="in",plot=plotsmatingternfunction(run,gen))
  ggsave(paste(path,"training_couples_plot_PUR_gen_",sprintf("%03d", gen),"_run_",sprintf("%03d", run),".jpg",sep=""),device="jpg",width = 7,height=7,units="in",plot=plotsmatingbarsfunction(run,gen))
  ggsave(paste(path,"training_couples_plot_PUR_gen_",sprintf("%03d", gen),"_run_",sprintf("%03d", run),".pdf",sep=""),device="pdf",width = 7,height=7,units="in",plot=plotsmatingbarsfunction(run,gen))
}
ggsave(paste(path,"training_matingprobs_plot_PUR_gen_",sprintf("%03d", gen),"_run_",sprintf("%03d", run),".pdf",sep=""),device="pdf",width = 7,height=7,units="in",plot=plot_mating_prob)








