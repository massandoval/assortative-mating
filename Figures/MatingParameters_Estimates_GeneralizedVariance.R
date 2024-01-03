library(tidyr)
library(dplyr)
library(stats)
library(ggplot2)
library(truncnorm)
path<-"./"
calc_param <- function(x,y) {
  det=var(x)*var(y)-(cov(x,y))^2
  return(det)
}
all_pred_false<-c()
all_mse_false<-c()
for (pop in c("CLM","MXL","PEL","PUR","ASW","ACB")) {
  for (bins in c("21","22")){
    for (option in c("0","2")){
      for (LA in c("gnomix","rfmix")){  
        pred<-read.table(paste0(path,"predicted_pulses_FALSE_",pop,"_",LA,"_maxanc_bins_",bins,"_normoption_",option,"_onlyfemales.txt"))
        pred<-pred[,c(1,3,4)]
        names(pred)<-c("param","value","NN")
        pred$population<-pop
        pred$LA<-LA
        pred$bins<-bins
        pred$option<-option
        all_pred_false<-rbind(all_pred_false,pred)
      }
      mse<-read.table(paste0(path,"MSE_R2_test_pulses_FALSE_",pop,"_maxanc_bins_",bins,"_normoption_",option,"_onlyfemales.txt"))
      names(mse)<-c("param","opt","value","NN")
      mse <- mse %>% filter(opt=="MSE") %>% select(-opt)
      mse$population<-pop
      mse$bins<-bins
      mse$option<-option
      all_mse_false<-rbind(all_mse_false,mse)
    }
  }
}



all_pred_false <- all_pred_false %>% group_by(population, NN, bins, option) %>% filter(n() == 10) %>% 
  group_by(population, bins, option, LA) %>% top_n(1000, NN) %>%ungroup()
all_mse_false <- all_mse_false %>% group_by(population, NN, bins, option) %>% filter(n() == 5) %>% 
  group_by(population, bins, option) %>% top_n(1000, NN) %>%ungroup()

value_rfmix_column <- all_pred_false %>% filter (LA=="rfmix")  %>% rename(value_rfmix=value)%>% select(value_rfmix)
all_pred_false_plot <- all_pred_false %>% filter (LA=="gnomix") %>% select(-LA) %>% rename(value_gnomix=value) %>% mutate(value_rfmix=value_rfmix_column$value_rfmix) %>% mutate(ancestry=gsub(".*([0-9]$)", "\\1", param))
all_pred_false_plot<- all_pred_false_plot %>% mutate(param=factor(param),population=factor(population),bins=factor(bins),option=factor(option),ancestry=factor(ancestry)) %>%group_by(bins,option,population,param)
param_GV_df<-all_pred_false_plot %>%group_by(bins,option,population,param) %>% summarize(param_GV_value = calc_param(value_rfmix, value_gnomix)) 
all_pred_false_plot<-all_pred_false_plot %>% mutate(value_gnomix=ifelse(param%in%c("SB1","SB2"),value_gnomix*2-1,value_gnomix),
                                                    value_rfmix=ifelse(param%in%c("SB1","SB2"),value_rfmix*2-1,value_rfmix)) %>%
  mutate( value_gnomix=ifelse(param%in%c("AM1","AM2","AM3"),1-value_gnomix,value_gnomix),
          value_rfmix=ifelse(param%in%c("AM1","AM2","AM3"),1-value_rfmix,value_rfmix)) %>%
  group_by(bins,option,population,NN)
all_pred_false_plot <- merge(
  all_pred_false_plot%>%ungroup()%>%select(-value_gnomix,-ancestry)%>%pivot_wider(names_from=param,values_from=value_rfmix) %>% mutate(SB3=0-(SB1+SB2)) %>% pivot_longer(cols = c("AM1","AM2","AM3","SB1","SB2","SB3"),names_to="param",values_to = "value_rfmix"),
  all_pred_false_plot%>%ungroup()%>%select(-value_rfmix,-ancestry)%>%pivot_wider(names_from=param,values_from=value_gnomix) %>% mutate(SB3=0-(SB1+SB2)) %>% pivot_longer(cols = c("AM1","AM2","AM3","SB1","SB2","SB3"),names_to="param",values_to = "value_gnomix"),
  by=c("param","NN","population","bins","option")) %>% mutate(ancestry=gsub(".*([0-9]$)", "\\1", param)) 

param_GV_mse_df<- merge(
  param_GV_df %>% select("param","population","bins","option","param_GV_value"),
  all_mse_false %>% group_by(param,population,bins,option) %>% summarise(mean_mse=mean(value)) %>% select("param","population","bins","option","mean_mse"),
  by=c("param","population","bins","option")) %>% 
  mutate(ancestry=gsub(".*([0-9]$)", "\\1", param),option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
         bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows"))

param_GV_mse_df %>% mutate(standard_GV=(param_GV_value-min(param_GV_value))/(max(param_GV_value)-min(param_GV_value)),standard_mse=(mean_mse - min(mean_mse))/(max(mean_mse)-min(mean_mse))) %>%  group_by(bins,option)%>% summarise(mean_param_GV=mean(param_GV_value),mean_mean_mse=mean(mean_mse),mean_standard_GV=mean(standard_GV),mean_standard_mse=mean(standard_mse)) %>% mutate(prod_standard_GV_mse=mean_standard_GV*mean_standard_mse) 
table_GV_mse_false<-param_GV_mse_df %>% mutate(standard_GV=(param_GV_value-min(param_GV_value))/(max(param_GV_value)-min(param_GV_value)),standard_mse=(mean_mse - min(mean_mse))/(max(mean_mse)-min(mean_mse))) %>%  group_by(bins,option)%>% summarise(mean_GV=mean(param_GV_value),mean_MSE=mean(mean_mse))#,prod_GV_MSE=mean(standard_GV)*mean(standard_mse))
#table_GV_mse_true<-param_GV_mse_df %>% group_by(bins,option)  %>% summarise(mean_GV=mean(param_GV_value),mean_mse=mean(mean_mse)) 
table_GV_mse_false<-table_GV_mse_false %>% mutate(Windows=gsub(pattern="\n",replacement=" ",x=bins),
                                                Scaled=gsub(pattern="\n",replacement="",x=option),
                                                mean_GV=formatC(mean_GV,format="e",digits=2),
                                                mean_MSE=as.character(round(mean_MSE,4))) %>% #,
                                                #prod_GV_MSE=as.character(round(prod_GV_MSE,4))) %>%
  ungroup() %>% select(-bins,-option) %>% relocate(c("mean_GV","mean_MSE"), .after = last_col())


write.table(table_GV_mse_false,file=paste0(path,"Table_07_GnomixVsRFMix_pulsesFALSE_table_GV_mse_Reviews_oct2023.txt"),row.names=FALSE,quote=FALSE,col.names=TRUE,sep=",")


plot_FALSE_AM_all_options<- ggplot(all_pred_false_plot %>% mutate(option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
                                                                  bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows"))%>%
                                     filter(param%in%c("AM1","AM2","AM3")),
                                   aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.5)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(population+param~bins+option)+
  geom_text(data = param_GV_mse_df%>% filter(param%in%c("AM1","AM2","AM3")), aes(x = 0.025, y = 0.975, label = paste0("GV=",formatC(param_GV_value,format="e",digits=2),"\nMSE=",round(mean_mse,4))),
            hjust=0, vjust = 1, size = 2)+
  xlim(0,1)+
  ylim(0,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the One Pulse model\nfrom continuous ancestry tract lengths profile\nwith Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5,angle=90))

ggsave(filename=paste0(path,"SuppFigure_MatingProbs_03_GnomixVsRFMix_pulsesFALSE_AM_alloptions_Reviews_oct2023.pdf"),plot=plot_FALSE_AM_all_options,height=25.0,width=6.0,units="in")

plot_FALSE_SB_all_options<- ggplot(all_pred_false_plot %>% mutate(option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
                                                                  bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows")) %>%
                                     filter(param%in%c("SB1","SB2","SB3")),
                                   aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_hline(yintercept = 0,size=0.5,color="gray80")+
  geom_vline(xintercept = 0,size=0.5,color="gray80")+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.5)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(population+param~bins+option)+
  geom_text(data = param_GV_mse_df%>% filter(param%in%c("SB1","SB2","SB3")), aes(x = -0.95, y = 0.95, label = paste0("GV=",formatC(param_GV_value,format="e",digits=2),"\nMSE=",round(mean_mse,4))),
            hjust=0, vjust = 1, size = 2)+
  xlim(-1,1)+
  ylim(-1,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the One Pulse model\nfrom continuous ancestry tract lengths profile\nwith Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5,angle=90))

ggsave(filename=paste0(path,"SuppFigure_MatingProbs_04_GnomixVsRFMix_pulsesFALSE_SB_alloptions_Reviews_oct2023.pdf"),plot=plot_FALSE_SB_all_options,height=25.0,width=6.0,units="in")

plot_FALSE_AM_bins_22_option_0<- ggplot(all_pred_false_plot %>% filter(param%in%c("AM1","AM2","AM3")&(option=="0")&(bins=="22")),aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.2,size=2,stroke=0)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(param~population)+
  #geom_text(data = param_GV_mse_df%>% filter(param%in%c("AM1","AM2","AM3")&(option=="0")&(bins=="22")), aes(x = 0.3, y = 0.9, label = paste("GV=",round(param_GV_value, 5))), vjust = 1, size = 3)+
  xlim(0,1)+
  ylim(0,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the One Pulse model\nfrom tract length profile with Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 8,angle=90))

ggsave(filename=paste0(path,"GnomixVsRFMix_pulsesFALSE_AM_bins_22_option_0_onlyfemales.pdf"),plot=plot_FALSE_AM_bins_22_option_0,height=5.0,width=9.0,units="in")

plot_FALSE_SB_bins_22_option_0<- ggplot(all_pred_false_plot %>% filter(param%in%c("SB1","SB2","SB3")&(option=="0")&(bins=="22")),aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_hline(yintercept = 0,size=0.5,color="gray80")+
  geom_vline(xintercept = 0,size=0.5,color="gray80")+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.2,size=2,stroke=0)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(param~population)+
  #geom_text(data = param_GV_mse_df%>% filter(param%in%c("SB1","SB2","SB3")&(option=="0")&(bins=="22")), aes(x = -0.7, y = 0.9, label = paste("GV=",round(param_GV_value, 5))),vjust = 1, size = 3)+
  xlim(-1,1)+
  ylim(-1,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the One Pulse model\nfrom tract length profile with Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 8,angle=90))

ggsave(filename=paste0(path,"GnomixVsRFMix_pulsesFALSE_SB_bins_22_option_0_onlyfemales.pdf"),plot=plot_FALSE_SB_bins_22_option_0,height=5.0,width=9.0,units="in")

#View(all_pred_false_plot %>%group_by(bins,option,population,param) %>% summarize(param_GV_value = calc_param(value_rfmix, value_gnomix))%>%group_by(option,bins)  %>% summarize(sum_param=sum(param_GV_value)) %>% arrange(sum_param))


###TRUE

all_pred_true<-c()
all_mse_true<-c()
for (pop in c("CLM","MXL","PEL","PUR","ASW","ACB")) {
  for (bins in c("21","22")){
    for (option in c("0","2")){
      for (LA in c("gnomix","rfmix")){
        pred<-read.table(paste0(path,"predicted_pulses_TRUE_",pop,"_",LA,"_maxanc_bins_",bins,"_normoption_",option,"_onlyfemales.txt"))
        pred<-pred[,c(1,3,4)]
        names(pred)<-c("param","value","NN")
        pred$population<-pop
        pred$LA<-LA
        pred$bins<-bins
        pred$option<-option
        all_pred_true<-rbind(all_pred_true,pred)
      }
      mse<-read.table(paste0(path,"MSE_R2_test_pulses_TRUE_",pop,"_gnomix_maxanc_bins_",bins,"_normoption_",option,"_onlyfemales.txt"))
      names(mse)<-c("param","opt","value","NN")
      mse <- mse %>% filter(opt=="MSE") %>% select(-opt)
      mse$population<-pop
      mse$bins<-bins
      mse$option<-option
      all_mse_true<-rbind(all_mse_true,mse)
    }
  }
}

all_pred_true <- all_pred_true %>% group_by(population, NN, bins, option) %>% filter(n() == 16) %>% 
group_by(population, bins, option, LA) %>% top_n(1000, NN) %>%ungroup()
all_mse_true <- all_mse_true %>% group_by(population, NN, bins, option) %>% filter(n() == 8) %>% 
group_by(population, bins, option) %>% top_n(1000, NN) %>%ungroup()

a<-all_pred_true %>% filter (LA=="rfmix") %>% rename(value_rfmix=value)
b<-all_pred_true %>% filter (LA=="gnomix") %>% rename(value_gnomix=value)
all_pred_true_plot <-merge(a,b,by=c("param","NN","population","bins","option")) %>% mutate(ancestry=gsub(".*([0-9]$)", "\\1", param)) %>% select(-LA.x,-LA.y)
all_pred_true_plot<- all_pred_true_plot %>% mutate(param=factor(param),population=factor(population),bins=factor(bins),option=factor(option),ancestry=factor(ancestry)) %>%group_by(bins,option,population,param)
param_GV_df<-all_pred_true_plot  %>% mutate(param=recode_factor(param,prop_pop1_t1="GFR1",prop_pop2_t1="GFR2",prop_pop3_t1="GFR3")) %>%group_by(bins,option,population,param)%>% summarize(ancestry=ancestry, param_GV_value = calc_param(value_rfmix, value_gnomix))

#sb from -1 to 1, am from 0 to 1
all_pred_true_plot<-all_pred_true_plot %>% mutate(value_gnomix=ifelse(param%in%c("SB1","SB2"),value_gnomix*2-1,value_gnomix),
                                                    value_rfmix=ifelse(param%in%c("SB1","SB2"),value_rfmix*2-1,value_rfmix)) %>%
  mutate( value_gnomix=ifelse(param%in%c("AM1","AM2","AM3"),1-value_gnomix,value_gnomix),
          value_rfmix=ifelse(param%in%c("AM1","AM2","AM3"),1-value_rfmix,value_rfmix)) %>%
  group_by(bins,option,population,NN) 
all_pred_true_plot <- merge(
  all_pred_true_plot%>%select(-value_gnomix,-ancestry)%>%pivot_wider(names_from=param,values_from=value_rfmix) %>% mutate(SB3=0-(SB1+SB2)) %>% rename(GFR1=prop_pop1_t1,GFR2=prop_pop2_t1,GFR3=prop_pop3_t1) %>% pivot_longer(cols = c("AM1","AM2","AM3","SB1","SB2","SB3","GFR1","GFR2","GFR3"),names_to="param",values_to = "value_rfmix"),
  all_pred_true_plot%>%select(-value_rfmix,-ancestry)%>%pivot_wider(names_from=param,values_from=value_gnomix) %>% mutate(SB3=0-(SB1+SB2)) %>% rename(GFR1=prop_pop1_t1,GFR2=prop_pop2_t1,GFR3=prop_pop3_t1) %>% pivot_longer(cols = c("AM1","AM2","AM3","SB1","SB2","SB3","GFR1","GFR2","GFR3"),names_to="param",values_to = "value_gnomix"),
  by=c("param","NN","population","bins","option")) %>% mutate(ancestry=gsub(".*([0-9]$)", "\\1", param)) %>%
  mutate(value_rfmix = ifelse(param %in% c("GFR1","GFR2","GFR3"), 1-value_rfmix,value_rfmix),
         value_gnomix = ifelse(param %in% c("GFR1","GFR2","GFR3"), 1-value_gnomix,value_gnomix))

param_GV_mse_df<- merge(
  param_GV_df %>% select("param","population","bins","option","param_GV_value"),
  all_mse_true %>% mutate(param=recode_factor(param,prop_pop1_t1="GFR1",prop_pop2_t1="GFR2",prop_pop3_t1="GFR3"))%>% group_by(param,population,bins,option) %>% summarise(mean_mse=mean(value)) %>% select("param","population","bins","option","mean_mse"),
  by=c("param","population","bins","option")) %>% 
  mutate(ancestry=gsub(".*([0-9]$)", "\\1", param),option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
         bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows"))

table_GV_mse_true<-param_GV_mse_df %>% mutate(standard_GV=(param_GV_value-min(param_GV_value))/(max(param_GV_value)-min(param_GV_value)),standard_mse=(mean_mse - min(mean_mse))/(max(mean_mse)-min(mean_mse))) %>%  group_by(bins,option)%>% summarise(mean_GV=mean(param_GV_value),mean_MSE=mean(mean_mse))#,prod_GV_MSE=mean(standard_GV)*mean(standard_mse))
table_GV_mse_true<-table_GV_mse_true %>% mutate(Windows=gsub(pattern="\n",replacement=" ",x=bins),
                                                Scaled=gsub(pattern="\n",replacement="",x=option),
                                                mean_GV=formatC(mean_GV,format="e",digits=2),
                                                mean_MSE=as.character(round(mean_MSE,4))) %>%
                                                #prod_GV_MSE=as.character(round(prod_GV_MSE,4))) %>%
  ungroup() %>% select(-bins,-option) %>% relocate(c("mean_GV","mean_MSE"), .after = last_col())
write.table(table_GV_mse_true,file=paste0(path,"Table_08_GnomixVsRFMix_pulsesTRUE_table_GV_mse_Reviews_oct2023.txt"),row.names=FALSE,quote=FALSE,col.names=TRUE,sep=",")

plot_TRUE_AM_all_options<- ggplot(all_pred_true_plot %>% mutate(option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
                                                                bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows")) %>%
                                    filter(param%in%c("AM1","AM2","AM3")),
                                  aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.5)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(population+param~bins+option)+
  geom_text(data = param_GV_mse_df%>% filter(param%in%c("AM1","AM2","AM3")), aes(x = 0.025, y = 0.975, label = paste0("GV=",formatC(param_GV_value,format="e",digits=2),"\nMSE=",round(mean_mse,4))),
            hjust=0, vjust = 1, size = 2)+
  xlim(0,1)+
  ylim(0,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the Two Pulses model\nfrom tract length profile with Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5,angle=90))

ggsave(filename=paste0(path,"SuppFigure_MatingProbs_05_GnomixVsRFMix_pulsesTRUE_AM_alloptions_Reviews_oct2023.pdf"),plot=plot_TRUE_AM_all_options,height=25.0,width=6.0,units="in")

plot_TRUE_SB_all_options<- ggplot(all_pred_true_plot %>% mutate(option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
                                                                bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows")) %>%
                                    filter(param%in%c("SB1","SB2","SB3")),
                                  aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_hline(yintercept = 0,size=0.5,color="gray80")+
  geom_vline(xintercept = 0,size=0.5,color="gray80")+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.5)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(population+param~bins+option)+
  geom_text(data = param_GV_mse_df%>% filter(param%in%c("SB1","SB2","SB3")), aes(x = -0.95, y = 0.95, label = paste0("GV=",formatC(param_GV_value,format="e",digits=2),"\nMSE=",round(mean_mse,4))),
            hjust=0, vjust = 1, size = 2)+
  xlim(-1,1)+
  ylim(-1,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the One Pulse model\nfrom continuous ancestry tract lengths profile\nwith Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5,angle=90))

plot_TRUE_SB_all_options
ggsave(filename=paste0(path,"SuppFigure_MatingProbs_06_GnomixVsRFMix_pulsesTRUE_SB_alloptions_Reviews_oct2023.pdf"),plot=plot_TRUE_SB_all_options,height=25.0,width=6.0,units="in")

plot_TRUE_GFR_all_options<- ggplot(all_pred_true_plot %>% mutate(option=ifelse(option=="0","Divided by total\n sum of tracts","Raw"),
                                                                 bins=ifelse(bins=="21","Without shortest tract\nwindow (<0.2cM)","All windows")) %>%
                                     filter(param%in%c("GFR1","GFR2","GFR3")),
                                   aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_hline(yintercept = 0,size=0.5,color="gray80")+
  geom_vline(xintercept = 0,size=0.5,color="gray80")+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.5)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(population+param~bins+option)+
  geom_text(data = param_GV_mse_df%>% filter(param%in%c("GFR1","GFR2","GFR3")), aes(x = 0.025, y = 0.05, label = paste0("GV=",formatC(param_GV_value,format="e",digits=2),"\nMSE=",round(mean_mse,4))),
            hjust=0, vjust = 0, size = 2)+
  xlim(0,1)+
  ylim(0,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the One Pulse model\nfrom continuous ancestry tract lengths profile\nwith Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 5,angle=90))

ggsave(filename=paste0(path,"SuppFigure_MatingProbs_07_GnomixVsRFMix_pulsesTRUE_GFR_alloptions_Reviews_oct2023.pdf"),plot=plot_TRUE_GFR_all_options,height=25.0,width=6.0,units="in")



plot_TRUE_AM_bins_22_option_0<- ggplot(all_pred_true_plot %>% filter(param%in%c("AM1","AM2","AM3")&(option=="0")&(bins=="22")),aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.2,size=2,stroke=0)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(param~population)+
  #geom_text(data = param_GV_mse_df%>% filter(param%in%c("AM1","AM2","AM3")&(option=="0")&(bins=="22")), aes(x = 0.3, y = 0.9, label = paste("GV=",round(param_GV_value, 5))), vjust = 1, size = 3)+
  xlim(0,1)+
  ylim(0,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the Two Pulses model\nfrom tract length profile with Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 8,angle=90))

ggsave(filename=paste0(path,"GnomixVsRFMix_pulsesTRUE_AM_bins_22_option_0_onlyfemales.pdf"),plot=plot_TRUE_AM_bins_22_option_0,height=5.0,width=9.0,units="in")




plot_TRUE_SB_bins_22_option_0<- ggplot(all_pred_true_plot %>% filter(param%in%c("SB1","SB2","SB3")&(option=="0")&(bins=="22")),aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_hline(yintercept = 0,size=0.5,color="gray80")+
  geom_vline(xintercept = 0,size=0.5,color="gray80")+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.2,size=2,stroke=0)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(param~population)+
  #geom_text(data = param_GV_mse_df%>% filter(param%in%c("SB1","SB2","SB3")&(option=="0")&(bins=="22")), aes(x = -0.7, y = 0.9, label = paste("GV=",round(param_GV_value, 5))),vjust = 1, size = 3)+
  xlim(-1,1)+
  ylim(-1,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the Two Pulses model\nfrom tract length profile with Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 8,angle=90))

ggsave(filename=paste0(path,"GnomixVsRFMix_pulsesTRUE_SB_bins_22_option_0_onlyfemales.pdf"),plot=plot_TRUE_SB_bins_22_option_0,height=5.0,width=9.0,units="in")

plot_TRUE_GFR_bins_22_option_0<- ggplot(all_pred_true_plot %>% filter(param%in%c("GFR1","GFR2","GFR3")&(option=="0")&(bins=="22")),aes(x=value_rfmix,y=value_gnomix,color=ancestry,fill=ancestry))+
  geom_abline(intercept = 0, slope = 1,size=0.2,linetype="dashed")+
  geom_point(alpha=0.2,size=2,stroke=0)+
  #stat_cor(size=2,label.x=0,label.y=0.9)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(param~population)+
  #geom_text(data = param_GV_mse_df%>% filter(param%in%c("GFR1","GFR2","GFR3")&(option=="0")&(bins=="22")), aes(x = -0.7, y = 0.9, label = paste("GV=",round(param_GV_value, 5))),vjust = 1, size = 3)+
  xlim(-1,1)+
  ylim(-1,1)+
  ylab("Gnomix")+
  xlab("RFMix")+
  ggtitle(("Prediction of mating parameters in the Two Pulses model\nfrom tract length profile with Gnomix or RFMix for the local ancestry inference "))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_text(size = 8,angle=90))

ggsave(filename=paste0(path,"GnomixVsRFMix_pulsesTRUE_GFR_bins_22_option_0_onlyfemales.pdf"),plot=plot_TRUE_GFR_bins_22_option_0,height=5.0,width=9.0,units="in")




######

quantile05<-function(x) {
  quantile(x,0.05)
}
quantile95<-function(x) {
  quantile(x,0.95)
}



full_allpops_summarised_false <-c()
for (pop in c("CLM","MXL","PEL","PUR","ASW","ACB")) {
  full_pop_summarised <- all_pred_false_plot %>%
    filter(population==pop,bins==22,option==0) %>%
  select(-population, -ancestry, -bins, -option) %>% rename(rfmix=value_rfmix,gnomix=value_gnomix)%>%
  pivot_wider (names_from = param, values_from = c(rfmix,gnomix)) %>% select(-NN) %>%
  summarize_all(.funs = list(
      mean = ~ mean(.),
      low_CI95 = ~ quantile(.x, 0.05),
      high_CI95 = ~ quantile(.x, 0.95))) %>% mutate(population=pop)
  full_allpops_summarised_false <-  bind_rows(full_allpops_summarised_false,full_pop_summarised)
}

full_allpops_summarised_true <-c()
for (pop in c("CLM","MXL","PEL","PUR","ASW","ACB")) {
  full_pop_summarised <- all_pred_true_plot %>%
    filter(population==pop,bins==22,option==0) %>%
    select(-population, -ancestry, -bins, -option) %>% rename(rfmix=value_rfmix,gnomix=value_gnomix)%>%
    pivot_wider (names_from = param, values_from = c(rfmix,gnomix)) %>% select(-NN) %>%
    summarize_all(.funs = list(
      mean = ~ mean(.),
      low_CI95 = ~ quantile(.x, 0.05),
      high_CI95 = ~ quantile(.x, 0.95))) %>% mutate(population=pop)
  #full_pop_summarised <- full_pop %>% select(-V4,-population,-pulses)%>% summarise(across(everything(), list(mean = mean))) %>% mutate(population=pop, pulses=pulses)
  full_allpops_summarised_true <-  bind_rows(full_allpops_summarised_true,full_pop_summarised)
}




AM1<-paste0(round(full_allpops_summarised_false$rfmix_AM1_mean,2),"(",round(full_allpops_summarised_false$rfmix_AM1_low_CI95,2),",",round(full_allpops_summarised_false$rfmix_AM1_high_CI95,2),")")
AM2<-paste0(round(full_allpops_summarised_false$rfmix_AM2_mean,2),"(",round(full_allpops_summarised_false$rfmix_AM2_low_CI95,2),",",round(full_allpops_summarised_false$rfmix_AM2_high_CI95,2),")")
AM3<-paste0(round(full_allpops_summarised_false$rfmix_AM3_mean,2),"(",round(full_allpops_summarised_false$rfmix_AM3_low_CI95,2),",",round(full_allpops_summarised_false$rfmix_AM3_high_CI95,2),")")
SB1<-paste0(round(full_allpops_summarised_false$rfmix_SB1_mean,2),"(",round(full_allpops_summarised_false$rfmix_SB1_low_CI95,2),",",round(full_allpops_summarised_false$rfmix_SB1_high_CI95,2),")")
SB2<-paste0(round(full_allpops_summarised_false$rfmix_SB2_mean,2),"(",round(full_allpops_summarised_false$rfmix_SB2_low_CI95,2),",",round(full_allpops_summarised_false$rfmix_SB2_high_CI95,2),")")
SB3<-paste0(round(full_allpops_summarised_false$rfmix_SB3_mean,2),"(",round(full_allpops_summarised_false$rfmix_SB3_low_CI95,2),",",round(full_allpops_summarised_false$rfmix_SB3_high_CI95,2),")")
population<-c("CLM","MXL","PEL","PUR","ASW","ACB")
estimates_pulses_false_rfmix<-rbind(population,AM1,AM2,AM3,SB1,SB2,SB3)
write.table(estimates_pulses_false_rfmix, file=paste0(path,"Table_03_Mating_parameters_estimates_pulses_false_rfmix_Reviews_oct2023.txt"), row.names=TRUE, col.names=FALSE,quote=FALSE,sep=";")

AM1<-paste0(round(full_allpops_summarised_false$gnomix_AM1_mean,2),"(",round(full_allpops_summarised_false$gnomix_AM1_low_CI95,2),",",round(full_allpops_summarised_false$gnomix_AM1_high_CI95,2),")")
AM2<-paste0(round(full_allpops_summarised_false$gnomix_AM2_mean,2),"(",round(full_allpops_summarised_false$gnomix_AM2_low_CI95,2),",",round(full_allpops_summarised_false$gnomix_AM2_high_CI95,2),")")
AM3<-paste0(round(full_allpops_summarised_false$gnomix_AM3_mean,2),"(",round(full_allpops_summarised_false$gnomix_AM3_low_CI95,2),",",round(full_allpops_summarised_false$gnomix_AM3_high_CI95,2),")")
SB1<-paste0(round(full_allpops_summarised_false$gnomix_SB1_mean,2),"(",round(full_allpops_summarised_false$gnomix_SB1_low_CI95,2),",",round(full_allpops_summarised_false$gnomix_SB1_high_CI95,2),")")
SB2<-paste0(round(full_allpops_summarised_false$gnomix_SB2_mean,2),"(",round(full_allpops_summarised_false$gnomix_SB2_low_CI95,2),",",round(full_allpops_summarised_false$gnomix_SB2_high_CI95,2),")")
SB3<-paste0(round(full_allpops_summarised_false$gnomix_SB3_mean,2),"(",round(full_allpops_summarised_false$gnomix_SB3_low_CI95,2),",",round(full_allpops_summarised_false$gnomix_SB3_high_CI95,2),")")
population<-c("CLM","MXL","PEL","PUR","ASW","ACB")
estimates_pulses_false_gnomix<-rbind(population,AM1,AM2,AM3,SB1,SB2,SB3)
write.table(estimates_pulses_false_gnomix, file=paste0(path,"Table_04_Mating_parameters_estimates_pulses_false_gnomix_Reviews_oct2023.txt"), row.names=TRUE, col.names=FALSE,quote=FALSE,sep=";")




AM1<-paste0(round(full_allpops_summarised_true$rfmix_AM1_mean,2),"(",round(full_allpops_summarised_true$rfmix_AM1_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_AM1_high_CI95,2),")")
AM2<-paste0(round(full_allpops_summarised_true$rfmix_AM2_mean,2),"(",round(full_allpops_summarised_true$rfmix_AM2_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_AM2_high_CI95,2),")")
AM3<-paste0(round(full_allpops_summarised_true$rfmix_AM3_mean,2),"(",round(full_allpops_summarised_true$rfmix_AM3_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_AM3_high_CI95,2),")")
SB1<-paste0(round(full_allpops_summarised_true$rfmix_SB1_mean,2),"(",round(full_allpops_summarised_true$rfmix_SB1_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_SB1_high_CI95,2),")")
SB2<-paste0(round(full_allpops_summarised_true$rfmix_SB2_mean,2),"(",round(full_allpops_summarised_true$rfmix_SB2_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_SB2_high_CI95,2),")")
SB3<-paste0(round(full_allpops_summarised_true$rfmix_SB3_mean,2),"(",round(full_allpops_summarised_true$rfmix_SB3_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_SB3_high_CI95,2),")")
GFR1<-paste0(round(full_allpops_summarised_true$rfmix_GFR1_mean,2),"(",round(full_allpops_summarised_true$rfmix_GFR1_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_GFR1_high_CI95,2),")")
GFR2<-paste0(round(full_allpops_summarised_true$rfmix_GFR2_mean,2),"(",round(full_allpops_summarised_true$rfmix_GFR2_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_GFR2_high_CI95,2),")")
GFR3<-paste0(round(full_allpops_summarised_true$rfmix_GFR3_mean,2),"(",round(full_allpops_summarised_true$rfmix_GFR3_low_CI95,2),",",round(full_allpops_summarised_true$rfmix_GFR3_high_CI95,2),")")
population<-c("CLM","MXL","PEL","PUR","ASW","ACB")
estimates_pulses_true_rfmix<-rbind(population,AM1,AM2,AM3,SB1,SB2,SB3,GFR1,GFR2,GFR3)
write.table(estimates_pulses_true_rfmix, file=paste0(path,"Table_05_Mating_parameters_estimates_pulses_true_rfmix_Reviews_oct2023.txt"), row.names=TRUE, col.names=FALSE,quote=FALSE,sep=";")


AM1<-paste0(round(full_allpops_summarised_true$gnomix_AM1_mean,2),"(",round(full_allpops_summarised_true$gnomix_AM1_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_AM1_high_CI95,2),")")
AM2<-paste0(round(full_allpops_summarised_true$gnomix_AM2_mean,2),"(",round(full_allpops_summarised_true$gnomix_AM2_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_AM2_high_CI95,2),")")
AM3<-paste0(round(full_allpops_summarised_true$gnomix_AM3_mean,2),"(",round(full_allpops_summarised_true$gnomix_AM3_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_AM3_high_CI95,2),")")
SB1<-paste0(round(full_allpops_summarised_true$gnomix_SB1_mean,2),"(",round(full_allpops_summarised_true$gnomix_SB1_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_SB1_high_CI95,2),")")
SB2<-paste0(round(full_allpops_summarised_true$gnomix_SB2_mean,2),"(",round(full_allpops_summarised_true$gnomix_SB2_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_SB2_high_CI95,2),")")
SB3<-paste0(round(full_allpops_summarised_true$gnomix_SB3_mean,2),"(",round(full_allpops_summarised_true$gnomix_SB3_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_SB3_high_CI95,2),")")
GFR1<-paste0(round(full_allpops_summarised_true$gnomix_GFR1_mean,2),"(",round(full_allpops_summarised_true$gnomix_GFR1_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_GFR1_high_CI95,2),")")
GFR2<-paste0(round(full_allpops_summarised_true$gnomix_GFR2_mean,2),"(",round(full_allpops_summarised_true$gnomix_GFR2_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_GFR2_high_CI95,2),")")
GFR3<-paste0(round(full_allpops_summarised_true$gnomix_GFR3_mean,2),"(",round(full_allpops_summarised_true$gnomix_GFR3_low_CI95,2),",",round(full_allpops_summarised_true$gnomix_GFR3_high_CI95,2),")")
population<-c("CLM","MXL","PEL","PUR","ASW","ACB")
estimates_pulses_true_gnomix<-rbind(population,AM1,AM2,AM3,SB1,SB2,SB3,GFR1,GFR2,GFR3)
write.table(estimates_pulses_true_gnomix, file=paste0(path,"Table_06_Mating_parameters_estimates_pulses_true_gnomix_Reviews_oct2023.txt"), row.names=TRUE, col.names=FALSE,quote=FALSE,sep=";")
table_GV_mse_false
table_GV_mse_true



#####
##mating prob plots

predicted_summary_sd<-rbind(full_allpops_summarised_false[,c(grep("mean",names(full_allpops_summarised_false)),ncol(full_allpops_summarised_false))] %>% pivot_longer(cols=1:(ncol(.)-1),names_to="param",values_to="value") %>% separate(param,c("LA","param","mean"),"_") %>% select(-mean) %>% extract(param,c("parameter","anc"),"([A-Z].*)([0-9])") %>% mutate(pulses="false"),
full_allpops_summarised_true[,c(grep("mean",names(full_allpops_summarised_true)),ncol(full_allpops_summarised_true))] %>% pivot_longer(cols=1:(ncol(.)-1),names_to="param",values_to="value") %>% separate(param,c("LA","param","mean"),"_") %>% select(-mean) %>% extract(param,c("parameter","anc"),"([A-Z].*)([0-9])") %>% mutate(pulses="true"))

trunc_distribution_allpops<-matrix(nrow=0,ncol=9)

for (LA in c("gnomix","rfmix")){
for (pop in c("CLM","MXL","PEL","PUR","ACB","ASW")){
  for (pulses in c("false","true")){
    for (anc in c(1:3)){
      AM<-as.numeric(predicted_summary_sd$value[((predicted_summary_sd$parameter=="AM")&(predicted_summary_sd$anc==anc)&(predicted_summary_sd$population==pop)&(predicted_summary_sd$LA==LA)&(predicted_summary_sd$pulses==pulses))])
      SB<-as.numeric(predicted_summary_sd$value[((predicted_summary_sd$parameter=="SB")&(predicted_summary_sd$anc==anc)&(predicted_summary_sd$population==pop)&(predicted_summary_sd$LA==LA)&(predicted_summary_sd$pulses==pulses))])
      trunc_distribution<-cbind(seq(-1,1,0.001),dtruncnorm(seq(-1,1,0.001),-1,1,SB,3^((1-AM)*7-4)),0,AM,SB,anc,pulses,pop,LA)   
      #trunc_distribution<-cbind(seq(-1,1,0.001),dtruncnorm(seq(-1,1,0.001),-1,1,SB,1-AM),0,AM,SB,anc,pulses,pop,LA)   
      if (SB>=0){
        trunc_distribution[,3][(as.numeric(trunc_distribution[,1])>=0)&(as.numeric(trunc_distribution[,1])<SB)]<-"7"
      }else{ 
        trunc_distribution[,3][(as.numeric(trunc_distribution[,1])<0)&(as.numeric(trunc_distribution[,1])>SB)]<-"7"
      }
      
      trunc_distribution_allpops<-rbind(trunc_distribution_allpops,trunc_distribution)
      
}}}}
trunc_distribution_df<-as.data.frame(trunc_distribution_allpops)
names(trunc_distribution_df)[1]<-"x"
names(trunc_distribution_df)[2]<-"value"
names(trunc_distribution_df)[3]<-"SB_square"
names(trunc_distribution_df)[8]<-"population"
trunc_distribution_df$value<-as.numeric(as.character(trunc_distribution_df$value))
trunc_distribution_df$SB_square<-as.numeric(as.character(trunc_distribution_df$SB_square))
trunc_distribution_df$x<-as.numeric(as.character(trunc_distribution_df$x))
trunc_distribution_df$AM<-as.numeric(as.character(trunc_distribution_df$AM))
trunc_distribution_df$SB<-as.numeric(as.character(trunc_distribution_df$SB))
trunc_distribution_df$pulses[trunc_distribution_df$pulses=="false"]<-"One pulse 19 generations ago"
trunc_distribution_df$pulses[trunc_distribution_df$pulses=="true"]<-"Two pulses 19 and 9 generations ago"


v.line.data<-predicted_summary_sd[predicted_summary_sd$parameter=="SB",]
names(v.line.data)[5]<-"SB"
v.line.data$pulses[v.line.data$pulses=="false"]<-"One pulse 19 generations ago"
v.line.data$pulses[v.line.data$pulses=="true"]<-"Two pulses 19 and 9 generations ago"
v.line.data$pulses<-as.factor(v.line.data$pulses)
v.line.data$anc<-as.factor(v.line.data$anc)
v.line.data$AMvalue <- predicted_summary_sd$value[predicted_summary_sd$parameter=="AM"]
v.line.data$AM<-3^((1-predicted_summary_sd$value[predicted_summary_sd$parameter=="AM"])*7-3)
v.line.data<-v.line.data[,!names(v.line.data)=="parameter"]

text_true_rfmix<-merge(as.data.frame(t(estimates_pulses_true_rfmix)) %>% select(-SB1,-SB2,-SB3,-GFR1,-GFR2,-GFR3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="AMtext") %>% extract(param,"anc",".*([0-9])"),
                       as.data.frame(t(estimates_pulses_true_rfmix)) %>% select(-AM1,-AM2,-AM3,-GFR1,-GFR2,-GFR3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="SBtext") %>% extract(param,"anc",".*([0-9])")) %>% mutate(pulses="Two pulses 19 and 9 generations ago",LA="rfmix")
text_false_rfmix<-merge(as.data.frame(t(estimates_pulses_false_rfmix)) %>% select(-SB1,-SB2,-SB3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="AMtext") %>% extract(param,"anc",".*([0-9])"),
                        as.data.frame(t(estimates_pulses_false_rfmix)) %>% select(-AM1,-AM2,-AM3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="SBtext") %>% extract(param,"anc",".*([0-9])")) %>% mutate(pulses="One pulse 19 generations ago",LA="rfmix")
text_true_gnomix<-merge(as.data.frame(t(estimates_pulses_true_gnomix)) %>% select(-SB1,-SB2,-SB3,-GFR1,-GFR2,-GFR3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="AMtext") %>% extract(param,"anc",".*([0-9])"),
                        as.data.frame(t(estimates_pulses_true_gnomix)) %>% select(-AM1,-AM2,-AM3,-GFR1,-GFR2,-GFR3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="SBtext") %>% extract(param,"anc",".*([0-9])")) %>% mutate(pulses="Two pulses 19 and 9 generations ago",LA="gnomix")
text_false_gnomix<-merge(as.data.frame(t(estimates_pulses_false_gnomix)) %>% select(-SB1,-SB2,-SB3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="AMtext") %>% extract(param,"anc",".*([0-9])"),
                         as.data.frame(t(estimates_pulses_false_gnomix)) %>% select(-AM1,-AM2,-AM3) %>% pivot_longer(cols=c(2:4),names_to = "param",values_to="SBtext") %>% extract(param,"anc",".*([0-9])")) %>% mutate(pulses="One pulse 19 generations ago",LA="gnomix")
v.line.data<- merge(v.line.data,
                    rbind(text_true_rfmix,
                          text_false_rfmix,
                          text_true_gnomix,
                          text_false_gnomix))

plot_mating_prob<-ggplot()+
  geom_vline(xintercept=0,color="black",size=0.1)+
  geom_line(data=trunc_distribution_df,aes(x=x,y=value,color=anc),size=1.2)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  facet_grid(population~pulses+LA)+
  ylim(0,5)+
  ylab("Density")+
  xlab("Sex Bias (ancestry proportion of male - ancestry proportion of female)")+
  ggtitle(("Inferred mating probabilities"))+
  theme_minimal()+
  theme(legend.position='none')
plot(plot_mating_prob)

top_y<-4.5
plot_mating_prob<-ggplot()+
  geom_vline(xintercept=0,color="black",size=0.1)+
  geom_line(data=trunc_distribution_df %>% filter(LA=="rfmix"),aes(x=x,y=value,color=anc),size=1.2)+
  geom_text(data=v.line.data %>% filter(LA=="rfmix"),aes(x=-1,y=top_y-2-0.1*top_y*as.numeric(anc),color=anc,label=paste("AM",as.character(anc),"=",as.character(round(AMvalue,digits=2)),sep="")),size=3,hjust=0)+
  geom_text(data=v.line.data %>% filter(LA=="rfmix"),aes(x=SB+(SB/abs(SB)*0.05),y=top_y+0.2-0.1*top_y*as.numeric(anc),color=anc,label=paste("SB",as.character(anc),"=",as.character(round(SB,digits=2)),sep=""),hjust=-(SB/abs(SB)-1)/2),size=3)+
  coord_fixed(ratio=1/top_y)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"))+
  scale_color_manual(values=c("#b84a62","#bdb246","#437c90"))+
  geom_vline(data=v.line.data%>% filter(LA=="rfmix"),mapping=aes(xintercept=SB,color=anc),size=0.6,linetype="dotted")+
  facet_grid(population~pulses)+
  ylim(0,top_y)+
  ylab("Density")+
  xlab("Sex Bias (ancestry proportion of male - ancestry proportion of female)")+
  ggtitle(("Inferred mating probabilities"))+
  theme_minimal()+
  theme(strip.text.y.right=element_text(angle=0))+
  theme(legend.position='none')
plot(plot_mating_prob)

ggsave(filename=paste0(path,"Mating_Probs_bins_22_normoption_0_onlyfemales.pdf"),plot=plot_mating_prob,height=10.0,width=6.0,units="in")




####
###
###
####




gen_t2<-10
gen_t3<-19

adm_trajectories_allpop<-as.data.frame(matrix(nrow=0,ncol=5))
for (pop_index in c(1:6)) {
  POPname<-c("CLM","MXL","PEL","PUR","ACB","ASW")[pop_index]
  POP1<-c(0.082,0.043,0.032,0.138,0.881,0.765)[pop_index]
  POP2<-c(0.277,0.5,0.761,0.15,0.005,0.041)[pop_index]
  POP3<-1-(POP1+POP2)
  predicted_true <- all_pred_true %>% filter(LA=="rfmix",population==POPname,bins==22,option==0) %>% select(param,value,NN) 
  names(predicted_true)<-c("stat","value","run")
  predicted_false <- all_pred_false %>% filter(LA=="rfmix",population==POPname,bins==22,option==0) %>% select(param,value,NN) 
  names(predicted_false)<-c("stat","value","run")
  
  predicted_true_summary<-as.data.frame(predicted_true %>% group_by(stat)%>% summarize(mean=mean(value)))
  predicted_true_summary<-rbind(predicted_true_summary,c("SB3",1-(predicted_true_summary$mean[predicted_true_summary$stat=="SB1"]+(predicted_true_summary$mean[predicted_true_summary$stat=="SB2"]))))
  predicted_true_summary$stat<-c(rep("AM",3),rep("prop_t1",3),rep("SB",3))
  predicted_true_summary$anc<-c(1,2,3)
  predicted_true_summary$pop<-POPname
  predicted_true_summary$pulses<-"true"
  predicted_false_summary<-as.data.frame(predicted_false %>% group_by(stat)%>% summarize(mean=mean(value)))
  predicted_false_summary<-rbind(predicted_false_summary,c("SB3",1-(predicted_false_summary$mean[predicted_false_summary$stat=="SB1"]+(predicted_false_summary$mean[predicted_false_summary$stat=="SB2"]))))
  predicted_false_summary$stat<-c(rep("AM",3),rep("SB",3))
  predicted_false_summary$anc<-c(1,2,3)
  predicted_false_summary$pop<-POPname
  predicted_false_summary$pulses<-"false"
  predicted_summary<-rbind(predicted_false_summary)#,predicted_true_summary)
  
  sin_mod<-function(x,x_width,x_shift,y_height,y_shift){(((sin(x*pi/(x_width)-(pi/2*(1+x_shift/x_width*2))))/2+1/2)*y_height)+y_shift}
  xmin<- -3
  xmax<-9
  gen_t2<-10
  gen_t3<-19
  
  admixture_curves_pulses_false<-function(x,anc){
    return_vector<-rep(NA,length(x))
    propanc<-c(POP1,POP2,POP3)[anc]
    AM<-as.numeric(predicted_false %>% group_by(stat)%>% summarize(mean=mean(value)) %>% filter(stat==paste("AM",as.character(anc),sep="")) %>% select(mean))
    condition_1<-sin_mod(x,x_width=gen_t3*(1-AM),x_shift=0,y_height=propanc,y_shift=0)
    condition_2<-rep(propanc,length(x))
    return_vector[x<(1-AM)*gen_t3]<-condition_1[x<(1-AM)*gen_t3]
    return_vector[x>=(1-AM)*gen_t3]<-condition_2[x>=(1-AM)*gen_t3]
    return(return_vector)
  }
  
  admixture_curves_pulses_true<-function(x,anc){
    return_vector<-rep(NA,length(x))
    propanc<-c(POP1,POP2,POP3)[anc]
    AM<-as.numeric(predicted_true %>% group_by(stat)%>% summarize(mean=mean(value)) %>% filter(stat==paste("AM",as.character(anc),sep="")) %>% select(mean))
    prop_t1<-as.numeric(predicted_true %>% group_by(stat)%>% summarize(mean=mean(value)) %>% filter(stat==paste("prop_pop",as.character(anc),"_t1",sep="")) %>% select(mean))
    condition_1<-sin_mod(x,x_width=gen_t2*(1-AM),x_shift=0,y_height=prop_t1*propanc,y_shift=0)
    condition_2<-rep(prop_t1*propanc,length(x))
    condition_3<-sin_mod(x,x_width=(gen_t3-gen_t2)*(1-AM),x_shift=gen_t2,y_height=(1-prop_t1)*propanc,y_shift=prop_t1*propanc)
    condition_4<-rep(propanc,length(x))
    return_vector[x<(1-AM)*gen_t2]<-condition_1[x<(1-AM)*gen_t2]
    return_vector[((x>=(1-AM)*gen_t2)&(x<gen_t2))]<-condition_2[((x>=(1-AM)*gen_t2)&(x<gen_t2))]
    return_vector[((x>=gen_t2)&(x<(gen_t2+(1-AM)*(gen_t3-gen_t2))))]<-condition_3[((x>=gen_t2)&(x<(gen_t2+(1-AM)*(gen_t3-gen_t2))))]
    return_vector[(x>=(gen_t2+(1-AM)*(gen_t3-gen_t2)))]<-condition_4[(x>=(gen_t2+(1-AM)*(gen_t3-gen_t2)))]
    return(return_vector)
  }
  admixture_curves_pulses_true(seq(0,gen_t3,0.1),2)
  
  adm_trajectories<-matrix(ncol=4,nrow=0)
  for (anc in c(1,2,3)){
    adm_trajectories<-rbind(adm_trajectories,cbind(seq(0,gen_t3,0.1),admixture_curves_pulses_false(seq(0,gen_t3,0.1),anc),anc,0))
  }
  for (anc in c(1,2,3)){
    adm_trajectories<-rbind(adm_trajectories,cbind(seq(0,gen_t3,0.1),admixture_curves_pulses_true(seq(0,gen_t3,0.1),anc),anc,1))
  }
  
  adm_trajectories<-as.data.frame(adm_trajectories,stringsAsFactors=FALSE)
  names(adm_trajectories)<-c("generation","anc_prop","anc","pulses")
  adm_trajectories$anc<-as.factor(adm_trajectories$anc)
  adm_trajectories$pulses[adm_trajectories$pulses==0]<-"One pulse 19 generations ago"
  adm_trajectories$pulses[adm_trajectories$pulses==1]<-"Two pulses 19 and 9 generations ago"
  adm_trajectories$pulses<-as.factor(adm_trajectories$pulses)
  adm_trajectories$pop<-POPname
  adm_trajectories_allpop<-rbind(adm_trajectories_allpop,adm_trajectories)
}

plot_adm_traj<-ggplot(adm_trajectories_allpop,aes(x=generation,y=anc_prop,color=anc))+
  geom_line(size=1.2)+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"),guide='none')+
  scale_color_manual(name="Ancestry",labels=c("sub-Saharan African","Native American","European"),values=c("#b84a62","#bdb246","#437c90"),guide='none')+
  ylab("Proportion of the present gene pool")+
  xlab("Generation")+
  ggtitle(paste("Migration pulses"))+
  theme_minimal()+
  facet_grid(pop~pulses)
plot(plot_adm_traj)

adm_trajectories_allpop$cum_anc_prop<- adm_trajectories_allpop$anc_prop
adm_trajectories_allpop$cum_anc_propadm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==1]<-adm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==1]+
  adm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==2]+
  adm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==3]
adm_trajectories_allpop$cum_anc_propadm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==2]<-adm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==1]+
  adm_trajectories_allpop$anc_prop[adm_trajectories_allpop$anc==2]



plot_adm_traj<-ggplot(adm_trajectories_allpop,aes(x=generation-19,y=cum_anc_prop,color=anc,fill=anc))+
  geom_area()+
  scale_fill_manual(values=c("#b84a62","#bdb246","#437c90"),guide='none')+
  scale_color_manual(name="Ancestry",labels=c("sub-Saharan African","Native American","European"),values=c("#b84a62","#bdb246","#437c90"),guide='none')+
  ylab("Proportion of the present gene pool")+
  xlab("Generation")+
  coord_fixed(ratio=9)+
  ggtitle(paste("Migration pulses"))+
  theme_minimal()+
  geom_vline(xintercept=-19,col="black",linetype="dotted")+
  geom_vline(xintercept=-9,col="black",linetype="dotted")+
  facet_grid(pop~pulses, switch="y")+
  scale_y_continuous(position="right")+
  theme(strip.text.y.left=element_text(angle=0))
plot(plot_adm_traj)


ggsave(filename=paste0(path,"Adm_Traj_bins_22_normoption_0_onlyfemales.pdf"),plot=plot_adm_traj,height=10.0,width=6.0,units="in") 



#######



all_pred_false_plot %>% filter(bins==22,option==0) %>% 
  select( -ancestry, -bins, -option, -value_gnomix) %>% 
  rename(value=value_rfmix,pop=population,parameter=param) %>% 
  group_by(pop,parameter) %>% summarise(value=mean(value)) %>%
  mutate(value=ifelse(((parameter=="AM1")|(parameter=="AM2")|(parameter=="AM3")),3^((1-value)*7-4),value)) %>%
  write.table(., file=paste0(path,"mean_parameter_allpop_false_3oct2023.txt"), row.names=FALSE, col.names=TRUE,quote=FALSE,sep=" ")


all_pred_true_plot %>% filter(bins==22,option==0) %>% 
  select( -ancestry, -bins, -option, -value_gnomix) %>% 
  rename(value=value_rfmix,pop=population,parameter=param) %>% 
  group_by(pop,parameter) %>% summarise(value=mean(value)) %>%
  mutate(value=ifelse(((parameter=="AM1")|(parameter=="AM2")|(parameter=="AM3")),3^((1-value)*7-4),value))%>%
  write.table(., file=paste0(path,"mean_parameter_allpop_true_3oct2023.txt"), row.names=FALSE, col.names=TRUE,quote=FALSE,sep=" ")

