############################################################################################################
# Version: 22 January 2024
# 
# Author: Shuli Chen
#
# When using these data please cite the original article:
#      
#      
# See: http://
#
############################################################################################################
gc()
rm(list=ls(all=TRUE))

#Packages required:  
library(mgcv); library(qgam); library(mgcViz);library(gam);library(plm);library(LaplacesDemon);library(gtable);library(grid)
library(dplyr);library(gratia);library(rlang);require(fields);library(MASS);library(smatr);library(ggplot2);library(ggpmisc);
library(Rmisc) ;library(lattice);library(plyr);library(data.table);library(ggsignif);library(dplyr);library(scales);
library(grid);library(ggthemes);library(RColorBrewer);library(plotly);library(MuMIn);library(Hotelling);
library(mvdalab);library(gamlss);library(gamlss.dist);library(gamlss.add);library(fitdistrplus);library(logspline);
library(pracma);library(lubridate);library(r2d2);library(paletteer);library(viridis);library(patchwork)

# set your own working directory
setwd("F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests")

############################################################################################################
# DATASETS REQUIRED:
############################################################################################################

## These are the remote sensed 1km data; see Readme documentation for definition of variables 
## functions are listed at the end of the page
## Read EVI anomaly and climate data 
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
Drought_year='2005'                                                                            
file_name=paste('Droughtof2005_EVI_Climate_WTD_simple','.csv',sep='')#_simple
file_ful_path=paste(file_path,"/",file_name,sep='')
file_ful_path
Drought.data_2005<- read.csv(file=file_ful_path,header=T) 
summary(Drought.data_2005)


####-------------------------------------------main code for Figure 2-------------------------------------------
##------------------------------------Panel B Original EVI Anomaly with Severity-------------------------------- 
options(digits=7) 
Drought_ad.data_2005<-data.frame(Drought.data_2005)
Drought_ad.data_2005 <- Drought_ad.data_2005[which(is.finite(Drought_ad.data_2005$EVI_anomaly) ),]
#calculate mean value of EVI anomaly for shallow and deep water tables
SWTD_MeanEVI<-mean(Drought_ad.data_2005[which(Drought_ad.data_2005$WTD <=8 ),]$EVI_anomaly)
DWTD_MeanEVI<-mean(Drought_ad.data_2005[which(Drought_ad.data_2005$WTD > 22 ),]$EVI_anomaly)
Drought_ad.data_2005<-Drought_ad.data_2005[which(Drought_ad.data_2005$WTD <= 40 & Drought_ad.data_2005$HAND_CLASS <=40 & Drought_ad.data_2005$HAND_CLASS >=0 ),] 
Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >= 1),]$HAND_CLASS=Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >= 1),]$HAND_CLASS-1
# summarise EVI anomaly binned by HAND
Summary_2005=summary_group_severity(Drought_ad.data_2005,0) #calculate for original EVI
# for drawing the right order from Modest to Severe
Summary_2005$Drought_severity=""
Summary_2005[which((Summary_2005$Drought_Condition =="Modest Drought") ),]$Drought_severity="Meodest_Ori"
Summary_2005[which((Summary_2005$Drought_Condition =="Medium Drought") ),]$Drought_severity="Moedium_Ori"
Summary_2005[which((Summary_2005$Drought_Condition  =="Severe Drought")),]$Drought_severity="Severe_Ori"

rangeMin=-0.5
rangeMax=rangeMin+1
BY_increase=0.5
Figure2_PanleB<- ggplot_adjustment_withSeverity(Summary_2005,rangeMin,rangeMax,BY_increase,Drought_year,SWTD_MeanEVI,DWTD_MeanEVI) 
Figure2_PanleB
##--------------------------------------------------------------------------------------------------------------

##------------------------------------Panel C Original EVI Anomaly VS Par class--------------------------------- 
Drought_ad.data_2005<-data.frame(Drought.data_2005)
Drought_ad.data_2005<-Drought_ad.data_2005[which(is.finite(Drought_ad.data_2005$EVI_anomaly) & 
                                        is.finite(Drought_ad.data_2005$PAR_anomaly) ),]
# Grouped by par 
Drought_ad.data_2005$PAR_class<-0
Par_class_start<-(-1.25)
group_interval=0.5
group_number=7
for (i in 1:group_number) {
  if (i == 1) { Drought_ad.data_2005[which(Drought_ad.data_2005$PAR_anomaly<=Par_class_start+(i-1)*group_interval),]$PAR_class=(Par_class_start+(i-1.5)*group_interval)}
  else{ if (i== group_number ) {Drought_ad.data_2005[which(Drought_ad.data_2005$PAR_anomaly>(Par_class_start+(i-2)*group_interval)),]$PAR_class=Par_class_start+(i-1.5)*group_interval}
    else{Drought_ad.data_2005[which(Drought_ad.data_2005$PAR_anomaly>(Par_class_start+(i-2)*group_interval) & Drought_ad.data_2005$PAR_anomaly<=(Par_class_start+(i-1)*group_interval) ),]$PAR_class=Par_class_start+(i-1.5)*group_interval}
  }
}
Drought_ad.data_2005_1<-Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >=0 &Drought_ad.data_2005$HAND_CLASS <8),]
Drought_ad.data_2005_2<-Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >=8 & Drought_ad.data_2005$HAND_CLASS<16),]
Drought_ad.data_2005_3<-Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >=16 & Drought_ad.data_2005$HAND_CLASS<24),]
Drought_ad.data_2005_4<-Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >=24 &Drought_ad.data_2005$HAND_CLASS <32),]
Drought_ad.data_2005_5<-Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >=32 &Drought_ad.data_2005$HAND_CLASS <40),]
Drought_ad.data_2005_6<-Drought_ad.data_2005[which(Drought_ad.data_2005$HAND_CLASS >=40 &Drought_ad.data_2005$HAND_CLASS <120),]

tgc_1 <- summarySE(Drought_ad.data_2005_1, measurevar="EVI_anomaly",  groupvars=c("PAR_class"))  #
tgc_1$Type="WTD_08m"
tgc_2<- summarySE(Drought_ad.data_2005_2, measurevar="EVI_anomaly", groupvars=c("PAR_class")) 
tgc_2$Type="WTD_16m"
tgc_3 <- summarySE(Drought_ad.data_2005_3, measurevar="EVI_anomaly", groupvars=c("PAR_class")) 
tgc_3$Type="WTD_24m"
tgc_4 <- summarySE(Drought_ad.data_2005_4, measurevar="EVI_anomaly", groupvars=c("PAR_class")) 
tgc_4$Type="WTD_32m"
tgc_5 <- summarySE(Drought_ad.data_2005_5, measurevar="EVI_anomaly", groupvars=c("PAR_class")) 
tgc_5$Type="WTD_40m"
tgc_6 <- summarySE(Drought_ad.data_2005_6, measurevar="EVI_anomaly", groupvars=c("PAR_class")) 
tgc_6$Type="WTD_D>40m"

tgc_combine<- rbind(tgc_1,tgc_2,tgc_3,tgc_4,tgc_5,tgc_6)
rangeMin=-0.7
rangeMax=rangeMin+1.5
BY_increase=0.5
cbPalette <- c("#26456EFF", "#1F74B1FF", "#6AB1D6FF", "#BCCACFFF", "#EEB78DFF", "#FD8E3FFF","#FD8E3FFF")
Figure2_PanleC<-ggplot(tgc_combine, aes(x=PAR_class, y=EVI_anomaly,colour=factor(Type), group=factor(Type)))+ #, colour=Drought_Condition, group=Drought_Condition+
  geom_errorbar(aes(ymin=EVI_anomaly-ci, ymax=EVI_anomaly+ci,group=factor(Type),colour=factor(Type),fill=factor(Type)),  size=0.5, width=.1,face="bold") + #
  geom_path(size=1.0, aes(group=factor(Type),color=factor(Type),fill=factor(Type)))+
  geom_point( size=4.5, shape=23,aes(group=factor(Type),colour=factor(Type),fill=factor(Type))) + 
  geom_hline(yintercept=0)+
  xlab("PAR Anomaly") + #
  ylab("EVI Anomaly") +
  #ggtitle(Title_name) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  scale_y_continuous(limits =c(-0.78,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  scale_x_continuous(limits =c(-1.75,1.75), breaks=seq(-1.5,1.5,  by=0.5))+ #-1.5,1.5 1,5
  theme_bw() +
  theme( legend.position=c(1,0))+
  theme(legend.justification=c(1,0))+
  theme_classic() + # 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
  theme(axis.title=element_text(size=24))+ #,face="bold"
  theme(plot.title=element_text(size=24))
Figure2_PanleC

##--------------------------------------------------------------------------------------------------------------

##------------------------------------Panel C PAR Distribution--------------------------------------------------
# PAR distribution at 1km categorised by drought severity
Drought_ad.data_2005<-data.frame(Drought.data_2005)
Drought_ad.data_2005<-Drought_ad.data_2005[which( Drought_ad.data_2005$WTD<=40 & Drought_ad.data_2005$WTD >= 0),] # Severe Drought Medium Drought
rangeMin=0
rangeMax=200000#
BY_increase=100000#
group_number=14
cbPalette1 <- c( "#FEE08B" ) #"#D53E4F",  "#F46D43", "#FDAE61", "#FEE08B", 
# choose either 'Modest Drought', 'Medium Drought', 'Severe Drought'
Fig2_PanelCDistribution<-Draw_Fig2_PanelCDistribution(Drought_ad.data_2005,DroughtSeverity = 'Modest Drought',
                                                      x0.var = 'PAR_anomaly',x1.var = 'Drought_Condition',
                                                      rangeMin,rangeMax,BY_increase,group_number,cbPalette1)
  
Fig2_PanelCDistribution
##---------------------------------------------------------------------------------------------------------------





##------------------------Functions for figure 2-----------------------------------------------------------------
summary_group_severity<-function(Input_model.data, group_flag){
 
  if (group_flag ==1) {
    Residual_summary <- summarySE(Input_model.data, measurevar="Residual_R", groupvars=c("HAND_CLASS","Drought_Condition"))  
  } else{ if (group_flag ==0) {
    Residual_summary <- summarySE(Input_model.data, measurevar="EVI_anomaly", groupvars=c("HAND_CLASS","Drought_Condition")) 
  }else {
    Residual_summary <- summarySE(Input_model.data, measurevar="EVI_Prediction", groupvars=c("HAND_CLASS","Drought_Condition")) 
  }
  }
  return(Residual_summary)
}

# function for  Figure 2 Panel B 
ggplot_adjustment_withSeverity<-function(tgca,rangeMin,rangeMax,BY_increase,Drought_year,SWTD_MeanEVI,DWTD_MeanEVI){
  
  Title_name=paste("       ",Drought_year,"Drought")
  
    cbPalette <- c("#FED976", "#FE9929","#D94801","#F0F0F0", "#D9D9D9","#969696") 

  ggplot(tgca, aes(x=HAND_CLASS, y=EVI_anomaly, colour = factor(Drought_severity),group=factor(Drought_severity)))+ #, colour=Drought_Condition, group=Drought_Condition))+ , group=tgc1.type tgc1.Droughtlength
    geom_errorbar(aes(ymin=EVI_anomaly-ci, ymax=EVI_anomaly+ci,colour = factor(Drought_severity),group=factor(Drought_severity)), size=0.5, width=.3,face="bold") +
    geom_point( size=4.5, shape=23, aes(group=factor(Drought_severity),fill=factor(Drought_severity))) + # 21is filled circle , fill="white"
    scale_fill_manual(values=c( "#FED976", "#FE9929","#D94801", "#F0F0F0","#D9D9D9", "#969696")) + 
    geom_hline(yintercept = 0,  size=1, linetype="dashed") +
    geom_hline(yintercept = SWTD_MeanEVI,  size=1, linetype="dashed", color="dark green") +
    geom_hline(yintercept =DWTD_MeanEVI,  size=1, linetype="dashed",color="orange") +
    
    xlab("HAND (meter)") +
    ylab("EVI anomaly   ") +
    
    scale_colour_manual(values=cbPalette)  +                     
    ggtitle(Title_name) +
    scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
    scale_x_continuous(breaks=seq(0, 40,10))+
    theme_few() %+replace% 
    theme(panel.background = element_rect(fill = NA,colour = "black", 
                                          size =0.1)) +
    theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
    theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
    theme(axis.title=element_text(size=24))+ #,face="bold"
    theme(plot.title=element_text(size=24)) #face="bold"
}


# function for  Figure 2 Panel C Distributions 
Draw_Fig2_PanelCDistribution<-function(Data,DroughtSeverity = 'Modest Drought',x0.var = 'PAR_anomaly',x1.var = 'Drought_Condition',rangeMin,rangeMax,BY_increase,group_number, cbPalette){
  
  ix0 <- match(x0.var, names(Data))
  ix1 <- match(x1.var, names(Data)) 
  print(ix0)
  Data<-Data[which(is.finite(Data[,ix0]) & Data[,ix1]==DroughtSeverity),] # Severe Drought, Medium Drought,Modest Drought
  
  Data$PAR_class<-0
  Par_class_start<-(-1.5)
  group_interval=0.25
  for (i in 1:group_number) {
    if (i == 1) { Data[which(Data[,ix0]<=Par_class_start+(i-1)*group_interval),]$PAR_class=(Par_class_start+(i-1.5)*group_interval)}
    else{ if (i== group_number ) {Data[which(Data[,ix0]>(Par_class_start+(i-2)*group_interval)),]$PAR_class=Par_class_start+(i-1.5)*group_interval}
          else{Data[which(Data[,ix0]>(Par_class_start+(i-2)*group_interval) & Data[,ix0]<=(Par_class_start+(i-1)*group_interval) ),]$PAR_class=Par_class_start+(i-1.5)*group_interval}
      }
  }

  grid.newpage()
  p_Distribution <- ggplot(data=Data, mapping=aes(x=PAR_class     ))+
    geom_bar(stat="count",colour="#737373",fill = cbPalette[1], width=0.25,size=0.6, position =  position_dodge())+  #
    scale_fill_manual(values =cbPalette[1] )+ #cbPalette[1]
    scale_x_continuous(limits =c(-1.75,1.75), breaks=seq(-1.5,1.5,  by=0.50))+
    scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +    theme_few() %+replace% 
    theme(panel.background = element_rect(fill = NA,colour = "black", 
                                          size =1))+
    ylab("Area (km2)") + 
    xlab("PAR Anomaly") + 
    theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
    theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
    theme(axis.title=element_text(size=24))+ #,face="bold"
    theme(plot.title=element_text(size=24))
  return(p_Distribution)
}
##---------------------------------------------------------------------------------------------------------------                                                                 

