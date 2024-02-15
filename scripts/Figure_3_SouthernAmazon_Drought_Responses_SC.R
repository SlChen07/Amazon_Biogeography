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
library(mgcv); library(qgam); library(mgcViz);library(gam);library(plm);library(LaplacesDemon);library(gtable);library(grid);
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

## These are the remote sensed aggregated 0.4 degree data; see documentation for definition of variables
## functions are listed at the end of the page
## Read EVI anomaly and climate data
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
file_name=paste('SouthernAmazon_EVI_Climate_04De','.csv',sep='')
file_ful_path=paste(file_path,"/",file_name,sep='')
Drought_04De.data<- read.csv(file=file_ful_path,header=T) 


####-------------------------------------------main code for Figure 3-------------------------------------------
##------------------------------------Panel A Corrected EVI Anomaly in three droughts---------------------------
## establish a unified model for three drought; see model selection below
seg_target=12
Drought_04De_ad.data<-Drought_04De.data[which(Drought_04De.data$SegGeo_number >= 9 & Drought_04De.data$SegGeo_number <=21  & Drought_04De.data$WTD <=120),]
Drought_04De_ad_new2.data<-Drought_04De_ad.data[which(is.finite(Drought_04De_ad.data$EVI_anomaly) & is.finite(Drought_04De_ad.data$PAR_anomaly)
                                                      & is.finite(Drought_04De_ad.data$VPD_anomaly) &  (Drought_04De_ad.data$EVI_anomaly != 0)  
                                                      & is.finite(Drought_04De_ad.data$MCWD_anomaly) & is.finite(Drought_04De_ad.data$Pre_anomaly) 
                                                      &  is.finite(Drought_04De_ad.data$SoilFertility)  & Drought_04De_ad.data$SegGeo_number ==seg_target 
                                                      & Drought_04De_ad.data$WTD <=60 ),] # 


Drought_04De_ad_new1.data<-data.frame(Drought_04De_ad_new2.data$PAR_anomaly,Drought_04De_ad_new2.data$VPD_anomaly,+
                                       Drought_04De_ad_new2.data$Pre_anomaly,Drought_04De_ad_new2.data$MCWD_anomaly, +
                                       Drought_04De_ad_new2.data$SoilFertility,Drought_04De_ad_new2.data$SoilClay_content, +
                                        Drought_04De_ad_new2.data$WTD,Drought_04De_ad_new2.data$Tree_Height,+
                                       Drought_04De_ad_new2.data$Drought_Length,Drought_04De_ad_new2.data$MCWD_STD,+
                                       Drought_04De_ad_new2.data$WTD, Drought_04De_ad_new2.data$EVI_anomaly, Drought_04De_ad_new2.data$year,+
                                       Drought_04De_ad_new2.data$SegGeo_number,Drought_04De_ad_new2.data$HAND_CLASS)#

names(Drought_04De_ad_new1.data)<-c('PAR_anomaly','VPD_anomaly','Pre_anomaly','MCWD_anomaly',  
                                   'SoilFertility','SoilClay_content','HAND',   
                                   'TreeHeight','Drought_Length','MCWD_STD', 'WTD',
                                   'EVI_anomaly','year','SegGeo_number','HAND_CLASS')
summary(Drought_04De_ad_new1.data)

Drought_04De_ad_new1.data<-Drought_04De_ad_new1.data[which(Drought_04De_ad_new1.data$HAND >=0 & Drought_04De_ad_new1.data$HAND <=63 &    Drought_04De_ad_new1.data$HAND_CLASS >=-0 ),] #

#----------------------------------------------------------
Drought_04De_ad_new1.data$SegGeo_numberf<-factor(Drought_04De_ad_new1.data$SegGeo_number)
Drought_04De_ad_new.data2<-Drought_04De_ad_new1.data  #before scaling
Drought_04De_ad_new.data_before_scaling<-Drought_04De_ad_new1.data 

#--Scale the data---------------------------
Scaled.Drought_04De_ad_new1.data<-apply(Drought_04De_ad_new1.data[,1:10],2, scale)
Drought_04De_ad_new.data <- as.data.frame( cbind(Scaled.Drought_04De_ad_new1.data,Drought_04De_ad_new1.data[,11:15]))

#threshold_scale=10
#Drought_04De_ad_new.data <- Drought_04De_ad_new.data[which(abs(Drought_04De_ad_new.data$EVI_anomaly) <= 20 &   
                                                               #abs(Drought_04De_ad_new.data$PAR_anomaly)<= threshold_scale & 
                                                               #abs(Drought_04De_ad_new.data$VPD_anomaly)<= threshold_scale & 
                                                               #abs(Drought_04De_ad_new.data$Pre_anomaly)<= threshold_scale & 
                                                               #abs(Drought_04De_ad_new.data$MCWD_anomaly)<= threshold_scale & 
                                                               #abs(Drought_04De_ad_new.data$SoilFertility)<= threshold_scale & 
                                                               #abs(Drought_04De_ad_new.data$HAND)<=10 & 
                                                               #abs(Drought_04De_ad_new.data$Drought_Length)<= threshold_scale & 
                                                               #abs(Drought_04De_ad_new.data$MCWD_STD)<= threshold_scale),] #  


Ancova_model.lm.slope_model_climate=ancova_establish_slope_Guiana(Drought_04De_ad_new.data) # see model selection below
#Drought_04De_adj_new.data$Prediction=fitted(Ancova_model.lm.slope_model_climate)
summary(Ancova_model.lm.slope_model_climate)
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
gam.check(Ancova_model.lm.slope_model_climate)
AIC(Ancova_model.lm.slope_model_climate)
#concurvity(Ancova_model.lm.slope_model_climate, full = FALSE)
#b <- getViz(Ancova_model.lm.slope_model_climate)
#print(plot(b, allTerms = T,xlim=range(-1.2:3),ylim=range(-2:2) ), pages = 1)
##--------------------------------------------------------------------------------------------------------------

##------------------------------------Panel A observations continued--------------------------------------------
##---------------------If you want to recalculate, please start it here-----------------------------------------
seg=seg_target
Drought_04De_ad_new.data_ha<-Drought_04De_ad_new.data

#calculate the HAND Class
Drought_04De_ad_new.data_ha<-Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$WTD<=60 ),]
Drought_04De_ad_new.data_ha$HAND_CLASS<-floor(Drought_04De_ad_new.data_ha$WTD)+1 ## 0.4 Degree
Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS >= 1),]$HAND_CLASS=Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS >= 1),]$HAND_CLASS-1
Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS=(Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS%/%2+1)*2
Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS== -1),]$HAND_CLASS=Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS == -1),]$HAND_CLASS-0
Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS=Drought_04De_ad_new.data_ha[which(Drought_04De_ad_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS-1


#Calculate the mean by year
Drought_04De_ad_new.data_ha_mean<-Drought_04De_ad_new.data_ha
#set mean for 2005, 2010, 2015 separately
Drought_04De_ad_new.data_ha_mean<-set_value_byyear(Drought_04De_ad_new.data_ha_mean,x0.var='year',x1.var='PAR_anomaly',x2.var='VPD_anomaly',x3.var='MCWD_anomaly',
                                                   x4.var='Pre_anomaly',x5.var='Drought_Length',drought_year='2005',moment=median)
Drought_04De_ad_new.data_ha_mean<-set_value_byyear(Drought_04De_ad_new.data_ha_mean,x0.var='year',x1.var='PAR_anomaly',x2.var='VPD_anomaly',x3.var='MCWD_anomaly',
                                                   x4.var='Pre_anomaly',x5.var='Drought_Length',drought_year='2010',moment=median)
Drought_04De_ad_new.data_ha_mean<-set_value_byyear(Drought_04De_ad_new.data_ha_mean,x0.var='year',x1.var='PAR_anomaly',x2.var='VPD_anomaly',x3.var='MCWD_anomaly',
                                                   x4.var='Pre_anomaly',x5.var='Drought_Length',drought_year='2015',moment=mean)
#calculate the prediction according to mean climate of each drought
Drought_04De_ad_new.data_ha_mean<-add_fitted(Drought_04De_ad_new.data_ha_mean,Ancova_model.lm.slope_model_climate)
colnames(Drought_04De_ad_new.data_ha_mean)[16] = 'EVI_anomaly_fit_mean'#'EVI_anomaly_fit_mean_DL' 

#calculate the prediction 
Drought_04De_ad_new.data_ha<-add_fitted(Drought_04De_ad_new.data_ha,Ancova_model.lm.slope_model_climate)
colnames(Drought_04De_ad_new.data_ha)[16] = 'EVI_anomaly_fit'
#corrected EVI anomaly for 2005 drought
Drought_04De_ad_new.data_ha_2005<-Drought_EVI_correction(Drought_04De_ad_new.data_ha,Drought_04De_ad_new.data_ha_mean,
                                                         x0.var='year',x1.var='EVI_anomaly',x2.var='HAND_CLASS',data1.var='EVI_anomaly_fit',
                                                         data2.var='EVI_anomaly_fit_mean',drought_year='2005')
arr_observation_2005<-summary_group_full(Drought_04De_ad_new.data_ha_2005,x0.var='EVI_anomaly',
                                         x1.var='EVI_anomaly_corrected', flag=1)
arr_observation_2005$combinetype="2005_year"

#corrected EVI anomaly for 2010 drought
Drought_04De_ad_new.data_ha_2010<-Drought_EVI_correction(Drought_04De_ad_new.data_ha,Drought_04De_ad_new.data_ha_mean,
                                                         x0.var='year',x1.var='EVI_anomaly',x2.var='HAND_CLASS',data1.var='EVI_anomaly_fit',
                                                         data2.var='EVI_anomaly_fit_mean',drought_year='2010')
arr_observation_2010<-summary_group_full(Drought_04De_ad_new.data_ha_2010,x0.var='EVI_anomaly',
                                         x1.var='EVI_anomaly_corrected', flag=1)
arr_observation_2010$combinetype="2010_year"

#corrected EVI anomaly for 2015 drought
Drought_04De_ad_new.data_ha_2015<-Drought_EVI_correction(Drought_04De_ad_new.data_ha,Drought_04De_ad_new.data_ha_mean,
                                                         x0.var='year',x1.var='EVI_anomaly',x2.var='HAND_CLASS',data1.var='EVI_anomaly_fit',
                                                         data2.var='EVI_anomaly_fit_mean',drought_year='2015')
arr_observation_2015<-summary_group_full(Drought_04De_ad_new.data_ha_2015,x0.var='EVI_anomaly',
                                         x1.var='EVI_anomaly_corrected', flag=1)
arr_observation_2015$combinetype="2015_year"

#calculate the trend lines for each drought
lm.HD_Line_2005_new <-(lm( EVI_anomaly ~ HAND_CLASS,data=arr_observation_2005))
lm.GAM_Line_2005_new<-coef(lm.HD_Line_2005_new)
lm.HD_Line_2010_new <-(lm( EVI_anomaly ~ HAND_CLASS,data=arr_observation_2010))
lm.GAM_Line_2010_new<-coef(lm.HD_Line_2010_new)
lm.HD_Line_2015_new <-(lm( EVI_anomaly ~ HAND_CLASS,data=arr_observation_2015))
lm.GAM_Line_2015_new<-coef(lm.HD_Line_2015_new)
arr_ori<-rbind(arr_observation_2005,arr_observation_2010,arr_observation_2015)# arr_ori_1, arr_ori_2, arr_ori_3
#colnames(arr_ori)[5] = 'combinetype'

#Draw the figure
rangeMin=-2.5
range_mm=3.6
rangeMax=rangeMin+range_mm
BY_increase=0.5
cbPalette <- c("#006837",  "#88419D",  "#08519c", "#08519c","#D9D9D9","#969696")
size_number=1
Figure3_PanelA_observations<-ggplot(arr_ori, aes(x=HAND_CLASS, y=EVI_anomaly, group=factor(combinetype)))+ 
  
  geom_errorbar(aes(ymin=EVI_anomaly-ci, ymax=EVI_anomaly+ci,colour = factor(combinetype),group=factor(combinetype)), size=0.5, width=.3) +  #,face="bold"
  geom_point( size=3.5, shape=23,aes(group=factor(combinetype),colour=factor(combinetype),fill=factor(combinetype))) + 
  geom_hline(yintercept = 0,width=0.6) +
  geom_abline(intercept =lm.GAM_Line_2005_new[1], slope = lm.GAM_Line_2005_new[2],colour=c("#006837"),size=size_number)+
  geom_abline(intercept =lm.GAM_Line_2010_new[1], slope = lm.GAM_Line_2010_new[2],colour=c("#88419D"),size=size_number)+
  geom_abline(intercept = lm.GAM_Line_2015_new[1], slope = lm.GAM_Line_2015_new[2],colour=c("#08519c"),size=size_number)+
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  xlab("HAND (meter)") +
  ylab("      EVI anomaly  ") +
  scale_y_continuous(limits =c(-2.6,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  scale_x_continuous(limits =c(0,40), breaks=seq(0,40,  by=10))+
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

Figure3_PanelA_observations
#---------------------------------------------------------------------------------------------------------------


##------------------Panel A prediction continued----------------------------------------------------------------

#Draw prediction with for each drought (0.4De)
seg=12
Drought_04De_ad_new.data_2005<-Drought_04De_ad_new.data[which(Drought_04De_ad_new.data$SegGeo_number ==seg & Drought_04De_ad_new.data$year=="2005")  ,] #
WTD_arr=c((-180:360)/20) #
number_group=541
new0_2005<-Drought_EVI_ModelPrediction(Drought_04De_ad_new.data, Drought_04De_ad_new.data_before_scaling,Ancova_model.lm.slope_model_climate,
                                       WTD_arr,x0.var='year',x1.var='PAR_anomaly',
                                       x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                       x5.var='Drought_Length',x6.var='HAND',drought_year='2005',moment=median,number_group)
new0_2005$type="2005_year"  
new0_2010<-Drought_EVI_ModelPrediction(Drought_04De_ad_new.data, Drought_04De_ad_new.data_before_scaling,Ancova_model.lm.slope_model_climate,
                                       WTD_arr,x0.var='year',x1.var='PAR_anomaly',
                                       x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                       x5.var='Drought_Length',x6.var='HAND',drought_year='2010',moment=median,number_group)
new0_2010$type="2010_year"
new0_2015<-Drought_EVI_ModelPrediction(Drought_04De_ad_new.data, Drought_04De_ad_new.data_before_scaling,Ancova_model.lm.slope_model_climate,
                                       WTD_arr,x0.var='year',x1.var='PAR_anomaly',
                                       x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                       x5.var='Drought_Length',x6.var='HAND',drought_year='2015',moment=median,number_group)
new0_2015$type="2015_year"


tgca_combine<-rbind(new0_2005,new0_2010,new0_2015) #
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0.0 & tgca_combine$WTD_reverse<=40),]

rangeMin=-2.5#-1.5
range_mm=3.6#2.4
rangeMax=rangeMin+range_mm
BY_increase=0.5

lm.HD_Line_2005_new <-(lm(fit ~ WTD_reverse,data=tgca_combine[which(tgca_combine$type=="2005_year"),]))
lm.GAM_Line_2005_new<-coef(lm.HD_Line_2005_new)
lm.HD_Line_2010_new <-(lm(fit ~ WTD_reverse,data=tgca_combine[which(tgca_combine$type=="2010_year"),]))
lm.GAM_Line_2010_new<-coef(lm.HD_Line_2010_new)
lm.HD_Line_2015_new <-(lm(fit ~ WTD_reverse,data=tgca_combine[which(tgca_combine$type=="2015_year"),]))
lm.GAM_Line_2015_new<-coef(lm.HD_Line_2015_new)


cbPalette <- c("#006837","#88419D","#08519c", "#D9D9D9","#969696") 
size_number=1
Figure3_PanelAPrediction<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  # scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  scale_y_continuous(limits =c(-2.6,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  geom_line(color="red") +
  geom_hline(yintercept = 0,  width=.6) +
  
  geom_smooth(aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") + #, show.legend = F
  
  geom_abline(intercept =lm.GAM_Line_2005_new[1], slope = lm.GAM_Line_2005_new[2],linetype="dashed",colour=c("#006837"),size=size_number)+
  
  geom_abline(intercept =lm.GAM_Line_2010_new[1], slope = lm.GAM_Line_2010_new[2],linetype="dashed",colour=c("#88419D"),size=size_number)+
  
  geom_abline(intercept = lm.GAM_Line_2015_new[1], slope = lm.GAM_Line_2015_new[2],linetype="dashed",colour=c("#08519c"),size=size_number)+
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual( values=cbPalette)  +   
  scale_x_continuous(limits =c(0,40),breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))

Figure3_PanelAPrediction
##------------------Merge Panel A observation and prediction ---------------------------------------------------
Figure3_PanelA<-merge_figures_sameunits(Figure3_PanelA_observations,Figure3_PanelAPrediction)

grid.draw(Figure3_PanelA)
##--------------------------------------------------------------------------------------------------------------


##------------------------------------Model selection-----------------------------------------------------------
# set a full model
mod.SAFull<-mgcv::gam(EVI_anomaly ~ s(HAND,PAR_anomaly,bs = "tp")  + s(VPD_anomaly, Drought_Length, bs = "tp")+ s(MCWD_anomaly, Pre_anomaly, bs = "tp")
  + ti(HAND, Pre_anomaly, bs = "tp") + ti(HAND, MCWD_anomaly, bs = "tp") +  ti(HAND, VPD_anomaly, bs = "tp") + ti(HAND, Drought_Length, bs = "tp")  + 
  ti(PAR_anomaly, Drought_Length, bs = "tp") + ti(PAR_anomaly, Pre_anomaly, bs = "tp") +ti(MCWD_anomaly, Drought_Length, bs = "tp")
  ti(VPD_anomaly, PAR_anomaly, bs = "tp") + ti(VPD_anomaly,Pre_anomaly, bs = "tp")+  ti(PAR_anomaly, MCWD_anomaly, bs = "tp") +
                                                 ti(MCWD_anomaly,VPD_anomaly, bs = "tp")+ti(Pre_anomaly, Drought_Length, bs = "tp"),method = "ML",data = Drought_04De_ad_new.data) #,method = "REML"


  summary(mod.SAFull) 
  #getAllTerms(mod.SAFull)
  #AICc(mod.SAFull)
  AIC(mod.SAFull)
#using dredge to select the model based on AICc/AIC
  options(na.action = "na.fail")
  dd<-dredge(mod.SAFull, trace =2)
  subset(dd, delta < 10)
  #summary(dd[1,])
  #importance(dd)
  summary(mod.SAFull(dd, subset = delta < 10))
  
##--------------------------------------------------------------------------------------------------------------



##--------------------------------------Figure 3 Panel B--------------------------------------------------------
WTD_set_array=((0:6)+0.1)*6.8 # set preferred intervals
PAR_set_array=((-30:30)/15)#
number_group=61
WTD_inital<-0
byincrease<-8
for(i in seq(from=1, 6, by=1))
{
  Model_fit_array<-Drought_EVI_ModelPred_ClimateSensitivity(Drought_04De_ad_new.data, Drought_04De_ad_new.data_before_scaling,Ancova_model.lm.slope_model_climate,
                                                                  WTD_set_array[i],x1.var='PAR_anomaly',
                                                                  x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                                                  x5.var='Drought_Length',x6.var='HAND',key_variable='PAR_anomaly',PAR_set_array, moment=mean,number_group)
  newd0<-Model_fit_array
  if ((i-1)<10) {st_tag=paste(as.character(0),as.character(i-1))}else {
    st_tag=c(as.character(i-1))
  }
  st_tag<-gsub(" ", "", st_tag)
  print(st_tag)
  type_str<-paste('HAND',st_tag,as.character(WTD_inital+byincrease*(i-1.0)),sep = "_")
  newd0$type=type_str
  print(type_str)
  
  if (i == 1 ) {
    newd=rbind(newd0)
  }
  newd<-rbind(newd,newd0)
}

tgca_combine<-newd
BY_increase=0.5
rangeMin=-2.5#
range_mm=3.6#
rangeMax=rangeMin+range_mm
BY_increase=0.5

cbPalette <- c("#26456EFF", "#1C5F9EFF", "#4993C0FF", "#9EC9D9FF", "#D8C4B6FF", "#FE9C52FF","#FD8E3FFF")
Figure3_PanelB<-ggplot(tgca_combine, aes(x = PAR_anomaly_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(-2.6,rangeMax+0.1), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  scale_x_continuous(limits =c(-1.5,1.5), breaks=seq(-1.5,1.5,  by=0.5)) +
  geom_line(color="red") +
  geom_hline(yintercept = 0,  width=.6) +
  geom_smooth(aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") +
  scale_fill_manual(values=cbPalette)  +   
  scale_colour_manual(values=cbPalette)  +  theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
geom_point()
Figure3_PanelB
##--------------------------------------------------------------------------------------------------------------


##----------------------------------Figure 3 Panel C------------------------------------------------------------
WTD_set_array=((0:6)+0.1)*6.8 #set preferred intervals
Drought_Length_array=(0:60)/10#
number_group=61
WTD_inital<-0
byincrease<-8
for(i in seq(from=1, 6, by=1))
{
  Model_fit_array<-Drought_EVI_ModelPred_ClimateSensitivity(Drought_04De_ad_new.data, Drought_04De_ad_new.data_before_scaling,Ancova_model.lm.slope_model_climate,
                                                            WTD_set_array[i],x1.var='PAR_anomaly',
                                                            x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                                            x5.var='Drought_Length',x6.var='HAND',key_variable='Drought_Length',Drought_Length_array, moment=mean,number_group)
  newd0<-Model_fit_array
  if ((i-1)<10) {st_tag=paste(as.character(0),as.character(i-1))}else {
    st_tag=c(as.character(i-1))
  }
  st_tag<-gsub(" ", "", st_tag)
  print(st_tag)
  type_str<-paste('HAND',st_tag,as.character(WTD_inital+byincrease*(i-1.0)),sep = "_")
  newd0$type=type_str
  print(type_str)
  
  if (i == 1 ) {
    newd=rbind(newd0)
  }
  newd<-rbind(newd,newd0)
}

tgca_combine<-newd
BY_increase=0.5
rangeMin=-2.5#
range_mm=3.6#
rangeMax=rangeMin+range_mm
BY_increase=0.5
cbPalette <- c("#26456EFF", "#1C5F9EFF", "#4993C0FF", "#9EC9D9FF", "#D8C4B6FF", "#FE9C52FF","#FD8E3FFF")
Figure3_PanelC<-ggplot(tgca_combine, aes(x = Drought_Length_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(-2.6,rangeMax+0.1), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  scale_x_continuous(limits =c(0,6), breaks=seq(0,6,  by=1)) +
  geom_line(color="red") +
  geom_hline(yintercept = 0,  width=.6) +
  geom_smooth(aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") +
  scale_fill_manual(values=cbPalette)  +   
  scale_colour_manual(values=cbPalette)  +  theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
geom_point()
Figure3_PanelC

##--------------------------------------------------------------------------------------------------------------



##--------------------------------------Figure 3 Panel D--------------------------------------------------------
#Read Brienen related data files
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
Drought_year='2015'                                                                                                                       
file_name=paste('Brienen_Basin_validation','.csv',sep='')#  
file_ful_path=paste(file_path,"/",file_name,sep='')
Field.data<- read.csv(file=file_ful_path,header=T) 
Field_adjust.data<-Field.data

#calculate HAND class
Field_adjust.data$HAND_CLASS<-ceiling(Field_adjust.data$HAND_00_84) #HAND_00_84
Field_adjust.data[which(Field_adjust.data$HAND_CLASS >= 0),]$HAND_CLASS=(Field_adjust.data[which(Field_adjust.data$HAND_CLASS >= 0),]$HAND_00_84%/%2+1)*2
Field_adjust.data[which(Field_adjust.data$HAND_CLASS >= 0),]$HAND_CLASS=Field_adjust.data[which(Field_adjust.data$HAND_CLASS >= 0),]$HAND_CLASS-1
#calculate delta mortality for 2005
Field_adjust.data_Southern<-Field_adjust.data[which(Field_adjust.data$HAND_CLASS <= 40& Field_adjust.data$HAND_CLASS >0 & is.finite(Field_adjust.data$AGBMor_tot_2005)&
                                                     is.finite(Field_adjust.data$AGBMor_tot_mean) )  ,] #& 
Summary_Geo_2005<-Calculate_DeltaMortality(Field_adjust.data_Southern,x0.var='AGBMor_tot_2005',x1.var='AGBMor_tot_mean',x2.var='HAND_CLASS')
Summary_Geo_2005$combine_type="2005"
#calculate delta mortality for 2010
Field_adjust.data_Southern<-Field_adjust.data[which(Field_adjust.data$HAND_CLASS <= 40& Field_adjust.data$HAND_CLASS >0 & is.finite(Field_adjust.data$AGBMor_tot_2010)&
                                                     is.finite(Field_adjust.data$AGBMor_tot_mean) )  ,] 
Summary_Geo_2010<-Calculate_DeltaMortality(Field_adjust.data_Southern,x0.var='AGBMor_tot_2010',x1.var='AGBMor_tot_mean',x2.var='HAND_CLASS')
Summary_Geo_2010$combine_type="2010"
#calculate trend lines when HAND is less than 30 meter
Summary_Geo_2005_1<-Summary_Geo_2005[which(Summary_Geo_2005$HAND_CLASS<30),]
Summary_Geo_2010_1<-Summary_Geo_2010[which(Summary_Geo_2010$HAND_CLASS<30),]
lm.HD_Line_2005_ori_2 <-(lm(Delta_Mortality ~ HAND_CLASS,data=Summary_Geo_2005_1))
arr_ori_2_Line_2005<-coef(lm.HD_Line_2005_ori_2)
lm.HD_Line_2010_ori_2 <-(lm(Delta_Mortality ~ HAND_CLASS,data=Summary_Geo_2010_1))
arr_ori_2_Line_2010<-coef(lm.HD_Line_2010_ori_2)

#draw the plot
tgca_combine<-rbind(Summary_Geo_2005,Summary_Geo_2010)#
tgca_combine_sub<-tgca_combine[which(tgca_combine$HAND_CLASS >30),]
rangeMin=-0.5
rangeMax=rangeMin+1.7
BY_increase=0.5
size_number=0.7
cbPalette <- c( "#006837","#88419D", "#08519c")
Figure3_PanelD<-ggplot(tgca_combine, aes(x=HAND_CLASS, y=Delta_Mortality, colour = factor(combine_type),group=factor(combine_type),fill=factor(combine_type)))+ #, 
  geom_errorbar(aes(ymin=Delta_Mortality-se, ymax=Delta_Mortality+se), size=0.2, width=.1,face="bold") + #Drought_Condition
  geom_point( size=3.5, shape=23) + # 21is filled circle ,
  scale_fill_manual(values=cbPalette) + #
  geom_hline(yintercept = 0,width=.6, linetype="dashed") +
  geom_function(fun=Vectorize(function(x) {
    if(x > 30)
      return(NA)
    else
      return(arr_ori_2_Line_2005[2]*x+arr_ori_2_Line_2005[1])
  }),colour=c("#006837"),size=size_number) +
  
  geom_function(fun=Vectorize(function(x) {
    if(x > 30)
      return(NA)
    else
      return(arr_ori_2_Line_2010[2]*x+arr_ori_2_Line_2010[1])
  }),colour=c("#88419D"),size=size_number) +
  scale_colour_manual(values=cbPalette)  + 
  geom_point(data=tgca_combine_sub, aes(x=HAND_CLASS, y=Delta_Mortality), fill="white", shape=23, size=3.5, na.rm=TRUE)+
  xlab("HAND (meter)") +
  ylab("AGB mortality anomaly") +
  scale_y_continuous(limits =c(rangeMin-0.1,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  scale_x_continuous(limits =c(0,40),breaks=seq(0,40,10))+
  theme_bw() +
  theme( legend.position=c(1,0))+
  theme(legend.justification=c(1,0))+
  theme_few() %+replace% 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1))+
  theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
  theme(axis.title=element_text(size=24))+ #,face="bold"
  theme(plot.title=element_text(size=24))
Figure3_PanelD

#summary(lm.HD_Line_2005_ori_2)
#summary(lm.HD_Line_2010_ori_2)
##--------------------------------------------------------------------------------------------------------------

##--------------------------------------Figure 3 Panel E--------------------------------------------------------
#Read field data files
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
file_name=paste('Brienen_RAINFOR_Efforts_HAND','.csv',sep='')#  
file_ful_path=paste(file_path,"/",file_name,sep='')
Field.data<- read.csv(file=file_ful_path,header=T) 
Field_adjust.data<-Field.data
summary(Field_adjust.data)
Figure3_PanelE1<-ggplot(Field_adjust.data, aes(x=HAND, y=Effort_Ratio)) +
  geom_bar(stat="identity", color="#081D58", fill=c("#081D58"))+
  geom_hline(yintercept = 1,width=.6, linetype="dashed") +
  xlab("HAND (meter)") +
  ylab("Effort fraction:Area fraction") +
  theme_few() %+replace% 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1))+
  theme(axis.text.x = element_text( color="Black", size=20))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=20))+  #,face="bold"
  theme(axis.title=element_text(size=20))+ #,face="bold"
  theme(plot.title=element_text(size=20))
Figure3_PanelE1
##--------------------------------------Figure 3 Panel E continued----------------------------------------------
#Read field data files
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
file_name=paste('Brienen_RAINFOR_Efforts','.csv',sep='')#  
file_ful_path=paste(file_path,"/",file_name,sep='')
Field.data<- read.csv(file=file_ful_path,header=T) 
Field_adjust.data<-Field.data
summary(Field_adjust.data)
Figure3_PanelE2<-ggplot(Field_adjust.data, aes(x=HAND, y=Cumulative_HANDAreaProportion)) +
  geom_bar(stat="identity", colour="#969696", fill=c("#D9D9D9"))+
  scale_y_continuous(limits =c(0,1), breaks=seq(0,1,  by=0.2)) +
  xlab("HAND (meter)") +
  ylab("Area fraction") +
  theme_few() %+replace% 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1))+
  theme(axis.text.x = element_text( color="Black", size=14))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=20))+  #,face="bold"
  theme(axis.title=element_text(size=20))+ #,face="bold"
  theme(plot.title=element_text(size=20))
Figure3_PanelE2
##--------------------------------------------------------------------------------------------------------------


##-------------------------------------------Functions for Figure 3---------------------------------------------
#establish a selected GAM model
ancova_establish_slope_Guiana<- function(Input_model.data){
  names(Input_model.data)
  
  Input_model.data$SegGeo_numberf<-factor(Input_model.data$SegGeo_number)
  
  mod.IA<-mgcv::gam(EVI_anomaly ~ s(HAND, PAR_anomaly, k = 15, bs = "tp") + ti(HAND, Pre_anomaly, bs = "tp") + 
                      ti(HAND, PAR_anomaly, bs = "tp") + ti(HAND, VPD_anomaly, bs = "tp") + ti(HAND, Drought_Length, bs = "tp") +
                      s(VPD_anomaly, Drought_Length, bs = "tp") + s(MCWD_anomaly, Pre_anomaly, bs = "tp") +
                      ti(PAR_anomaly, Drought_Length, bs = "tp") + ti(PAR_anomaly, Pre_anomaly,  bs = "tp") + 
                      ti(VPD_anomaly, PAR_anomaly, bs = "tp") + ti(VPD_anomaly, Pre_anomaly,  bs = "tp")+ 
                      ti(HAND, MCWD_anomaly, bs = "tp")+ti(PAR_anomaly, MCWD_anomaly, bs = "tp") + 
                      ti(MCWD_anomaly, VPD_anomaly, bs = "tp"), method = "ML", data = Input_model.data) 

  
  return(mod.IA) 
}

#function to get mean/median values
set_value_byyear<- function(Data, x0.var='year',x1.var='PAR_anomaly',x2.var='VPD_anomaly',x3.var='MCWD_anomaly',
                          x4.var='Pre_anomaly',x5.var='Drought_Length',drought_year='2005',moment=median){
  ix0 <- match(x0.var, names(Data))
  ix1 <- match(x1.var, names(Data))
  ix2 <- match(x2.var, names(Data))
  ix3 <- match(x3.var, names(Data))
  ix4 <- match(x4.var, names(Data))
  ix5 <- match(x5.var, names(Data))
  
  print(ix0)
  Data[Data[,ix0]== drought_year,ix1]<-moment(Data[Data[,ix0]== drought_year,ix1]) # Severe Drought, Medium Drought,Modest Drought
  Data[Data[,ix0]== drought_year,ix2]<-moment(Data[Data[,ix0]== drought_year,ix2])
  Data[Data[,ix0]== drought_year,ix3]<-moment(Data[Data[,ix0]== drought_year,ix3])
  Data[Data[,ix0]== drought_year,ix4]<-moment(Data[Data[,ix0]== drought_year,ix4])
  Data[Data[,ix0]== drought_year,ix5]<-moment(Data[Data[,ix0]== drought_year,ix5])
  return(Data)
}

#function to correct EVI anomaly
Drought_EVI_correction<-function(Data1, Data2,x0.var='year',x1.var='EVI_anomaly',x2.var='HAND_CLASS',data1.var='EVI_anomaly_fit',
                                 data2.var='EVI_anomaly_fit_mean',drought_year='2015' ){
  D1_ix0 <- match(x0.var, names(Data1))
  D2_ix0 <- match(x0.var, names(Data2))
  Data11<-Data1[Data1[,D1_ix0] == drought_year, ]
  Data22<-Data2[Data2[,D2_ix0] == drought_year, ]
  
  D11_ix0 <- match(x0.var, names(Data11))
  D11_ix1 <- match(x1.var, names(Data11))
  D11_ix2 <- match(x2.var, names(Data11))

  D11_ix3 <- match(data1.var, names(Data11))
  D22_ix3 <- match(data2.var, names(Data22))
  
  correction<-Data11[,D11_ix3]-Data22[,D22_ix3]
  EVI_anomaly_corrected<-Data11[,D11_ix1] - correction#
  
  new_fit_array<-as.data.frame(cbind(Data11[,D11_ix1],correction,EVI_anomaly_corrected,Data11[,D11_ix2],Data11[,D11_ix0]))
  names(new_fit_array) <- c(x1.var,'correction','EVI_anomaly_corrected',x2.var,x0.var) 
  return(new_fit_array)
  
}

#calculate summary 
summary_group_full<-function(data,x0.var='EVI_anomaly',x1.var='EVI_anomaly_corrected', flag=1) # if using EVI_anomaly_corrected, flag=1 
{
  Summary_Geo=summary_group_simple(data,flag)
  Summary_Geo<-Summary_Geo[which(Summary_Geo$N >=4  & Summary_Geo$HAND_CLASS <=40  & Summary_Geo$HAND_CLASS >=2),]
  Summary_Geo<-data.frame(Summary_Geo$HAND_CLASS,Summary_Geo$EVI_anomaly_corrected,Summary_Geo$ci,Summary_Geo$se)
  names(Summary_Geo) <- c('HAND_CLASS',x0.var,'ci','se') 
  return(Summary_Geo)
}
#get summary information
summary_group_simple<-function(Input_model.data, group_flag){
  if (group_flag ==1) {
    Residual_summary <- summarySE(Input_model.data, measurevar="EVI_anomaly_corrected", groupvars=c("HAND_CLASS"))   #"Residual_R"
  } else{ if (group_flag ==0) {
    Residual_summary <- summarySE(Input_model.data, measurevar="EVI_anomaly", groupvars=c("HAND_CLASS")) 
  }else {
    Residual_summary <- summarySE(Input_model.data, measurevar="Prediction", groupvars=c("HAND_CLASS")) 
  }
  }
  return(Residual_summary)
}

#functions for Panel A prediction
Drought_EVI_ModelPrediction<-function(Data1, Data1_BeforeScale, model, WTD_array,x0.var='year',x1.var='PAR_anomaly',
                                      x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                      x5.var='Drought_Length',x6.var='HAND',drought_year='2005',moment=median,number_group){
  
  ix0 <- match(x0.var, names(Data1))
  
  data11<-Data1[Data1[,ix0] == drought_year, ]
  
  ix1 <- match(x1.var, names(data11))
  ix2 <- match(x2.var, names(data11))
  ix3 <- match(x3.var, names(data11))
  ix4 <- match(x4.var, names(data11))
  ix5 <- match(x5.var, names(data11))
  Dix6 <- match(x6.var, names(Data1_BeforeScale))
  
  x1.level <-moment(data11[,ix1]) # PAR_anomaly
  x2.level <-moment(data11[,ix2])
  x3.level <-moment(data11[,ix3])
  x4.level <-moment(data11[,ix4])
  x5.level <-moment(data11[,ix5])
  print(x2.level)
  x1.SD <-sd(data11[,ix1])
  
  new <- data.frame(WTD_array,rep(x1.level,number_group),+
                      rep(x2.level,number_group),+
                      rep(x3.level,number_group),+
                      rep(x4.level,number_group),+
                      rep(x5.level,number_group)) #,as.data.frame(WTD_array), Pre_anomaly=rep(pre_set,times=number_group) PAR_set_2005 -0.4204329
  names(new) <- c(x6.var,x1.var,x2.var,x3.var,x4.var,x5.var) 
  
  prediction_wtd00=predict(model,new,se.fit = TRUE)
  new$lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
  new$fit <- prediction_wtd00$fit
  new$uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
  nix6 <- match(x6.var, names(new))
  new$WTD_reverse=(new[,nix6]*sd(Data1_BeforeScale[, Dix6]))+mean(Data1_BeforeScale[, Dix6])
  names(new) <- c(x6.var,x1.var,x2.var,x3.var,x4.var,x5.var,'lci','fit','uci','WTD_reverse') 
  return(new)
}

#functions for merging two figures
merge_figures_sameunits<-function(base_plt,over_plt){
    
    plot_theme <- function(p) {
      plyr::defaults(p$theme, theme_get())
    }
    
    base_g = ggplot_gtable(ggplot_build(base_plt))
    overlay_g = ggplot_gtable(ggplot_build(over_plt))
    
    plt_panel = c(subset(base_g$layout, name == "panel", se = t:r))
    pnl_ind = which(overlay_g$layout$name == "panel")
    leg_ind = which(overlay_g$layout$name == "guide-box") 
    final_grob = gtable_add_grob(base_g,
                                 overlay_g$grobs[[pnl_ind]],
                                 plt_panel$t,
                                 plt_panel$l,
                                 plt_panel$b,
                                 plt_panel$r, name = "a")
    
    #final_grob = gtable_add_grob(final_grob,
                                 #overlay_g$grobs[[leg_ind]])
                                 #plt_panel$t,
                                 #plt_panel$l,
                                 #plt_panel$b,
                                 #plt_panel$r, name = "b") #
    return(final_grob)
  
}

#function for Panel B and C
Drought_EVI_ModelPred_ClimateSensitivity<-function(Data1, Data1_BeforeScale, model, WTD_value,x1.var='PAR_anomaly',
                                                   x2.var='VPD_anomaly',x3.var='MCWD_anomaly',x4.var='Pre_anomaly',
                                                   x5.var='Drought_Length',x6.var='HAND',key_variable='PAR_anomaly',key_var_array, moment=mean,number_group){
  
  
  ix1 <- match(x1.var, names(Data1))
  ix2 <- match(x2.var, names(Data1))
  ix3 <- match(x3.var, names(Data1))
  ix4 <- match(x4.var, names(Data1))
  ix5 <- match(x5.var, names(Data1))
  
  Dikey <- match(key_variable, names(Data1_BeforeScale))
  Dix6 <- match(x6.var, names(Data1_BeforeScale))
  
  x1.level <-moment(Data1[,ix1]) # PAR_anomaly
  x2.level <-moment(Data1[,ix2])
  x3.level <-moment(Data1[,ix3])
  x4.level <-moment(Data1[,ix4])
  x5.level <-moment(Data1[,ix5])
  print(x2.level)
  
  Dix6.SD <-sd(Data1_BeforeScale[,Dix6])
  WTD_value<-(WTD_value-mean(Data1_BeforeScale[,Dix6]))/Dix6.SD
  print(WTD_value)
  new <- data.frame(rep(WTD_value,number_group),rep(x1.level,number_group),+
                      rep(x2.level,number_group),+
                      rep(x3.level,number_group),+
                      rep(x4.level,number_group),+
                      rep(x5.level,number_group)) #,as.data.frame(WTD_array), Pre_anomaly=rep(pre_set,times=number_group) PAR_set_2005 -0.4204329
  
  names(new) <- c(x6.var,x1.var,x2.var,x3.var,x4.var,x5.var) 
  key_ix <- match(key_variable, names(new))
  new[,key_ix]<-(key_var_array-mean(Data1_BeforeScale[, Dikey]))/sd(Data1_BeforeScale[, Dikey]) # scale the key variable
  print( new[,key_ix])
  prediction_wtd00=predict(model,new,se.fit = TRUE)
  new$lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
  new$fit <- prediction_wtd00$fit
  new$uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
  new$key_variable_reverse=(new[,key_ix]*sd(Data1_BeforeScale[, Dikey]))+mean(Data1_BeforeScale[, Dikey])
  
  names(new) <- c(x6.var,x1.var,x2.var,x3.var,x4.var,x5.var,'lci','fit','uci',paste(key_variable,'_reverse',sep='')) 
  return(new)
}
#function for panel D (calculate the delta mortality)
Calculate_DeltaMortality<-function(Data,x0.var='AGBMor_tot_2005',x1.var='AGBMor_tot_mean',x2.var='HAND_CLASS'){
  ix0 <- match(x0.var, names(Data))
  ix1 <- match(x1.var, names(Data))
  ix2 <- match(x2.var, names(Data))
  
  Data$Delta_Mortality<-(Data[,ix0]-Data[,ix1])/Data[,ix1] #
  new_data<-as.data.frame(cbind(Data[,ix0],Data[,ix1],Data[,ix2],Data$Delta_Mortality))
  names(new_data) <- c(x0.var,x1.var,x2.var,'Delta_Mortality')
  
  new_data<-new_data[which(new_data[,4] != 0 | is.finite(new_data[,4]) )  ,]
  Summary_Geo<-summarySE(new_data, measurevar="Delta_Mortality", groupvars=c("HAND_CLASS")) 
  Summary_Geo<-Summary_Geo[which(Summary_Geo$N >1 ),]
  Summary_Geo<-data.frame(Summary_Geo$HAND_CLASS,Summary_Geo$Delta_Mortality,Summary_Geo$ci,Summary_Geo$se)
  names(Summary_Geo) <- c('HAND_CLASS','Delta_Mortality','ci','se')
  return(Summary_Geo)
}
##--------------------------------------------------------------------------------------------------------------




