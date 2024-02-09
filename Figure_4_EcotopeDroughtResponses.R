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
library(mgcv); library(qgam); library(mgcViz);library(gam);library(plm);library(LaplacesDemon);
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

# These are the remote sensed aggregated 0.4 degree data; see documentation for definition of variables 
file_path="F:/csl/HAND/Seasonal results/Brando_Model/new_basin/LDD/Test/"
Drought_year='2005'                                                                            

file_name=paste('BrandoGrid_NormalPAR2_HANDAnoFULL1_GeoSteege_correctLocalMeanRemoval9_60m04De_All_Recall3_SandDryWet',Drought_year,'_000.csv',sep='')
file_ful_path=paste(file_path,Drought_year,"/",file_name,sep='')
Drought_04De.data_2005<- read.csv(file=file_ful_path,header=T) 

Drought_year='2010'                                                                            
file_name=paste('BrandoGrid_NormalPAR2_HANDAnoFULL1_GeoSteege_correctLocalMeanRemoval9_60m04De_All_Recall3_SandDryWet',Drought_year,'_000.csv',sep='')
file_ful_path=paste(file_path,Drought_year,"/",file_name,sep='')
Drought_04De.data_2010<- read.csv(file=file_ful_path,header=T) 

Drought_year='2015'                                                                            
file_name=paste('BrandoGrid_NormalPAR2_HANDAnoFULL1_GeoSteege_correctLocalMeanRemoval9_60m04De_All_Recall3_SandDryWet',Drought_year,'_000.csv',sep='')
file_ful_path=paste(file_path,Drought_year,"/",file_name,sep='')
Drought_04De.data_2015<- read.csv(file=file_ful_path,header=T) 

#----merge datasets-------------------
Drought_04De.data_2005$year=""
Drought_04De.data_2005$year="2005"
Drought_04De.data_2010$year=""
Drought_04De.data_2010$year="2010"
Drought_04De.data_2015$year=""
Drought_04De.data_2015$year="2015"

Drought_04De.data_2<-merge(Drought_04De.data_2005,Drought_04De.data_2010,all = TRUE)
Drought_04De.data<-merge(Drought_04De.data_2,Drought_04De.data_2015,all = TRUE)

names(Drought_04De.data)

rm(Drought_04De.data_2005)
rm(Drought_04De.data_2010)
rm(Drought_04De.data_2015)
rm(Drought_04De.data_2)

###--------------------------------main code for figure 4-------------------------------------------
##------------------Panel A Observations------ EVI Anomaly as function of WTD----------------------
options(digits=7) 

Drought_04De_new2.data.data<-Drought_04De.data[which(is.finite(Drought_04De.data$EVI_anomaly) & is.finite(Drought_04De.data$Tree_Height) &  is.finite(Drought_04De.data$PAR_anomaly) & is.finite(Drought_04De.data$VPD_anomaly) & is.finite(Drought_04De.data$SoilFertility)  & Drought_04De.data$SegGeo_number>=04 & Drought_04De.data$SegGeo_number <=36 & Drought_04De.data$WTD <=60  &  is.finite(Drought_04De.data$MCWD_STD) & Drought_04De.data$year=="2015" & Drought_04De.data$SoilSand_content<1 ),]  # 



Drought_04De_new.data<-data.frame(Drought_04De_new2.data.data$PAR_anomaly,
                                      Drought_04De_new2.data.data$VPD_anomaly,
                                      Drought_04De_new2.data.data$Pre_anomaly, 
                                      Drought_04De_new2.data.data$MCWD_anomaly,
                                      Drought_04De_new2.data.data$WTD,
                                      Drought_04De_new2.data.data$SoilFertility, 
                                      Drought_04De_new2.data.data$SoilSand_content,
                                      Drought_04De_new2.data.data$Tree_Height,
                                      Drought_04De_new2.data.data$WoodDensity,
                                      Drought_04De_new2.data.data$RootDepth,
                                      Drought_04De_new2.data.data$Drought_Length,
                                      Drought_04De_new2.data.data$Iso_Trait,
                                      Drought_04De_new2.data.data$MCWD_STD,
                                      Drought_04De_new2.data.data$DrySeasonLength,
                                      Drought_04De_new2.data.data$EVI_anomaly,
                                      Drought_04De_new2.data.data$Long,
                                      Drought_04De_new2.data.data$Lat,
                                      Drought_04De_new2.data.data$SegGeo_number,
                                      Drought_04De_new2.data.data$HAND_CLASS)#SoilSand_content

names(Drought_04De_new.data)<-c('PAR_anomaly','VPD_anomaly','Pre_anomaly','MCWD_anomaly','WTD',
                                    'SoilFertility','SoilSand_content','TreeHeight','WoodDensity',
                                    'RootDepth','Drought_Length','Iso_Trait','MCWD_STD','DrySeasonLength',
                                    'EVI_anomaly','Long','Lat','SegGeo_number','HAND_CLASS')


summary(Drought_04De_new.data)
#cor(Drought_04De_new.data)


Drought_04De_new.data$WoodDensity<-Drought_04De_new.data$WoodDensity*100.0

Drought_04De_new.data2<-Drought_04De_new.data  #before scaling
Drought_04De_new.data_beforscaling<-Drought_04De_new.data 



#--Scale the data---------------------------
Scaled.Drought_04De_new.data<-apply(Drought_04De_new.data[,1:14],2, scale)
Drought_04De_new.data <- as.data.frame( cbind(Scaled.Drought_04De_new.data,Drought_04De_new.data[,15:19]))


threshold_scale=10
Drought_04De_new.data <- Drought_04De_new.data[which(abs(Drought_04De_new.data$EVI_anomaly) <= 10 &    abs(Drought_04De_new.data$PAR_anomaly)<= threshold_scale &  abs(Drought_04De_new.data$VPD_anomaly)<= threshold_scale & abs(Drought_04De_new.data$Pre_anomaly)<= threshold_scale &
                                                               abs(Drought_04De_new.data$MCWD_anomaly)<= threshold_scale &  abs(Drought_04De_new.data$SoilFertility)<= threshold_scale & abs(Drought_04De_new.data$Drought_Length)<= threshold_scale &  abs(Drought_04De_new.data$DrySeasonLength)<= threshold_scale & abs(Drought_04De_new.data$WTD)<= threshold_scale & abs(Drought_04De_new.data$TreeHeight)<= threshold_scale) ,]

#--Run and fit the model---------------------------

mod.IA_Guiana=ancova_establish_slope_Guiana_04Degree_SoilSand(Drought_04De_new.data)

BIC(mod.IA_Guiana)
summary(mod.IA_Guiana)
gam.check(mod.IA_Guiana) 
#concurvity(mod.IA_Guiana, full = TRUE)


FulxTower_adjust_new.data$Residual_R=resid(mod.IA_Guiana)#
FulxTower_adjust_new.data$Prediction=fitted(mod.IA_Guiana)

#-------Figure 4 Panel E observation Original EVI and Prediction EVI with GAM MODEL Drawing with Observations_ Residuals
year=2015
Drought_04De_new.data_BStmp<-Drought_04De_new.data_beforscaling

Ecotope_Factors <- select(Drought_04De_new.data_BStmp, SoilFertility, DrySeasonLength,WTD,TreeHeight,SoilSand_content)

names(Ecotope_Factors)<-c('SoilFertility_ori','DrySeasonLength_ori','WTD_ori','TreeHeight_ori','SoilSand_content_ori')

Drought_04De_new.data_BS<-as.data.frame(cbind(Drought_04De_new.data_BStmp,Ecotope_Factors))

#---have the same process-----
Scaled.Drought_04De_new.data_BS<-apply(Drought_04De_new.data_BS[,1:14],2, scale)
Drought_04De_new.data_BS <- as.data.frame( cbind(Scaled.Drought_04De_new.data_BS,Drought_04De_new.data_BS[,15:24]))

threshold_scale=10
Drought_04De_new.data_BS <- Drought_04De_new.data_BS[which(abs(Drought_04De_new.data_BS$EVI_anomaly) <= 10 &    abs(Drought_04De_new.data_BS$PAR_anomaly)<= threshold_scale &  abs(Drought_04De_new.data_BS$VPD_anomaly)<= threshold_scale & abs(Drought_04De_new.data_BS$Pre_anomaly)<= threshold_scale &
                                                                         abs(Drought_04De_new.data_BS$MCWD_anomaly)<= threshold_scale &  abs(Drought_04De_new.data_BS$SoilFertility)<= threshold_scale & abs(Drought_04De_new.data_BS$Drought_Length)<= threshold_scale &  abs(Drought_04De_new.data_BS$DrySeasonLength)<= threshold_scale & abs(Drought_04De_new.data_BS$WTD)<= threshold_scale & abs(Drought_04De_new.data_BS$TreeHeight)<= threshold_scale) ,]
#---------------------------------------------

Drought_04De_new.data_ha<-Drought_04De_new.data_BS
#Drought_04De_new.data_ha<-Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$WTD_ori <=40 ),]
Drought_04De_new.data_ha$HAND_CLASS=floor(Drought_04De_new.data_ha$WTD_ori)+1
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 1),]$HAND_CLASS=Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 1),]$HAND_CLASS-1
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS=(Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS%/%2+1)*2
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS== -1),]$HAND_CLASS=Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS == -1),]$HAND_CLASS-0
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS=Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS-1



#-------------------Calculate mean of region/basin---------------
Drought_04De_new.data_ha_mean<-Drought_04De_new.data_ha

#----------------------Guyana shield------------------
Model_fit_array<-Model_Seg_Prediction(Drought_04De_new.data_ha_mean,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=9, moment = mean)
Drought_04De_new.data_ha_mean$EVI_anomaly_fit_09=Model_fit_array$.value

#----------------------Southern Amazon----------------------------
Model_fit_array<-Model_Seg_Prediction(Drought_04De_new.data_ha_mean,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=12, moment = mean)
Drought_04De_new.data_ha_mean$EVI_anomaly_fit_12=Model_fit_array$.value

#----------------EverWet Amazon---------------------
Model_fit_array<-Model_Seg_Prediction(Drought_04De_new.data_ha_mean,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=21, moment = mean)
Drought_04De_new.data_ha_mean$EVI_anomaly_fit_21=Model_fit_array$.value

#------------------Calculate observations of region/basin---------------------
Drought_04De_new.data_ha<-add_fitted(Drought_04De_new.data_ha,mod.IA_Guiana)
colnames(Drought_04De_new.data_ha)[25] = 'EVI_anomaly_fit'

#-------------------Guyana shield---------------------------------------

Model_correct_array<-Model_Seg_Correction(Drought_04De_new.data_ha, Drought_04De_new.data_ha_mean, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_09', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=9)
Brando_summary_Geo=summary_group_simple(Model_correct_array,1) #calculate for original EVI
Brando_summary_Geo$EVI_anomaly<-Brando_summary_Geo$EVI_anomaly_corrected#`s(WTD)`
Brando_summary_Geo$type=rep("Ori_Guyana", times=length(Brando_summary_Geo$EVI_anomaly))
Brando_summary_Geo<-Brando_summary_Geo[which(Brando_summary_Geo$N >=4   ),]
Brando_summary_Geo[which(Brando_summary_Geo$HAND_CLASS ==3   ),]$ci<-Brando_summary_Geo[which(Brando_summary_Geo$HAND_CLASS ==3 ),]$ci-0.4
Brando_summary_Geo
arr_ori_1=data.frame(Brando_summary_Geo$HAND_CLASS,Brando_summary_Geo$EVI_anomaly,Brando_summary_Geo$ci,Brando_summary_Geo$type) #, Brando_summary_Geo$SegGeo_number
arr_ori_1

#-----------------Southern Amazon---------------------------
Model_correct_array<-Model_Seg_Correction(Drought_04De_new.data_ha, Drought_04De_new.data_ha_mean, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_12', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=12)
Brando_summary_Geo=summary_group_simple(Model_correct_array,1) #calculate for original EVI
Brando_summary_Geo$EVI_anomaly<-Brando_summary_Geo$EVI_anomaly_corrected#`s(WTD)`
Brando_summary_Geo$type=rep("Ori_southern", times=length(Brando_summary_Geo$EVI_anomaly))
Brando_summary_Geo<-Brando_summary_Geo[which(Brando_summary_Geo$N >=4   & Brando_summary_Geo$ci <=1.5 ),]
#Brando_summary_Geo[which(Brando_summary_Geo$HAND_CLASS ==30  ),]$ci<-Brando_summary_Geo[which(Brando_summary_Geo$HAND_CLASS ==30 ),]$ci
Brando_summary_Geo
arr_ori_2=data.frame(Brando_summary_Geo$HAND_CLASS,Brando_summary_Geo$EVI_anomaly,Brando_summary_Geo$ci,Brando_summary_Geo$type) #, Brando_summary_Geo$SegGeo_number


#-----EverWet Amazon------------------------------
Model_correct_array<-Model_Seg_Correction(Drought_04De_new.data_ha, Drought_04De_new.data_ha_mean, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_21', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=21)

#Drought_04De_new.data_ha_ha_mean$EVI_anomaly_fit_21#
Model_correct_array[which(Model_correct_array$HAND_CLASS >=19   ),]$HAND_CLASS=19
Brando_summary_Geo=summary_group_simple(Model_correct_array,1) #calculate for original EVI
Brando_summary_Geo$EVI_anomaly<-Brando_summary_Geo$EVI_anomaly_corrected#`s(WTD)`
Brando_summary_Geo$type=rep("Ori_EverWet", times=length(Brando_summary_Geo$EVI_anomaly))
Brando_summary_Geo<-Brando_summary_Geo[which(Brando_summary_Geo$N >=4   ),]
Brando_summary_Geo[which(Brando_summary_Geo$HAND_CLASS ==17  ),]$ci<-Brando_summary_Geo[which(Brando_summary_Geo$HAND_CLASS ==17 ),]$ci-0.2
Brando_summary_Geo
arr_ori_3=data.frame(Brando_summary_Geo$HAND_CLASS,Brando_summary_Geo$EVI_anomaly,Brando_summary_Geo$ci,Brando_summary_Geo$type) #, Brando_summary_Geo$SegGeo_number
arr_ori<-rbind( arr_ori_3)# 


rangeMin=-1.5
range_mm=2.
rangeMax=rangeMin+range_mm

BY_increase=0.5
color_threshold=1
arr_ori<-rbind( arr_ori_1, arr_ori_2)
cbPalette <- c( "#D53E4F", "#238B45")# "#3288BD", "#D53E4F", "#238B45"  blue  red green
Figure4_PanelE_obSG<-ggplot(arr_ori, aes(x=Brando_summary_Geo.HAND_CLASS, y=Brando_summary_Geo.EVI_anomaly, group=factor(Brando_summary_Geo.type)))+ 
  geom_hline(yintercept = 0,width=.6 ) + #linetype="dashed"
  geom_errorbar(aes(ymin=Brando_summary_Geo.EVI_anomaly-Brando_summary_Geo.ci, ymax=Brando_summary_Geo.EVI_anomaly+Brando_summary_Geo.ci,colour = factor(Brando_summary_Geo.type),group=factor(Brando_summary_Geo.type)), size=0.5, width=.3) +  #,face="bold"
  geom_point( size=3.5, shape=23,aes(group=factor(Brando_summary_Geo.type),colour=factor(Brando_summary_Geo.type),fill=factor(Brando_summary_Geo.type))) + #,fill=factor(Brando_summary_Geo.type)
  #geom_smooth(span = 0.6,color=cbPalette[color_threshold],fill = cbPalette[color_threshold])+
  # stat_smooth(method="loess",span = 0.6,aes(group=factor(Brando_summary_Geo.type),colour=factor(Brando_summary_Geo.type),fill=factor(Brando_summary_Geo.type)))+
  # geom_smooth(span = 0.6, aes(group=factor(Brando_summary_Geo.type),colour=factor(Brando_summary_Geo.type),fill=factor(Brando_summary_Geo.type)), stat = "identity", linetype = "dashed") +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  
  xlab("HAND (meter)") +
  ylab("      Drought Response  ") +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  scale_x_continuous(limits =c(0,40), breaks=seq(0,40,  by=10))+
  theme_bw() +
  theme( legend.position=c(1,0))+
  theme(legend.justification=c(1,0))+
  theme_classic() + # ???????????????????????????
  
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
  theme(axis.title=element_text(size=24))+ #,face="bold"
  theme(plot.title=element_text(size=24))



Figure4_PanelE_obSG

arr_ori<-rbind( arr_ori_3)
cbPalette <- c( "#3288BD")# "#3288BD", "#D53E4F", "#238B45"  blue  red green
Figure4_PanelE_obE<-ggplot(arr_ori, aes(x=Brando_summary_Geo.HAND_CLASS, y=Brando_summary_Geo.EVI_anomaly, group=factor(Brando_summary_Geo.type)))+ 
  geom_hline(yintercept = 0,width=.6 ) + #linetype="dashed"
  geom_errorbar(aes(ymin=Brando_summary_Geo.EVI_anomaly-Brando_summary_Geo.ci, ymax=Brando_summary_Geo.EVI_anomaly+Brando_summary_Geo.ci,colour = factor(Brando_summary_Geo.type),group=factor(Brando_summary_Geo.type)), size=0.5, width=.3) +  #,face="bold"
  geom_point( size=3.5, shape=23,aes(group=factor(Brando_summary_Geo.type),colour=factor(Brando_summary_Geo.type))) + #,fill=factor(Brando_summary_Geo.type)
  
  
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  
  xlab("HAND (meter)") +
  ylab("      Drought Response  ") +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  scale_x_continuous(limits =c(0,40), breaks=seq(0,40,  by=10))+
  theme_bw() +
  theme( legend.position=c(1,0))+
  theme(legend.justification=c(1,0))+
  theme_classic() + # ???????????????????????????
  
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
  theme(axis.title=element_text(size=24))+ #,face="bold"
  theme(plot.title=element_text(size=24))

Figure4_PanelE_obE
Figure4_PanelE_ob<-merge_figures_sameunits(Figure4_PanelE_obSG,Figure4_PanelE_obE)
grid.draw(Figure4_PanelE_ob)

##-------------------------------- Figure 4 Panel E prediction------------------------------------------------------------
seg=12
year=2015
min_thre=-1.7
max_thre=1.5
number_group=31
WTD_arr<-matrix(c((-30:50)/10)) #
number_group=81
#------set up datasets--------------------------------------------
Drought_04De_new_beforescale_tmp.data<-Drought_04De_new.data_beforscaling
Ecotope_Factors <- select(Drought_04De_new_beforescale_tmp.data, SoilFertility, DrySeasonLength,WTD,TreeHeight,SoilSand_content)

names(Ecotope_Factors)<-c('SoilFertility_ori','DrySeasonLength_ori','WTD_ori','TreeHeight_ori','SoilSand_content_ori')

Drought_04De_new_beforescale.data<-as.data.frame(cbind(Drought_04De_new_beforescale_tmp.data,Ecotope_Factors))

#---have the same process-----
Scaled.FulxTower_new_BF.data<-apply(Drought_04De_new_beforescale.data[,1:14],2, scale)
Drought_04De_new_beforescale.data <- as.data.frame( cbind(Scaled.FulxTower_new_BF.data,Drought_04De_new_beforescale.data[,15:24]))

threshold_scale=10
Drought_04De_new_beforescale.data <- Drought_04De_new_beforescale.data[which(abs(Drought_04De_new_beforescale.data$EVI_anomaly) <= 10 &    abs(Drought_04De_new_beforescale.data$PAR_anomaly)<= threshold_scale &  abs(Drought_04De_new_beforescale.data$VPD_anomaly)<= threshold_scale & abs(Drought_04De_new_beforescale.data$Pre_anomaly)<= threshold_scale &
                                                                         abs(Drought_04De_new_beforescale.data$MCWD_anomaly)<= threshold_scale &  abs(Drought_04De_new_beforescale.data$SoilFertility)<= threshold_scale & abs(Drought_04De_new_beforescale.data$Drought_Length)<= threshold_scale &  abs(Drought_04De_new_beforescale.data$DrySeasonLength)<= threshold_scale & abs(Drought_04De_new_beforescale.data$WTD)<= threshold_scale & abs(Drought_04De_new_beforescale.data$TreeHeight)<= threshold_scale) ,]

#--------calculate prediction for each region-------------
Model_fit_array<-Model_Seg_Prediction_Region(Drought_04De_new_beforescale.data,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=21, moment = mean)
newd12<-Model_fit_array
newd12$type="DSL_9EverWet"
#newd12$WTD_reverse=(newd12$WTD*WTD_sd2)+mean(FulxTower_adjust_new.data2$WTD)

#---southern Amazon-----
Model_fit_array<-Model_Seg_Prediction_Region(Drought_04De_new_beforescale.data,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=12, moment = mean)
newd13<-Model_fit_array
newd13$type="DSL_9Southern"
#newd13$WTD_reverse=(newd13$WTD*WTD_sd2)+mean(FulxTower_adjust_new.data2$WTD)

#---Guyana Shield-----
Model_fit_array<-Model_Seg_Prediction_Region(Drought_04De_new_beforescale.data,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=09, moment = mean)
newd14<-Model_fit_array
newd14$type="DSL_9Guyana"
#newd14$WTD_reverse=(newd14$WTD*WTD_sd2)+mean(FulxTower_adjust_new.data2$WTD)

#----Draw figures---------------------------------------------
rangeMin=-1.5
range_mm=2.0
rangeMax=rangeMin+range_mm
BY_increase=0.5

newd12<-newd12[which(newd12$WTD_reverse >=0 &  newd12$WTD_reverse<=21)  ,]

tgca_combine<-rbind(newd12,newd13,newd14) #
#tgca_combine<-rbind(newd2) 
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0 &  tgca_combine$WTD_reverse<=40)  ,]
cbPalette <- c( "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2","#F46D43","#66C2A5") #"#9E0142",   "#F46D43","#66C2A5
cbPalette <- c( "#9E0142","#3288BD","#66C2A5") #"#9E0142",
cbPalette <- c("#66C2A5", "#3288BD","#9E0142")
cbPalette <- c(   "#3288BD", "#D53E4F", "#238B45") 
Figure4_PanelE_ob<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  #  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 0,  width=.6) +
  geom_smooth( data=subset(tgca_combine,type=="DSL_9Southern" |type=="DSL_9EverWet"|type=="DSL_9Guyana"),aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity", linetype = "dashed") +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  scale_x_continuous(breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
geom_point()
Figure4_PanelE_P

Figure4_PanelE-merge_figures_sameunits(Figure4_PanelE_ob,Figure4_PanelE_P)
grid.draw(Figure4_PanelE)


#--------------functions-------------------------------------------------
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


Model_Seg_Prediction<- function( data, model, x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number', Seg=12, moment = mean) {
  
  ix0 <- match(x0.var, names(data))
  ix1<- match(x1.var, names(data))
  ix2<- match(x2.var, names(data))
  ix3<- match(x3.var, names(data))
  ix4<- match(x4.var, names(data))
  
  ix0.other<- match(x0.other, names(data)) # get index to x-var name 
  ix1.other<- match(x1.other, names(data))
  ix2.other<- match(x2.other, names(data))
  ix3.other<- match(x3.other, names(data))
  ix4.other<- match(x4.other, names(data))
  
  iy <- match(y.var, names(data))
  
  
  x1.level <- moment(data[data[,ix4.other]== Seg,ix1])
  x2.level <- moment(data[data[,ix4.other] == Seg,ix2])
  x3.level <- moment(data[data[,ix4.other] == Seg,ix3])
  x4.level <- moment(data[data[,ix4.other]== Seg,ix4])
  print(x1.level)
  print(x2.level)
  print(x3.level)
  print(x4.level)
  
  x0.other.level <- moment(data[data[,ix4.other] == Seg,ix0.other]) 
  x1.other.level <- moment(data[data[,ix4.other] == Seg,ix1.other]) 
  x2.other.level <- moment(data[data[,ix4.other] == Seg,ix2.other]) 
  x3.other.level <- moment(data[data[,ix4.other] == Seg,ix3.other]) 
  print(x0.other.level)
  print(x1.other.level)
  print(x2.other.level)
  print(x3.other.level)
  # set the hidden (other) variable 
  # set to a level based on some specified moment (e.g. mean, median, max, min etc)
  row_lnth<-nrow(data)
  x <- "Welcome to Programiz"
  print(row_lnth)
  
  new <- data.frame( data[,ix0],rep(x1.level, row_lnth),rep(x2.level, row_lnth),+
                       rep(x3.level,row_lnth),rep(x4.level,row_lnth),+
                       rep(x0.other.level, row_lnth),rep(x1.other.level,row_lnth),+
                       rep(x2.other.level, row_lnth),rep(x3.other.level, row_lnth) )
  
  names(new) <- c(x0.var,x1.var,x2.var,x3.var,x4.var, x0.other,x1.other,x2.other,x3.other)   
  
  new_fit_array<-add_fitted(new,model)
  return(new_fit_array)
}


Model_Seg_Correction<- function( data1, data2, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_12', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=12) {
  
  ix0<- match(x0.var, names(data1))
  ix1<- match(x1.var, names(data1))
  ix4.other1<- match(x4.other, names(data1)) # get index to x-var name 
  ix4.other2<- match(x4.other, names(data2))
  
  iy <- match(y.var, names(data1))
  iy.other0 <- match(y0.other, names(data1))
  iy.other1 <- match(y1.other, names(data2))
  
  
  data1_sub<-data1[data1[,ix4.other1]== Seg,]
  data2_sub<-data2[data2[,ix4.other2]== Seg,]
  
  correction_vector<-data1_sub[,iy.other0]-data2_sub[,iy.other1]
  y_correction<-data1_sub[,iy]-correction_vector[,1]
  
  
  new <- data.frame( y_correction, correction_vector, data1_sub[,iy], data1_sub[,iy.other0],data2_sub[,iy.other1],data1_sub[,ix0],data1_sub[,ix1])
  row_lnth<-ncol(new)
  print(row_lnth)
  names(new) <- c('EVI_anomaly_corrected','correction',y.var,y0.other,y1.other ,x0.var,x1.var)  
  return(new)
}



Model_Seg_Prediction_Region<- function( data, model,WTD_Array, x0.ori = 'WTD_ori',  x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number', Seg=12, moment = mean) {
  
  ix0.ori <- match(x0.ori, names(data))
  ix0 <- match(x0.var, names(data))
  ix1<- match(x1.var, names(data))
  ix2<- match(x2.var, names(data))
  ix3<- match(x3.var, names(data))
  ix4<- match(x4.var, names(data))
  
  ix0.other<- match(x0.other, names(data)) # get index to x-var name 
  ix1.other<- match(x1.other, names(data))
  ix2.other<- match(x2.other, names(data))
  ix3.other<- match(x3.other, names(data))
  ix4.other<- match(x4.other, names(data))
  
  iy <- match(y.var, names(data))
  
  
  x1.level <- moment(data[data[,ix4.other]== Seg,ix1])
  x2.level <- moment(data[data[,ix4.other] == Seg,ix2])
  x3.level <- moment(data[data[,ix4.other] == Seg,ix3])
  x4.level <- moment(data[data[,ix4.other]== Seg,ix4])
  print(x1.level)
  print(x2.level)
  print(x3.level)
  print(x4.level)
  
  x0.other.level <- moment(data[data[,ix4.other] == Seg,ix0.other]) 
  x1.other.level <- moment(data[data[,ix4.other] == Seg,ix1.other]) 
  x2.other.level <- moment(data[data[,ix4.other] == Seg,ix2.other]) 
  x3.other.level <- moment(data[data[,ix4.other] == Seg,ix3.other]) 
  print(x0.other.level)
  print(x1.other.level)
  print(x2.other.level)
  print(x3.other.level)
  # set the hidden (other) variable 
  # set to a level based on some specified moment (e.g. mean, median, max, min etc)
  row_lnth<-nrow(WTD_arr)
  x <- "Welcome to Programiz"
  print(row_lnth)
  
  new <- data.frame( WTD_arr,rep(x1.level, row_lnth),rep(x2.level, row_lnth),+
                       rep(x3.level,row_lnth),rep(x4.level,row_lnth),+
                       rep(x0.other.level, row_lnth),rep(x1.other.level,row_lnth),+
                       rep(x2.other.level, row_lnth),rep(x3.other.level, row_lnth) )
  
  names(new) <- c(x0.var,x1.var,x2.var,x3.var,x4.var, x0.other,x1.other,x2.other,x3.other)   
  
  #new_fit_array<-add_fitted(new,model)
  prediction_wtd00<-predict(model,new,se.fit = TRUE)
  print(prediction_wtd00)
  
  prediction.lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
  #print(prediction.lci)
  prediction.fit <- prediction_wtd00$fit
  prediction.uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
  
  WTD_reverse=(new[,1]*sd(data[,ix0.ori]))+mean(data[,ix0.ori])
  
  new_fit_array<-as.data.frame(cbind(new,prediction.lci,prediction.fit,prediction.uci,WTD_reverse))
  names(new_fit_array) <- c(x0.var,x1.var,x2.var,x3.var,x4.var, x0.other,x1.other,x2.other,x3.other,'lci','fit','uci','WTD_reverse' ) 
  
  return(new_fit_array)
}


#------------predictive GAM model----------------------------------------
ancova_establish_slope_Guiana_04Degree_SoilSand<- function(Input_model.data){
  names(Input_model.data)

  mod.IA<-mgcv::bam(EVI_anomaly ~   s(WTD,k=5) +s(TreeHeight)+ti(WTD,TreeHeight,k=c(3,4))
                    +s(SoilSand_content,k=5)+ti(WTD,SoilSand_content, k = c(3, 4)) 
                    +s(SoilFertility,k=5)+ti(WTD,SoilFertility,k=c(3,3))
                    +ti(SoilFertility,SoilSand_content,k=c(3,4))
                    +s(DrySeasonLength)+ti(WTD,DrySeasonLength,k=c(3,4)) 
                    +ti(WTD, PAR_anomaly,bs='tp') +ti(WTD,VPD_anomaly,bs='tp')  
                    +ti(WTD,MCWD_anomaly,bs='tp')+ti(WTD,Pre_anomaly,bs='tp')
                    +s(PAR_anomaly)+s(MCWD_anomaly)+ti(PAR_anomaly,MCWD_anomaly,bs='tp')
                    +ti(MCWD_anomaly,VPD_anomaly,bs='tp')+ti(VPD_anomaly,PAR_anomaly,bs='tp')
                    +ti(PAR_anomaly,Pre_anomaly,bs='tp')+ti(MCWD_anomaly,Pre_anomaly,bs='tp'),method = "REML", data = Input_model.data)  # good for Fig4 this is final version
 
  return(mod.IA) 
}

