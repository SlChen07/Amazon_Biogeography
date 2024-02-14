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

# These are the remote sensed aggregated 0.4 degree data; see documentation for definition of variables 
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
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

#---------merge datasets-------------------
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


#-------------------------------------Figure 5 Panel A  Resilience-------------------------------------------
# Biogeography of drought --- simplied  whole basin 2015 drought _04Degree
##-------------- Read Ecotope data files----------------------------------------------------------------------
gc()
file_path="F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests/data/"
Drought_year='2005'                                                                            
file_name=paste('BrandoGrid_NormalPAR2_HANDAnoFULL1_GeoSteege_correctLocalMeanRemoval9_60m04De_All_H_BiogeographyLandcoverSoilSand',Drought_year,'_000.csv',sep='') 
file_ful_path=paste(file_path,Drought_year,"/",file_name,sep='')
Ecotope.data_2005<- read.csv(file=file_ful_path,header=T) 

Drought_year='2010'                                                                            
file_name=paste('BrandoGrid_NormalPAR2_HANDAnoFULL1_GeoSteege_correctLocalMeanRemoval9_60m04De_All_H_BiogeographyLandcoverSoilSand',Drought_year,'_000.csv',sep='') 
file_ful_path=paste(file_path,Drought_year,"/",file_name,sep='')
Ecotope.data_2010<- read.csv(file=file_ful_path,header=T) 

Drought_year='2015'                                                                            
file_name=paste('BrandoGrid_NormalPAR2_HANDAnoFULL1_GeoSteege_correctLocalMeanRemoval9_60m04De_All_H_BiogeographyLandcoverSoilSand',Drought_year,'_000.csv',sep='') 
file_ful_path=paste(file_path,Drought_year,"/",file_name,sep='')
Ecotope.data_2015<- read.csv(file=file_ful_path,header=T) 

Ecotope.data_2005$year=""
Ecotope.data_2005$year="2005"
Ecotope.data_2010$year=""
Ecotope.data_2010$year="2010"
Ecotope.data_2015$year=""
Ecotope.data_2015$year="2015"

Ecotope.data_2<-merge(Ecotope.data_2005,Ecotope.data_2010,all = TRUE)
Ecotope.data<-merge(Ecotope.data_2,Ecotope.data_2015,all = TRUE)

#names(Ecotope.data)

rm(Ecotope.data_2005)
rm(Ecotope.data_2010)
rm(Ecotope.data_2015)
rm(Ecotope.data_2)



#---calculate the values from model data-------
Ecotope.data_BStmp<-Ecotope.data


#----calculate the whole basin for all evergreen forest--
Ecotope.data_BS_ad<-Ecotope.data_BStmp
Ecotope.data_BS_ad<-Ecotope.data_BS_ad[which(Ecotope.data_BS_ad$SegGeo_number >= 0 & Ecotope.data_BS_ad$SegGeo_number <=36  & Ecotope.data_BS_ad$WTD <=10000),]

Drought_04De_new2.data_BS_ad<-Ecotope.data_BS_ad[which(  is.finite(Ecotope.data_BS_ad$SoilFertility) &  is.finite(Ecotope.data_BS_ad$Tree_Height) &  is.finite(Ecotope.data_BS_ad$SoilSand_content) & Ecotope.data_BS_ad$year =="2015" &  is.finite(Ecotope.data_BS_ad$WTD) & Ecotope.data_BS_ad$WTD <=100 & Ecotope.data_BS_ad$SoilSand_content <=0.9 & Ecotope.data_BS_ad$Long !=0  & Ecotope.data_BS_ad$Lat !=0 ),] # 

Ecotope.data_BS_ad<-data.frame(Drought_04De_new2.data_BS_ad$PAR_anomaly,Drought_04De_new2.data_BS_ad$VPD_anomaly,Drought_04De_new2.data_BS_ad$Pre_anomaly, Drought_04De_new2.data_BS_ad$MCWD_anomaly,Drought_04De_new2.data_BS_ad$WTD,Drought_04De_new2.data_BS_ad$WTD,Drought_04De_new2.data_BS_ad$SoilFertility,Drought_04De_new2.data_BS_ad$SoilSand_content,Drought_04De_new2.data_BS_ad$Tree_Height,Drought_04De_new2.data_BS_ad$MCWD_STD,Drought_04De_new2.data_BS_ad$DrySeasonLength,Drought_04De_new2.data_BS_ad$EVI_anomaly,Drought_04De_new2.data_BS_ad$Long,Drought_04De_new2.data_BS_ad$Lat,Drought_04De_new2.data_BS_ad$SegGeo_number,Drought_04De_new2.data_BS_ad$HAND_CLASS,Drought_04De_new2.data_BS_ad$year)#SoilSand_content

names(Ecotope.data_BS_ad)<-c('PAR_anomaly','VPD_anomaly','Pre_anomaly','MCWD_anomaly','WTD','HAND','SoilFertility','SoilSand_content','TreeHeight','MCWD_STD','DrySeasonLength','EVI_anomaly','Long','Lat','SegGeo_number','HAND_CLASS','year')




Ecotope.data_BS_ad_ha<-Ecotope.data_BS_ad

#---------------------------------------------
Ecotope.data_BS_ad_ha$HAND_CLASS<-floor(Ecotope.data_BS_ad_ha$HAND)+1 ## 0.4 Degree

Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS >= 1),]$HAND_CLASS=Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS >= 1),]$HAND_CLASS-1
Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS >= 0),]$HAND_CLASS=(Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS >= 0),]$HAND_CLASS%/%2+1)*2
Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS== -1),]$HAND_CLASS=Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS == -1),]$HAND_CLASS-0
Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS >= 0),]$HAND_CLASS=Ecotope.data_BS_ad_ha[which(Ecotope.data_BS_ad_ha$HAND_CLASS >= 0),]$HAND_CLASS-1



#-------------------Calculate mean of region/basin---------------
Ecotope.data_BS_ad_ha_mean<-Ecotope.data_BS_ad_ha

Model_fit_array<-Model_Seg_Resilience(Ecotope.data_BS_ad_ha_mean,Drought_04De_new.data_beforscaling,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', moment = mean)

Ecotope.data_BS_ad_ha_mean$EVI_anomaly_fit_mean<- Model_fit_array$.value
Ecotope.data_BS_ad_ha_mean$EVI_anomaly_fit_standardize<- (Ecotope.data_BS_ad_ha_mean$EVI_anomaly_fit_mean-mean(Ecotope.data_BS_ad_ha_mean$EVI_anomaly_fit_mean))/sd(Ecotope.data_BS_ad_ha_mean$EVI_anomaly_fit_mean)

currentTime <-Sys.Date()

file_ful_path=paste(file_path,Drought_year,"/",currentTime,"Seg_All_2015_MAIAC_Biogeograpy.csv"  ,sep='') #
write.csv(Ecotope.data_BS_ad_ha_mean,file =file_ful_path,row.names = F) 


##-------------------------------------------Figure 5 Panel C--------------------------------------
Amazon_biogeography<-Ecotope.data_BS_ad_ha_mean
Amazon_biogeography_data<-Amazon_biogeography[which( is.finite(Amazon_biogeography$EVI_anomaly_fit_standardize) ),]#& 

Amazon_biogeography_data$Type=""

Amazon_biogeography_data[which( Amazon_biogeography_data$WTD <= 10 & Amazon_biogeography_data$SoilFertility <=(-0.35) & Amazon_biogeography_data$TreeHeight<= 32.5),]$Type="67"#"77"#"Shallow_Infertile_Short" 
Amazon_biogeography_data[which( Amazon_biogeography_data$WTD > 10 &
                                  Amazon_biogeography_data$SoilFertility <=(-0.35) & Amazon_biogeography_data$TreeHeight<= 32.5),]$Type="52"#"Deep_Infertile_Short" 
Amazon_biogeography_data[which( Amazon_biogeography_data$WTD <= 10 &
                                  Amazon_biogeography_data$SoilFertility >(-0.35) & Amazon_biogeography_data$TreeHeight<= 32.5),]$Type="48"#"Shallow_Fertile_Short" 
Amazon_biogeography_data[which( Amazon_biogeography_data$WTD <= 10 &
                                  Amazon_biogeography_data$SoilFertility <=(-0.35) & Amazon_biogeography_data$TreeHeight> 32.5),]$Type="75"#"65"#"Shallow_Infertile_Tall" 
Amazon_biogeography_data[which( Amazon_biogeography_data$WTD <= 10 &
                                  Amazon_biogeography_data$SoilFertility >(-0.35) & Amazon_biogeography_data$TreeHeight> 32.5),]$Type="26"#"Shallow_Fertile_Tall" 

Amazon_biogeography_data[which( Amazon_biogeography_data$WTD > 10 &
                                  Amazon_biogeography_data$SoilFertility <=(-0.35) & Amazon_biogeography_data$TreeHeight> 32.5),]$Type="83"#"Deep_Infertile_Tall" 

Amazon_biogeography_data[which( Amazon_biogeography_data$WTD > 10 &
                                  Amazon_biogeography_data$SoilFertility >(-0.35) & Amazon_biogeography_data$TreeHeight<= 32.5),]$Type="11"#"Deep_Fertile_Short"

Amazon_biogeography_data[which( Amazon_biogeography_data$WTD > 10 &
                                  Amazon_biogeography_data$SoilFertility>(-0.35) & Amazon_biogeography_data$TreeHeight> 32.5),]$Type="34"#"Deep_Fertile_Tall" 

threshold_resilience=0.15
Amazon_biogeography_data$Resilience_Type=""
Amazon_biogeography_data[which( Amazon_biogeography_data$EVI_anomaly_fit_standardize <=(-1*threshold_resilience)),]$Resilience_Type="Vulnerable"
Amazon_biogeography_data[which( Amazon_biogeography_data$EVI_anomaly_fit_standardize >=(threshold_resilience)),]$Resilience_Type="Resilient"
Amazon_biogeography_data[which(
  Amazon_biogeography_data$EVI_anomaly_fit_standardize <(threshold_resilience) & Amazon_biogeography_data$EVI_anomaly_fit_standardize >(-1*threshold_resilience)),]$Resilience_Type="uNeutral"




grid.newpage()

# two plots HAND_HAND_Class
rangeMin=0
rangeMax=1
BY_increase=0.25
cbPalette <- c("#006837",  "#88419D",  "#08519c", "#08519c","#D9D9D9","#969696")
cbPalette <- c(  "#08519c",  "#88419D","#006837", "#08519c","#D9D9D9","#969696","#006837",  "#88419D")
cbPalette <- c( "#599F40FF","#F9E8A1FF", "#E46F00FF"  )

Figure5_PanelC_P <- ggplot(data=Amazon_biogeography_data, mapping=aes(x=factor(Type), group =  factor(Resilience_Type), fill = factor(Resilience_Type)    ))+ #,  group = factor(year), fill = factor(year)  
  geom_bar(alpha=I(0.3),colour="#737373", width=0.4,size=0.3,position = "fill") +  #colour="#737373", count  fill=cbPalette[3],
  
  scale_fill_manual(values = cbPalette)+
  xlab("Ecotope Stategies") +
  ylab("Proportion") +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  # scale_x_continuous(limits =c(0,60), breaks=seq(0,60,  by=20))
  theme_few() %+replace% 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1)) +
  theme(axis.text.x = element_text( color="Black", size=24))+ 
  theme(axis.text.y = element_text( color="Black", size=24))+
  theme(axis.title=element_text(size=24))+ 
  theme(plot.title=element_text(size=24))

Figure5_PanelC_P


Figure5_PanelC_P2 <- ggplot(data=Amazon_biogeography_data, mapping=aes(x=factor(Type)   ))+ #,  group = factor(year), fill = factor(year) 
  geom_bar(aes(group =  factor(Resilience_Type), fill = factor(Resilience_Type) ),alpha=I(0.3),colour="#737373", width=0.4,size=0.3,position = "fill") + 
  
  
  geom_bar(aes(y=..count../sum(..count..)),alpha=I(0.3),colour="#737373", width=0.4,size=0.3, position = "dodge") +  #colour="#737373", count  fill=cbPalette[3],

  scale_fill_manual(values = cbPalette)+
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  theme_few() %+replace% 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1)) +
  theme(axis.text.x = element_text( color="Black", size=24))+ 
  theme(axis.text.y = element_text( color="Black", size=24))+
  theme(axis.title=element_text(size=24))+ 
  theme(plot.title=element_text(size=24))

Figure5_PanelC_P2



rangeMin=-1.5
rangeMax=1.5
BY_increase=0.5
Figure5_PanelC_RR <- ggplot(data=Amazon_biogeography_data, mapping=aes(x=factor(Type),y=EVI_anomaly_fit_standardize   ))+ #,  group = factor(year), fill = factor(year) 
  
  geom_bar(stat = "summary", fun = "mean",alpha=I(0.5),fill="white",colour="blue")+
  #geom_errorbar(data=data_tmp,aes(ymin=mean_Resilience, ymax=mean_Resilience,colour="blue"))+
  scale_fill_manual(values = cbPalette)+
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  xlab("Ecotope Stategies") +
  ylab("Relative resilience") +
  #scale_y_continuous(sec.axis = sec_axis(~. , name = "secondary axis")) +
  theme_few() %+replace% 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1)) +
  theme(axis.text.x = element_text( color="Black", size=24))+ 
  theme(axis.text.y = element_text( color="Black", size=24))+
  theme(axis.title=element_text(size=24))+ 
  theme(plot.title=element_text(size=24))

Figure5_PanelC_RR


Figure5_PanelC<-merge_figures_sameunits(Figure5_PanelC_P,Figure5_PanelC_RR)
grid.draw(Figure5_PanelC)


#---------------------------------------function for Figure 5----------------------------


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
  #overlay_g$grobs[[leg_ind]],
 # plt_panel$t,
 # plt_panel$l,
 # plt_panel$b,
 # plt_panel$r, name = "b") #
  return(final_grob)
  
}



Model_Seg_Resilience<- function( data, data.training, model, x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly',  moment = mean) {
  
  ix0 <- match(x0.var, names(data))
  ix1<- match(x1.var, names(data))
  ix2<- match(x2.var, names(data))
  ix3<- match(x3.var, names(data))
  ix4<- match(x4.var, names(data))
  
  ix0.other<- match(x0.other, names(data)) # get index to x-var name 
  ix1.other<- match(x1.other, names(data))
  ix2.other<- match(x2.other, names(data))
  ix3.other<- match(x3.other, names(data))
  
  
  ix0.trn <- match(x0.var, names(data.training))
  ix1.trn<- match(x1.var, names(data.training))
  ix2.trn<- match(x2.var, names(data.training))
  ix3.trn<- match(x3.var, names(data.training))
  ix4.trn<- match(x4.var, names(data.training))
  
  iy <- match(y.var, names(data))
  
  data[,ix0]<-(data[,ix0]-mean(data.training[,ix0.trn]))/sd(data.training[,ix0.trn]) #using the same scale method in training datasets
  print(sd(data.training[,ix0.trn])) 
  data[,ix1]<-(data[,ix1]-mean(data.training[,ix1.trn]))/sd(data.training[,ix1.trn])
  print(sd(data.training[,ix1.trn])) 
  data[,ix2]<-(data[,ix2]-mean(data.training[,ix2.trn]))/sd(data.training[,ix2.trn])
  print(sd(data.training[,ix2.trn])) 
  data[,ix3]<-(data[,ix3]-mean(data.training[,ix3.trn]))/sd(data.training[,ix3.trn])
  print(sd(data.training[,ix3.trn])) 
  x1.level <- 0#moment(data[,ix1])
  x2.level <- 0#moment(data[,ix2])
  x3.level <- 0#moment(data[,ix3])
  x4.level <- 0#moment(data[,ix4])
  
  
  x0.other.level <- 0#moment(data[data[,ix4.other] == Seg,ix0.other]) 
  x1.other.level <- 0#moment(data[data[,ix4.other] == Seg,ix1.other]) 
  x2.other.level <- 0#moment(data[data[,ix4.other] == Seg,ix2.other]) 
  x3.other.level <- 0#moment(data[data[,ix4.other] == Seg,ix3.other]) 
  print(x0.other.level)
  print(x1.other.level)
  print(x2.other.level)
  print(x3.other.level)
  # set the hidden (other) variable 
  # set to a level based on some specified moment (e.g. mean, median, max, min etc)
  row_lnth<-nrow(data)
  x <- "Welcome to Programiz"
  print(row_lnth)
  
  new <- data.frame( data[,ix0],data[,ix1],data[,ix2],+
                       data[,ix3], rep(x4.level, row_lnth),+
                       rep(x0.other.level, row_lnth),rep(x1.other.level,row_lnth),+
                       rep(x2.other.level, row_lnth),rep(x3.other.level, row_lnth) )
  
  names(new) <- c(x0.var,x1.var,x2.var,x3.var,x4.var, x0.other,x1.other,x2.other,x3.other)   
  
  new_fit_array<-add_fitted(new,model)
  return(new_fit_array)
}
