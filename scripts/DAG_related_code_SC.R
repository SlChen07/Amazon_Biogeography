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
library(dagitty);library(ggdag);library(readxl);library(lavaan)#library(lavaan) library(flowCore); library(premessa);
# set your own working directory
setwd("F:/csl/HAND/Code_Sum/R/Biogeography_of_Amazon_forests")
############################################################################################################
# DATASETS REQUIRED:
############################################################################################################

## These are the remote sensed aggregated 0.4 degree data; see documentation for definition of variables
## functions are listed at the end of the page
## Read EVI anomaly and climate data
file_path="data/"
file_name=paste('ThreeDroughts_EVI_Ecotope_Climate_04De','.csv',sep='')
file_ful_path=paste(file_path,"/",file_name,sep='')
Drought_04De.data<- read.csv(file=file_ful_path,header=T) 
############################################################################################################
#organize the dataset 
options(digits=7) 
Drought_04De_new2.data.data<-Drought_04De.data[which(is.finite(Drought_04De.data$EVI_anomaly) & is.finite(Drought_04De.data$Tree_Height) &  is.finite(Drought_04De.data$PAR_anomaly) & is.finite(Drought_04De.data$VPD_anomaly) & is.finite(Drought_04De.data$SoilFertility)  & Drought_04De.data$SegGeo_number>=04 & Drought_04De.data$SegGeo_number <=36 & Drought_04De.data$WTD <=5000  &  is.finite(Drought_04De.data$MCWD_STD) & Drought_04De.data$SoilSand_content<1 ),]  # 
Drought_04De_new1.data<-data.frame(Drought_04De_new2.data.data$PAR_anomaly,
                                   Drought_04De_new2.data.data$VPD_anomaly,
                                   Drought_04De_new2.data.data$Pre_anomaly, 
                                   Drought_04De_new2.data.data$MCWD_anomaly,
                                   Drought_04De_new2.data.data$WTD,
                                   Drought_04De_new2.data.data$SoilFertility, 
                                   Drought_04De_new2.data.data$SoilSand_content,
                                   Drought_04De_new2.data.data$Tree_Height,
                                   Drought_04De_new2.data.data$WoodDensity,
                                   Drought_04De_new2.data.data$Drought_Length,
                                   Drought_04De_new2.data.data$MCWD_STD,
                                   Drought_04De_new2.data.data$DrySeasonLength,
                                   Drought_04De_new2.data.data$EVI_anomaly,
                                   Drought_04De_new2.data.data$WTD,
                                   Drought_04De_new2.data.data$Long,
                                   Drought_04De_new2.data.data$Lat,
                                   Drought_04De_new2.data.data$year,
                                   Drought_04De_new2.data.data$SegGeo_number,
                                   Drought_04De_new2.data.data$HAND_CLASS)
names(Drought_04De_new1.data)<-c('PAR_anomaly','VPD_anomaly','Pre_anomaly','MCWD_anomaly','WTD',
                                 'SoilFertility','SoilSand_content','TreeHeight','WoodDensity',
                                 'Drought_Length','MCWD_STD','DrySeasonLength',
                                 'EVI_anomaly','HAND','Long','Lat','year','SegGeo_number','HAND_CLASS')
summary(Drought_04De_new1.data)
Drought_04De_new.data2<-Drought_04De_new1.data  #
Drought_04De_new.data_beforscaling<-Drought_04De_new1.data 


DAG_ObservationData <- Drought_04De_new1.data

DAG_Observation=data.frame(DAG_ObservationData$EVI_anomaly,DAG_ObservationData$PAR_anomaly,
                           DAG_ObservationData$VPD_anomaly,DAG_ObservationData$Pre_anomaly, 
                           DAG_ObservationData$MCWD_anomaly,DAG_ObservationData$WTD,
                           DAG_ObservationData$SoilFertility,DAG_ObservationData$SoilSand_content,
                           DAG_ObservationData$TreeHeight,DAG_ObservationData$WoodDensity,
                           DAG_ObservationData$Drought_Length,DAG_ObservationData$DrySeasonLength)

rm(DAG_ObservationData)
names(DAG_Observation)<-c('DroughtResponses','PARAnomaly','VPDAnomaly','PrecipitationAnomaly','MCWDAnomaly',
                                 'WTD','SoilFertility','SoilTexture','ForestHeight',
                                 'WoodDensity','DroughtLength','DrySeasonLength')


#test the DAG-data consistency
#corr <- lavCor(DAG_Observation )
res <-localTests( dag, DAG_Observation, type="cis.loess",
                  R=100,max.conditioning.variables=4 ) #"cis"
res <- res[order( abs( res$estimate ) ),]
#print( res )
plotLocalTestResults( tail( res, 20), xlim = range(-0.8, 1.0) )



##-------------------------the initial DAG----------------------------------------
dag <- dagitty('dag {
  bb="-6.562,-6.574,5.67,7.958"
  DroughtLength [pos="-4.506,2.672"]
  DroughtResponses [outcome,pos="4.583,-0.391"]
  DrySeasonLength [pos="-3.253,1.251"]
  MCWDAnomaly [pos="1.167,3.539"]
  PARAnomaly [pos="-2.638,6.159"]
  PrecipitationAnomaly [pos="-0.109,6.176"]
  RootDepth [latent,pos="0.139,-3.399"]
  SoilFertility [pos="-4.624,-3.399"]
  SoilTexture [pos="-4.600,-1.585"]
  ForestHeight [pos="-1.811,-4.062"]
  VPDAnomaly [pos="0.789,2.445"]
  WTD [pos="-4.565,0.648"]
  WoodDensity [pos="-0.133,-1.439"]
  DroughtLength -> DroughtResponses
  DrySeasonLength -> DroughtLength
  MCWDAnomaly -> DroughtResponses
  PARAnomaly -> DroughtResponses
  PARAnomaly -> VPDAnomaly
  PrecipitationAnomaly -> MCWDAnomaly
  RootDepth -> DroughtResponses
  SoilFertility -> DroughtResponses
  SoilTexture -> DroughtResponses
  ForestHeight -> RootDepth
  VPDAnomaly -> DroughtResponses
  WTD -> DroughtResponses
  WoodDensity -> DroughtResponses
}')

##-------------------------the final DAG----------------------------------------

#DAG 1
dag <- dagitty('dag{
  bb="-6.562,-6.574,5.67,7.958"
  DroughtLength [pos="-4.506,2.672"]
  DroughtResponses [outcome,pos="4.583,-0.391"]
  DrySeasonLength [pos="-3.253,1.251"]
  MCWDAnomaly [pos="1.167,3.539"]
  OtherLongtermClimateVariability [latent,pos="-1.433,1.251"]
  PARAnomaly [pos="-2.638,6.159"]
  PrecipitationAnomaly [pos="-0.109,6.176"]
  RootDepth [latent,pos="0.139,-3.399"]
  SoilFertility [pos="-4.624,-3.399"]
  SoilTexture [pos="-4.600,-1.585"]
  ForestHeight [pos="-1.811,-4.062"]
  VPDAnomaly [pos="0.789,2.445"]
  WTD [pos="-4.565,0.648"]
  WoodDensity [pos="-0.133,-1.439"]
  DroughtLength -> DroughtResponses
  DroughtLength -> MCWDAnomaly
  DroughtLength -> PARAnomaly
  DroughtLength -> PrecipitationAnomaly
  DrySeasonLength -> DroughtLength
  DrySeasonLength -> DroughtResponses
  DrySeasonLength -> MCWDAnomaly
  DrySeasonLength -> PARAnomaly
  DrySeasonLength -> PrecipitationAnomaly
  DrySeasonLength -> RootDepth
  DrySeasonLength -> SoilFertility
  DrySeasonLength -> ForestHeight
  DrySeasonLength -> VPDAnomaly
  DrySeasonLength -> WoodDensity
  MCWDAnomaly -> DroughtResponses
  OtherLongtermClimateVariability -> PARAnomaly
  OtherLongtermClimateVariability -> WoodDensity
  PARAnomaly -> DroughtResponses
  PARAnomaly -> MCWDAnomaly
  PARAnomaly -> PrecipitationAnomaly
  PARAnomaly -> VPDAnomaly
  PrecipitationAnomaly -> DroughtResponses
  PrecipitationAnomaly -> MCWDAnomaly
  PrecipitationAnomaly -> VPDAnomaly
  RootDepth -> DroughtResponses
  SoilFertility -> DroughtResponses
  SoilFertility -> RootDepth
  SoilFertility -> ForestHeight
  SoilFertility -> WoodDensity
  SoilTexture -> DroughtResponses
  SoilTexture -> RootDepth
  SoilTexture -> SoilFertility
  SoilTexture -> ForestHeight
  SoilTexture -> WoodDensity
  ForestHeight -> DroughtResponses
  ForestHeight -> RootDepth
  ForestHeight -> WoodDensity
  VPDAnomaly -> DroughtResponses
  WTD -> DroughtResponses
  WTD -> RootDepth
  WTD -> ForestHeight
  WTD -> WoodDensity
  WoodDensity -> DroughtResponses
}')


#DAG2
dag <- dagitty('dag {
bb="-6.562,-6.574,5.67,7.958"
DroughtLength [pos="-4.506,2.672"]
DroughtResponses [outcome,pos="4.583,-0.391"]
DrySeasonLength [pos="-3.253,1.251"]
MCWDAnomaly [pos="1.167,3.539"]
OtherLongtermClimateVariability [latent,pos="-1.433,1.251"]
PARAnomaly [pos="-2.638,6.159"]
PrecipitationAnomaly [pos="-0.109,6.176"]
RootDepth [latent,pos="0.139,-3.399"]
SoilFertility [pos="-4.624,-3.399"]
SoilTexture [pos="-4.600,-1.585"]
ForestHeight [pos="-1.811,-4.062"]
VPDAnomaly [pos="0.789,2.445"]
WTD [pos="-4.565,0.648"]
WoodDensity [pos="-0.133,-1.439"]
DroughtLength -> DroughtResponses
DroughtLength -> MCWDAnomaly
DroughtLength -> PARAnomaly
DroughtLength -> PrecipitationAnomaly
DrySeasonLength -> DroughtLength
DrySeasonLength -> DroughtResponses
DrySeasonLength -> MCWDAnomaly
DrySeasonLength -> PARAnomaly
DrySeasonLength -> PrecipitationAnomaly
DrySeasonLength -> RootDepth
DrySeasonLength -> SoilFertility
DrySeasonLength -> ForestHeight
DrySeasonLength -> VPDAnomaly
DrySeasonLength -> WoodDensity
MCWDAnomaly -> DroughtResponses
OtherLongtermClimateVariability -> PARAnomaly
OtherLongtermClimateVariability -> WoodDensity
PARAnomaly -> DroughtResponses
PARAnomaly -> MCWDAnomaly
PARAnomaly -> PrecipitationAnomaly
PARAnomaly -> VPDAnomaly
PrecipitationAnomaly -> DroughtResponses
PrecipitationAnomaly -> MCWDAnomaly
PrecipitationAnomaly -> VPDAnomaly
RootDepth -> DroughtResponses
SoilFertility -> DroughtResponses
SoilFertility -> RootDepth
SoilFertility -> ForestHeight
SoilFertility -> WoodDensity
SoilTexture -> DroughtResponses
SoilTexture -> RootDepth
SoilTexture -> SoilFertility
SoilTexture -> ForestHeight
SoilTexture -> WTD
SoilTexture -> WoodDensity
ForestHeight -> DroughtResponses
ForestHeight -> RootDepth
ForestHeight -> WoodDensity
VPDAnomaly -> DroughtResponses
WTD -> DroughtResponses
WTD -> RootDepth
WTD -> ForestHeight
WTD -> WoodDensity
WoodDensity -> DroughtResponses
}')


##############################################################################
##DAG-derived model
Model_04Degree_DAG<- function(Input_model.data){
  names(Input_model.data)
  
  
  # mod.IA<-mgcv::gam(EVI_anomaly ~    s(WTD) ,method = "ML", data = Input_model.data) # DAG-derived GAM for WTD
  # mod.IA<-mgcv::gam(EVI_anomaly ~ s(PAR_anomaly) +  s(Drought_Length)+s(DrySeasonLength) ,method = "ML", data = Input_model.data) # DAG-derived GAM for PAR
  # mod.IA<-mgcv::gam(EVI_anomaly ~ s(TreeHeight) +s(SoilFertility)+s(SoilSand_content)+s(DrySeasonLength) ,method = "ML", data = Input_model.data) # DAG-derived GAM for Forest height after controlling WTD/HAND
  # mod.IA<-mgcv::gam(EVI_anomaly ~ s(Drought_Length)+s(DrySeasonLength) ,method = "ML", data = Input_model.data) # DAG-derived GAM for drought length
   mod.IA<-mgcv::gam(EVI_anomaly ~ s(SoilFertility,k=5)+s(DrySeasonLength)+s(SoilSand_content) ,method = "REML", data = Input_model.data) # DAG-derived GAM for drought length
  
  #mod.IA<-mgcv::gam(EVI_anomaly ~    s(SoilSand_content,k=10) ,method = "ML", data = Input_model.data) # DAG-derived GAM for soil texture
  return(mod.IA) 
}

#Scale the data
Drought_04De_new1.data2015<-Drought_04De_new1.data[which(Drought_04De_new1.data$year =="2015"),]
Scaled.Drought_04De_new1.data2015<-apply(Drought_04De_new1.data2015[,1:12],2, scale)
Drought_04De_new.data2015 <- as.data.frame( cbind(Scaled.Drought_04De_new1.data2015,Drought_04De_new1.data2015[,13:19]))

Drought_04De_new.data2015_sub<-Drought_04De_new.data2015[which(Drought_04De_new.data2015$HAND>=20  ),] #& Drought_04De_new.data2015$HAND<=40
model_DAG_deep<-Model_04Degree_DAG(Drought_04De_new.data2015_sub)
Drought_04De_new.data2015_sub<-Drought_04De_new.data2015[which(Drought_04De_new.data2015$HAND<10 ),]
model_DAG_shallow<-Model_04Degree_DAG(Drought_04De_new.data2015_sub)
##############################################################################

# draw the plots
SF_arr<-matrix(c((-40:40)/10))
number_group=81
newd0 <- data.frame(SoilFertility=SF_arr,DrySeasonLength=rep(0,times=number_group),SoilSand_content=rep(0,times=number_group)) #,Pre_anomaly=rep(pre_set,times=number_group) SF_set_reverse[1]
prediction_wtd00=predict(model_DAG_shallow,newd0,se.fit = TRUE)
newd0$lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
newd0$fit <- prediction_wtd00$fit
newd0$uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
newd0$type="shallow"
newd0$SoilFertility_reverse=(newd0$SoilFertility*sd(Drought_04De_new1.data2015$SoilFertility))+mean(Drought_04De_new1.data2015$SoilFertility)

newd1 <- data.frame(SoilFertility=SF_arr,DrySeasonLength=rep(0,times=number_group),SoilSand_content=rep(0,times=number_group)) #,Pre_anomaly=rep(pre_set,times=number_group) SF_set_reverse[1]
prediction_wtd00=predict(model_DAG_deep,newd1,se.fit = TRUE)
newd1$lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
newd1$fit <- prediction_wtd00$fit
newd1$uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
newd1$type="deep"
newd1$SoilFertility_reverse=(newd1$SoilFertility*sd(Drought_04De_new1.data2015$SoilFertility))+mean(Drought_04De_new1.data2015$SoilFertility)


newd0<-newd0[which(newd0$SoilFertility_reverse >=-1.0 &  newd0$SoilFertility_reverse<=0.2 )  ,]
newd1<-newd1[which(newd1$SoilFertility_reverse >=-1.0 &  newd1$SoilFertility_reverse<=0.2 )  ,]


#draw 2D figure
cbPalette<- c( "#416E9BFF","#C4263DFF")
#cbPalette<- c("#C4263DFF")
rangeMin=-1.5
range_mm=2
rangeMax=rangeMin+range_mm
BY_increase=0.5
DAG_fig_shallow<-Draw_Fig_DAG(newd0,x0.var = 'SoilFertility_reverse',y.var='fit',x0.range=c(-1.0, 0.1),x0.break=seq(-1.05,0.15, by=0.3),
             y.range=c(rangeMin,rangeMax), y.break=seq(rangeMin,rangeMax,  by=BY_increase),cbPalette[1])
DAG_fig_shallow
DAG_fig_deep<-Draw_Fig_DAG(newd1,x0.var = 'SoilFertility_reverse',y.var='fit',x0.range=c(-1.0, 0.1),x0.break=seq(-1.05,0.15, by=0.3),
                      y.range=c(rangeMin,rangeMax), y.break=seq(rangeMin,rangeMax,  by=BY_increase),cbPalette[2])
DAG_fig_deep
DAG_fig<-merge_figures_sameunits(DAG_fig_shallow,DAG_fig_deep)
grid.draw(DAG_fig)


Draw_Fig_DAG <- function(data,x0.var = 'SoilFertility_reverse',y.var='fit',x0.range=c(-1.0, 0.2),x0.break=seq(-1.0,0.2, by=0.4),
                        y.range=c(rangeMin,rangeMax), y.break=seq(rangeMin,rangeMax,  by=BY_increase),cbPalette)
{ 
  p_fig<-ggplot(data, aes(x = SoilFertility_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =y.range, breaks=y.break) + 
  #  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 0,  width=.6) +
  geom_smooth( aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity", linetype = "dashed") +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
    xlab("Soil fertility (cmol(+)/kg)") +
    ylab("EVI anomaly") +
  scale_x_continuous(limits = x0.range, breaks=x0.break)+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
  return(p_fig)
}

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

