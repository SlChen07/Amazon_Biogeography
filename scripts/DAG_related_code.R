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
                                   Drought_04De_new2.data.data$Long,
                                   Drought_04De_new2.data.data$Lat,
                                   Drought_04De_new2.data.data$SegGeo_number,
                                   Drought_04De_new2.data.data$HAND_CLASS)
names(Drought_04De_new1.data)<-c('PAR_anomaly','VPD_anomaly','Pre_anomaly','MCWD_anomaly','WTD',
                                 'SoilFertility','SoilSand_content','TreeHeight','WoodDensity',
                                 'Drought_Length','MCWD_STD','DrySeasonLength',
                                 'EVI_anomaly','Long','Lat','SegGeo_number','HAND_CLASS')
summary(Drought_04De_new1.data)
Drought_04De_new.data2<-Drought_04De_new1.data  #
Drought_04De_new.data_beforscaling<-Drought_04De_new1.data 
#Scale the data
#Scaled.Drought_04De_new1.data<-apply(Drought_04De_new1.data[,1:12],2, scale)
#Drought_04De_new.data <- as.data.frame( cbind(Scaled.Drought_04De_new1.data,Drought_04De_new1.data[,13:17]))


DAG_ObservationData <- Drought_04De_new1.data

DAG_ObservationData<-DAG_ObservationData[which(is.finite(DAG_ObservationData$EVI_anomaly) 
                                        & is.finite(DAG_ObservationData$PAR_anomaly) 
                                        & is.finite(DAG_ObservationData$VPD_anomaly) 
                                        & (DAG_ObservationData$EVI_anomaly != 0)  
                                        & is.finite(DAG_ObservationData$MCWD_anomaly) 
                                        & is.finite(DAG_ObservationData$SoilFertility)  
                                        & is.finite(DAG_ObservationData$Pre_anomaly) 
                                        & DAG_ObservationData$SegGeo_number >=9  ),] 

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

