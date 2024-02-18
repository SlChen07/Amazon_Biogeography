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
file_path="data/"
file_name=paste('ThreeDroughts_EVI_Ecotope_Climate_04De','.csv',sep='')
file_ful_path=paste(file_path,"/",file_name,sep='')
Drought_04De.data<- read.csv(file=file_ful_path,header=T) 
############################################################################################################
# Functions to run:
#Function for Figure 4 Panel A & C
Model_2D_Partial_Prediction<- function( data,model,x0.var = 'WTD', x1.var = 'TreeHeight', x3.smooth='smooth',y.smooth='est', x0.smooth='s(WTD)',x1.smooth='s(TreeHeight)',x2.smooth='ti(WTD,TreeHeight)', moment = mean) {
  
  ix0 <- match(x0.var, names(data))
  ix1<- match(x1.var, names(data))
  
  
  sm_ml <- as.data.frame(smooth_estimates(model, dist = 0.1))
  print(1)
  sm_ml.tmp<-sm_ml
  ix0_sm <- match(x0.var, names(sm_ml))
  
  ix1_sm<- match(x1.var, names(sm_ml))
  ix0.smooth_sm<- match(x0.smooth, names(sm_ml))
  ix1.smooth_sm<- match(x1.smooth, names(sm_ml))
  ix2.smooth_sm<- match(x2.smooth, names(sm_ml))
  ix3.smooth_sm<- match(x3.smooth, names(sm_ml))
  y.smooth_sm<- match(y.smooth, names(sm_ml))
  
  sm_ml.tmp[,ix0_sm]<- sm_ml[,ix0_sm]*sd(data[,ix0])+mean(data[,ix0])
  sm_ml.tmp[,ix1_sm]<- sm_ml[,ix1_sm]*sd(data[,ix1])+mean(data[,ix1])
  print(sd(data[,ix0]))
  print(mean(data[,ix1]))
  sm_ml.copy<-sm_ml.tmp
  x1_arr<-sm_ml.copy[sm_ml.copy[,ix3.smooth_sm] == x1.smooth,ix1_sm]
  #print(x1_arr)
  for (TH_val in x1_arr){
    sm_ml.copy[(sm_ml.copy[,ix3.smooth_sm] == x2.smooth & sm_ml.copy[,ix1_sm] == TH_val), y.smooth_sm]<- sm_ml.copy[(sm_ml.copy[,ix3.smooth_sm] == x2.smooth & sm_ml.copy[,ix1_sm] == TH_val), y.smooth_sm]+  sm_ml.copy[(sm_ml.copy[,ix3.smooth_sm] == x1.smooth & sm_ml.copy[,ix1_sm] == TH_val), y.smooth_sm] 
  }  
  
  x0_arr<-sm_ml.copy[sm_ml.copy[,ix3.smooth_sm] == x0.smooth,ix0_sm]
  for (WTD_val in x0_arr){
    sm_ml.copy[sm_ml.copy[,ix3.smooth_sm] == x2.smooth & sm_ml.copy[,ix0_sm] == WTD_val, y.smooth_sm]<- sm_ml.copy[sm_ml.copy[,ix3.smooth_sm] == x2.smooth & sm_ml.copy[,ix0_sm] == WTD_val, y.smooth_sm]+  sm_ml.copy[sm_ml.copy[,ix3.smooth_sm] == x0.smooth & sm_ml.copy[,ix0_sm] == WTD_val, y.smooth_sm] 
  }
  
  return(sm_ml.copy)
}


#draw 2D figure
Draw_Fig_2D <- function(data,x0.var = 'WTD', x1.var = 'TreeHeight',y.var='est',x0.range=c(0,40),x0.break=seq(0,40,  by=10),x1.range=c(20,40),x1.break=seq(21,39,  by=1.5),threshold_value=0.4,c_break=c(-0.2,-0.1,0,0.1,0.2))
{
  p_draw<- ggplot(data, aes_string(x =x0.var, y =x1.var)) +
    geom_raster(aes_string(fill = (y.var))) +
    geom_contour(aes_string(z = (y.var)), colour = "black",linetype="dashed") +
    scale_fill_gradientn(limits = c((-1)*threshold_value,threshold_value),
                         
                         colors = brewer.pal(10,"RdYlGn")[3:9],
                         breaks=(as.numeric(c_break)), labels=format(c_break))+ 
    scale_x_continuous(limits =as.numeric(x0.range))+
    scale_y_continuous(limits =as.numeric(x1.range),breaks=as.numeric(x1.break)) + theme_few() %+replace%
    theme(panel.background = element_rect(fill = NA,colour = "black", 
                                          size =1.0))+
    theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
    theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
    theme(axis.title=element_text(size=18))+ #,face="bold"
    theme(plot.title=element_text(size=18,face="bold"))
  
  return(p_draw)  
}


Data_threshold<-function(data,y.var='est',threshold_value=0.4)
{
  iy <- match(y.var, names(data))
  
  data[data[, iy]>=threshold_value & !is.na(data[, iy]),iy]<-threshold_value
  data[data[, iy]<=(-1)*threshold_value & !is.na(data[, iy]),iy]<-(-1)*threshold_value
  data[data[, iy]>=threshold_value & !is.na(data[, iy]),iy]<-threshold_value
  data[data[, iy]<=(-1)*threshold_value & !is.na(data[, iy]),iy]<-(-1)*threshold_value
  
  return(data)
}


#functions for figure 4 Panel C & D
Model_Seg_Prediction_Region_SoilFertility<- function( data, model,WTD_Array, x0.ori = 'WTD_ori',  x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number', SF.value=0, moment = mean) {
  
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
  
  
  x1.level <- SF.value #moment(data[data[,ix4.other]== Seg,ix1])
  x2.level <- moment(data[,ix2])
  x3.level <- moment(data[,ix3])
  x4.level <- moment(data[,ix4])
  print(x1.level)
  print(x2.level)
  print(x3.level)
  print(x4.level)
  
  x0.other.level <- moment(data[,ix0.other]) 
  x1.other.level <- moment(data[,ix1.other]) 
  x2.other.level <- moment(data[,ix2.other]) 
  x3.other.level <- moment(data[,ix3.other]) 
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
  #print(prediction_wtd00)
  
  prediction.lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
  #print(prediction.lci)
  prediction.fit <- prediction_wtd00$fit
  prediction.uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
  
  WTD_reverse=(new[,1]*sd(data[,ix0.ori]))+mean(data[,ix0.ori])
  
  new_fit_array<-as.data.frame(cbind(new,prediction.lci,prediction.fit,prediction.uci,WTD_reverse))
  names(new_fit_array) <- c(x0.var,x1.var,x2.var,x3.var,x4.var, x0.other,x1.other,x2.other,x3.other,'lci','fit','uci','WTD_reverse' ) 
  
  return(new_fit_array)
}



Model_Seg_Prediction_Region_TreeHeight<- function( data, model,WTD_Array, x0.ori = 'WTD_ori',  x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number', TH.value=0, moment = mean) {
  
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
  
  x1.level <- moment(data[,ix1])
  x2.level <- moment(data[,ix2])
  x3.level <- TH.value#moment(data[,ix3])
  x4.level <- moment(data[,ix4])
  print(x1.level)
  print(x2.level)
  print(x3.level)
  print(x4.level)
  
  x0.other.level <- moment(data[,ix0.other]) 
  x1.other.level <- moment(data[,ix1.other]) 
  x2.other.level <- moment(data[,ix2.other]) 
  x3.other.level <- moment(data[,ix3.other]) 
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
  #print(prediction_wtd00)
  
  prediction.lci <- prediction_wtd00$fit - 1.96 * prediction_wtd00$se.fit
  #print(prediction.lci)
  prediction.fit <- prediction_wtd00$fit
  prediction.uci <- prediction_wtd00$fit + 1.96 * prediction_wtd00$se.fit
  
  WTD_reverse=(new[,1]*sd(data[,ix0.ori]))+mean(data[,ix0.ori])
  
  new_fit_array<-as.data.frame(cbind(new,prediction.lci,prediction.fit,prediction.uci,WTD_reverse))
  names(new_fit_array) <- c(x0.var,x1.var,x2.var,x3.var,x4.var, x0.other,x1.other,x2.other,x3.other,'lci','fit','uci','WTD_reverse' ) 
  
  return(new_fit_array)
}


#functions for panel E
summary_group_full<-function(data,x0.var='EVI_anomaly',x1.var='EVI_anomaly_corrected', flag=1) # if using EVI_anomaly_corrected, flag=1 
{
  Summary_Geo=summary_group_simple(data,flag)
  Summary_Geo<-Summary_Geo[which(Summary_Geo$N >=4 & Summary_Geo$HAND_CLASS <=40  & Summary_Geo$HAND_CLASS >=2),]
  Summary_Geo<-data.frame(Summary_Geo$HAND_CLASS,Summary_Geo$EVI_anomaly_corrected,Summary_Geo$ci,Summary_Geo$se)
  names(Summary_Geo) <- c('HAND_CLASS',x0.var,'ci','se')  # Brando_summary_Geo$ci <=1.5 ),]
  return(Summary_Geo)
}

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


#establish a predictive selected GAM
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
############################################################################################################
############################################################################################################
####-------------------------------------------main code for Figure 4-------------------------------------------
##----------------------------------Panel A 2D-interactions among ecotope factors-------------------------------
options(digits=7) 
Drought_04De_new2.data.data<-Drought_04De.data[which(is.finite(Drought_04De.data$EVI_anomaly) & is.finite(Drought_04De.data$Tree_Height) &  is.finite(Drought_04De.data$PAR_anomaly) & is.finite(Drought_04De.data$VPD_anomaly) & is.finite(Drought_04De.data$SoilFertility)  & Drought_04De.data$SegGeo_number>=04 & Drought_04De.data$SegGeo_number <=36 & Drought_04De.data$WTD <=60  &  is.finite(Drought_04De.data$MCWD_STD) & Drought_04De.data$year=="2015" & Drought_04De.data$SoilSand_content<1 ),]  # 
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
                                      Drought_04De_new2.data.data$HAND_CLASS)#SoilSand_content

names(Drought_04De_new1.data)<-c('PAR_anomaly','VPD_anomaly','Pre_anomaly','MCWD_anomaly','WTD',
                                    'SoilFertility','SoilSand_content','TreeHeight','WoodDensity',
                                    'Drought_Length','MCWD_STD','DrySeasonLength',
                                    'EVI_anomaly','Long','Lat','SegGeo_number','HAND_CLASS')

summary(Drought_04De_new1.data)
#cor(Drought_04De_new1.data)

Drought_04De_new1.data$WoodDensity<-Drought_04De_new1.data$WoodDensity*100.0
Drought_04De_new.data2<-Drought_04De_new1.data  #before scaling
Drought_04De_new.data_beforscaling<-Drought_04De_new1.data 
#Scale the data
Scaled.Drought_04De_new1.data<-apply(Drought_04De_new1.data[,1:12],2, scale)
Drought_04De_new.data <- as.data.frame( cbind(Scaled.Drought_04De_new1.data,Drought_04De_new1.data[,13:17]))
#threshold_scale=10
#Drought_04De_new.data <- Drought_04De_new.data[which(abs(Drought_04De_new.data$EVI_anomaly) <= 10 &    abs(Drought_04De_new.data$PAR_anomaly)<= threshold_scale &  abs(Drought_04De_new.data$VPD_anomaly)<= threshold_scale & abs(Drought_04De_new.data$Pre_anomaly)<= threshold_scale &
                                                              # abs(Drought_04De_new.data$MCWD_anomaly)<= threshold_scale &  abs(Drought_04De_new.data$SoilFertility)<= threshold_scale & abs(Drought_04De_new.data$Drought_Length)<= threshold_scale &  abs(Drought_04De_new.data$DrySeasonLength)<= threshold_scale & abs(Drought_04De_new.data$WTD)<= threshold_scale & abs(Drought_04De_new.data$TreeHeight)<= threshold_scale) ,]

#Run and fit the selected GAM 
mod.IA_Guiana=ancova_establish_slope_Guiana_04Degree_SoilSand(Drought_04De_new.data)
BIC(mod.IA_Guiana)
summary(mod.IA_Guiana)
gam.check(mod.IA_Guiana) 
#concurvity(mod.IA_Guiana, full = TRUE)
FulxTower_adjust_new.data$Residual_R=resid(mod.IA_Guiana)#
FulxTower_adjust_new.data$Prediction=fitted(mod.IA_Guiana)

##----------------------------------Panel A 2D of soil Fertility & WTD continued--------------------------------
#get the fitted smooth
sm <- smooth_estimates(mod.IA_Guiana, dist = 0.1)
sm_copy<-sm
Drought_04De_new.data_BStmp<-Drought_04De_new.data_beforscaling
sm_copy2<-Model_2D_Partial_Prediction(Drought_04De_new.data_BStmp,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility', x3.smooth='smooth',y.smooth='est', x0.smooth='s(WTD)',x1.smooth='s(SoilFertility)',x2.smooth='ti(WTD,SoilFertility)', moment = median)
threshold_value=0.4
sm_copy2[which(sm_copy2$est >=threshold_value),]$est=threshold_value
sm_copy2[which(sm_copy2$est <=(-1)*threshold_value),]$est=(-1)*threshold_value
sm_copy[which(sm_copy$est >=threshold_value),]$est=threshold_value
sm_copy[which(sm_copy$est <=(-1)*threshold_value),]$est=(-1)*threshold_value
c_break<-c((-1)*threshold_value,-0.2,0,0.2,threshold_value)
Figure4_PanelA<-ggplot(sm_copy2, aes(x = WTD, y = SoilFertility)) +
  geom_raster(aes(fill = est)) +
  geom_contour(aes(z = est), colour = "black",linetype="dashed") + #, colour = "gray"
  scale_fill_gradientn(limits = c((-1)*threshold_value,threshold_value),
                       colors = brewer.pal(10,"RdYlGn")[3:9],
                       breaks=(c_break), labels=format(c_break))+
  xlab("HAND (meter)") +
  ylab("Soil fertility (cmol(+)/kg)") +
  scale_x_continuous(limits =c(0,40)) +
  scale_y_continuous(limits =c(-1.1,0.3),breaks=seq(-1.0,0.2,  by=0.1)) + theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
Figure4_PanelA
##--------------------------------------------------------------------------------------------------------------


##----------------------------------Panel B 2D of forest height & WTD continued---------------------------------
sm <- smooth_estimates(mod.IA_Guiana, dist = 0.1)
sm_copy<-sm
Drought_04De_new.data_BStmp<-Drought_04De_new.data_beforscaling
sm_copy2<-Model_2D_Partial_Prediction(Drought_04De_new.data_BStmp,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'TreeHeight', x3.smooth='smooth',y.smooth='est', x0.smooth='s(WTD)',x1.smooth='s(TreeHeight)',x2.smooth='ti(WTD,TreeHeight)', moment = median)
threshold_value=0.4
sm_copy2[which(sm_copy2$est >=threshold_value),]$est=threshold_value
sm_copy2[which(sm_copy2$est <=(-1)*threshold_value),]$est=(-1)*threshold_value
sm_copy[which(sm_copy$est >=threshold_value),]$est=threshold_value
sm_copy[which(sm_copy$est <=(-1)*threshold_value),]$est=(-1)*threshold_value
c_break<-c((-1)*threshold_value,-0.2,0,0.2,threshold_value)
paletteer_c("grDevices::RdYlGn", 30) #"grDevices::RdYlGn"
Figure4_PanelB<-ggplot(sm_copy2, aes(x = WTD, y = TreeHeight)) +
  geom_raster(aes(fill = est)) +
  geom_contour(aes(z = est), colour = "black",linetype="dashed") +
  scale_fill_gradientn(limits = c((-1)*threshold_value,threshold_value),
                       colors = brewer.pal(10,"RdYlGn")[3:9],
                       breaks=(c_break), labels=format(c_break))+ 
  xlab("HAND (meter)") +
  ylab("Tree hight (m)") +
  scale_x_continuous(limits =c(-0,40))+
  scale_y_continuous(limits =c(20,40),breaks=seq(21,39,  by=1.5)) + theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))

Figure4_PanelB
##--------------------------------------------------------------------------------------------------------------

##------------------------------------Model selection-----------------------------------------------------------
mod.AmzBiogeograhy<-mgcv::gam(EVI_anomaly ~   s(WTD,k=5) +s(SoilFertility,k=5)+ti(WTD,SoilFertility,k=c(3,3))
                  +s(TreeHeight,k=5)+ti(WTD,TreeHeight,k=c(3,3)) +s(DrySeasonLength)+ti(WTD,DrySeasonLength,k=c(3,3))
                  +ti(WTD,PAR_anomaly,bs='tp') +ti(WTD,VPD_anomaly,bs='tp') +ti(WTD,MCWD_anomaly,bs='tp')
                  +ti(WTD,Pre_anomaly,bs='tp')+s(PAR_anomaly)+s(Pre_anomaly)+s(VPD_anomaly)+s(MCWD_anomaly)
                  +ti(PAR_anomaly,MCWD_anomaly,bs='tp')+ti(MCWD_anomaly,VPD_anomaly,bs='tp')
                  +ti(VPD_anomaly,PAR_anomaly,bs='tp')+ti(PAR_anomaly,Pre_anomaly,bs='tp')
                  +ti(VPD_anomaly,Pre_anomaly,bs='tp')+ti(MCWD_anomaly,Pre_anomaly,bs='tp') ,method = "ML",data = FulxTower_adjust_new.data) #,method = "REML"

summary(mod.AmzBiogeograhy) 
#getAllTerms(mod.AmzBiogeograhy)
#AICc(mod.AmzBiogeograhy)
AIC(mod.AmzBiogeograhy)
#using dredge to select the model based on AICc/AIC
options(na.action = "na.fail")
dd<-dredge(mod.AmzBiogeograhy, trace =2)
#dd<-dredge(mod.AmzBiogeograhy, subset=  's(SoilFertility, k = 5)' & 's(TreeHeight, k = 5)' & 's(WTD, k = 5)' & 'ti(WTD, SoilFertility, k = c(3, 3))' & 'ti(WTD, TreeHeight, k = c(3, 3))' ,trace =2)
subset(dd, delta < 10)
#summary(dd[1,])
#importance(dd)
summary(mod.AmzBiogeograhy(dd, subset = delta < 10))

##--------------------------------------------------------------------------------------------------------------



##------------------------------------------ Panel C  soil Fertility--------------------------------------------
#set up the datasets
Drought_04De_new.data_BStmp<-Drought_04De_new.data_beforscaling
Ecotope_Factors <- select(Drought_04De_new.data_BStmp, SoilFertility, DrySeasonLength,WTD,TreeHeight,SoilSand_content)
names(Ecotope_Factors)<-c('SoilFertility_ori','DrySeasonLength_ori','WTD_ori','TreeHeight_ori','SoilSand_content_ori')
Drought_04De_new.data_BS<-as.data.frame(cbind(Drought_04De_new.data_BStmp,Ecotope_Factors))
Scaled.FulxTower_new.data_draw<-apply(Drought_04De_new.data_BS[,1:12],2, scale)
Drought_04De_new.data_BS <- as.data.frame( cbind(Scaled.FulxTower_new.data_draw,Drought_04De_new.data_BS[,15:22]))
#threshold_scale=10
#Drought_04De_new.data_BS <- Drought_04De_new.data_BS[which(abs(Drought_04De_new.data_BS$EVI_anomaly) <= 10 &    abs(Drought_04De_new.data_BS$PAR_anomaly)<= threshold_scale &  abs(Drought_04De_new.data_BS$VPD_anomaly)<= threshold_scale & abs(Drought_04De_new.data_BS$Pre_anomaly)<= threshold_scale &
                                                                         #abs(Drought_04De_new.data_BS$MCWD_anomaly)<= threshold_scale &  abs(Drought_04De_new.data_BS$SoilFertility)<= threshold_scale & abs(Drought_04De_new.data_BS$Drought_Length)<= threshold_scale &  abs(Drought_04De_new.data_BS$DrySeasonLength)<= threshold_scale & abs(Drought_04De_new.data_BS$WTD)<= threshold_scale & abs(Drought_04De_new.data_BS$TreeHeight)<= threshold_scale) ,]

#set parameters for soil fertility 
WTD_arr<-matrix(c((-20:60)/10)) 
number_group=81 
SF_set_reverse=((-10:4)*0.1-mean(Drought_04De_new.data_beforscaling$SoilFertility))/sd(Drought_04De_new.data_beforscaling$SoilFertility)
sf_inital=-1.0
byincrease=0.1
#calculate prediction for different levels of soil Fertility
for(i in seq(from=1, 13, by=1))
{
  Model_fit_array<-Model_Seg_Prediction_Region_SoilFertility(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',SF.value=SF_set_reverse[i], moment = mean)
  
  newd0<-Model_fit_array
  
  if ((i-1)<10) {st_tag=paste(as.character(0),as.character(i-1))}else {
    st_tag=c(as.character(i-1))
  }
  st_tag<-gsub(" ", "", st_tag)
  print(st_tag)
  type_str<-paste('SF',st_tag,as.character(sf_inital+byincrease*(i-1.0)),sep = "_")
  newd0$type=type_str
  print(type_str)
  
  if (i == 1 ) {
    newd=rbind(newd0)
  }
  newd<-rbind(newd,newd0)
}

tgca_combine<-newd
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0 &  tgca_combine$WTD_reverse<=40 )  ,]
rangeMin=-1.5
range_mm=2.0
rangeMax=rangeMin+range_mm
BY_increase=0.5
cbPalette<-c("#26456EFF", "#1F5691FF", "#1C6AA9FF", "#2D7DB4FF", "#4A93C1FF", "#7BB2D2FF", "#CACACAFF", "#A9C399FF", "#7CB070FF", "#4E9B51FF", "#29863DFF", "#1B7333FF", "#09622AFF") 
Figure4_PanelC1<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  geom_line(size=1) +
  geom_hline(yintercept = 0,  width=.6) +
  # geom_smooth(aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +  
  xlab("HAND (meter)") +
  ylab("EVI anomaly") +
  scale_x_continuous(breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
geom_point()
Figure4_PanelC1
##--------------------------------------------------------------------------------------------------------------

##----------------------------------Panel C soil fertility shade continued-------------------------------------- 
WTD_arr<-matrix(c((-20:60)/10)) 
number_group=81 
SF_set_reverse=(c(-0.9,-0.7,-1.05,0.125,0.3,-0.1)-mean(Drought_04De_new.data_beforscaling$SoilFertility))/sd(Drought_04De_new.data_beforscaling$SoilFertility)
SF_set<-c(-0.9,-0.7,-1.05,0.125,0.3,-0.1)
sf_inital=-1.0
byincrease=0.1
#calculate prediction for different levles of soil Fertility
for(i in seq(from=1, 6, by=1))
{
  Model_fit_array<-Model_Seg_Prediction_Region_SoilFertility(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',SF.value=SF_set_reverse[i], moment = mean)
  
  newd_tmp<-Model_fit_array
  
  if ((i-1)<10) {st_tag=paste(as.character(0),as.character(i-1))}else {
    st_tag=c(as.character(i-1))
  }
  st_tag<-gsub(" ", "", st_tag)
  print(st_tag)
  type_str<-paste('SF',st_tag,as.character(SF_set[i]),sep = "_")
  newd_tmp$type=type_str
  print(type_str)
  
  if (i == 1 ) {
    newd0<-newd_tmp
    newd=(newd_tmp)
  } else{ newd<-rbind(newd,newd_tmp)}
  if (i == 4 ) {newd3<-newd_tmp}
  
}
#find the confidential intervals
type_inx <- match('type', names(newd))
fit_inx <- match('fit', names(newd))
df1<-data.frame(newd[newd[,type_inx]=="SF_00_-0.9",fit_inx],newd[newd[,type_inx]=="SF_01_-0.7",fit_inx],newd[newd[,type_inx]=="SF_02_-1.05",fit_inx])
newd0$uci<-apply(X = df1[,], MARGIN = 1, FUN = max, na.rm = TRUE)
newd0$fit<-newd0$fit
newd0$lci<-apply(X = df1[,], MARGIN = 1, FUN = min, na.rm = TRUE)
df1<-data.frame(newd[newd[,type_inx]=="SF_03_0.125",fit_inx],newd[newd[,type_inx]=="SF_04_0.3",fit_inx],newd[newd[,type_inx]=="SF_05_-0.1",fit_inx])
newd3$uci<-apply(X = df1[,], MARGIN = 1, FUN = max, na.rm = TRUE)
newd3$fit<-newd3$fit
newd3$lci<-apply(X = df1[,], MARGIN = 1, FUN = min, na.rm = TRUE)

rangeMin=-1.5
range_mm=2.0
rangeMax=rangeMin+range_mm
BY_increase=0.5
tgca_combine<-rbind(newd0,newd3) #
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0 &  tgca_combine$WTD_reverse<=40)  ,]
cbPalette<- c("#1F5691FF", "#1B7333FF")
Figure4_PanelC_shde<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  #  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 0,  width=.6) +
  #geom_smooth( aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") +
  geom_smooth( aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity", linetype = "blank") + #dashed
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  xlab("HAND (meter)") +
  ylab("EVI anomaly") +
  scale_x_continuous(breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
Figure4_PanelC_shde
#merge figures
Figure4_PanelC<-merge_figures_sameunits(Figure4_PanelC1,Figure4_PanelC_shde)
grid.draw(Figure4_PanelC)
##--------------------------------------------------------------------------------------------------------------
##-------------------------------------Figure 4 Panel D  Forest Height--------------------------------------------
#set parameters for tree height 
WTD_arr<-matrix(c((-20:60)/10)) 
number_group=81 
TH_set_reverse=((14:26)*1.5-mean(Drought_04De_new.data_beforscaling$TreeHeight))/sd(Drought_04De_new.data_beforscaling$TreeHeight)
th_inital=21
byincrease=1.5
#calculate prediction for different levels of tree height
for(i in seq(from=1, 13, by=1))
{
  Model_fit_array<-Model_Seg_Prediction_Region_TreeHeight(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',TH.value=TH_set_reverse[i], moment = mean)
  
  newd0<-Model_fit_array
  
  if ((i-1)<10) {st_tag=paste(as.character(0),as.character(i-1))}else {
    st_tag=c(as.character(i-1))
  }
  st_tag<-gsub(" ", "", st_tag)
  print(st_tag)
  type_str<-paste('TH',st_tag,as.character(th_inital+byincrease*(i-1.0)),sep = "_")
  newd0$type=type_str
  print(type_str)
  
  if (i == 1 ) {
    newd=rbind(newd0)
  }
  
  newd<-rbind(newd,newd0)
  
}
#Draw the plot
tgca_combine<-newd
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0 &  tgca_combine$WTD_reverse<=40 )  ,]
rangeMin=-1.5
range_mm=2.0
rangeMax=rangeMin+range_mm
BY_increase=0.5
cbPalette<- c("#2E5A87FF", "#416E9BFF", "#5383AFFF", "#6998C1FF", "#7EAED3FF", "#B2C1D2FF", "#DFD4D1FF", "#F1AB9CFF", "#F87F69FF", "#ED6055FF", "#E03B42FF", "#C4263DFF", "#A90C38FF") 
Figure4_PanelD1<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  geom_line(size=1) +
  geom_hline(yintercept = 0,  width=.6) +
  # geom_smooth(aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +  
  xlab("HAND (meter)") +
  ylab("EVI anomaly") +
  scale_x_continuous(breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
Figure4_PanelD1
##--------------------------------------------------------------------------------------------------------------

##----------------------------------Panel C forest Height shade continued--------------------------------------- 
##set parameters for forest Height predictions
WTD_arr<-matrix(c((-20:60)/10)) 
number_group=81 
TH_set_reverse=(c(26.5,21,31.5,37.5,39.5,33)-mean(Drought_04De_new.data_beforscaling$TreeHeight))/sd(Drought_04De_new.data_beforscaling$TreeHeight)
TH_set<-c(26.5,21,31.5,37.5,39.5,33)
th_inital=21
byincrease=1.5
#calculate prediction for different levels of forest height
for(i in seq(from=1, 6, by=1))
{
  Model_fit_array<-Model_Seg_Prediction_Region_TreeHeight(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',TH.value=TH_set_reverse[i], moment = mean)
  
  newd_tmp<-Model_fit_array
  
  if ((i-1)<10) {st_tag=paste(as.character(0),as.character(i-1))}else {
    st_tag=c(as.character(i-1))
  }
  st_tag<-gsub(" ", "", st_tag)
  print(st_tag)
  type_str<-paste('TH',st_tag,as.character(TH_set[i]),sep = "_")
  newd_tmp$type=type_str
  print(type_str)
  
  if (i == 1 ) {
    newd0<-newd_tmp
    newd=(newd_tmp)
  } else{ newd<-rbind(newd,newd_tmp)}
  if (i == 4 ) {newd3<-newd_tmp}
  
}
#find the condidential intervals
type_inx <- match('type', names(newd))
fit_inx <- match('fit', names(newd))
df1<-data.frame(newd[newd[,type_inx]=="TH_00_26.5",fit_inx],newd[newd[,type_inx]=="TH_01_21",fit_inx],newd[newd[,type_inx]=="TH_02_31.5",fit_inx])
newd0$uci<-apply(X = df1[,], MARGIN = 1, FUN = max, na.rm = TRUE)
newd0$fit<-newd0$fit
newd0$lci<-apply(X = df1[,], MARGIN = 1, FUN = min, na.rm = TRUE)
df1<-data.frame(newd[newd[,type_inx]=="TH_03_37.5",fit_inx],newd[newd[,type_inx]=="TH_04_39.5",fit_inx],newd[newd[,type_inx]=="TH_05_33",fit_inx])
newd3$uci<-apply(X = df1[,], MARGIN = 1, FUN = max, na.rm = TRUE)
newd3$fit<-newd3$fit
newd3$lci<-apply(X = df1[,], MARGIN = 1, FUN = min, na.rm = TRUE)
#draw the plot
rangeMin=-1.5
range_mm=2.0
rangeMax=rangeMin+range_mm
BY_increase=0.5
tgca_combine<-rbind(newd0,newd3) 
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0 &  tgca_combine$WTD_reverse<=40)  ,]
cbPalette<- c("#7EAED3FF","#E03B42FF")
Figure4_PanelD_shde<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  #  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 0,  width=.6) +
  #geom_smooth( aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity") +
  geom_smooth( aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity", linetype = "blank") + #dashed
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  xlab("HAND (meter)") +
  ylab("EVI anomaly") +
  scale_x_continuous(breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))
Figure4_PanelD_shde
#merge figures
Figure4_PanelD<-merge_figures_sameunits(Figure4_PanelD1,Figure4_PanelD_shde)
grid.draw(Figure4_PanelD)
##--------------------------------------------------------------------------------------------------------------

##------------------------------------------Figure 4 Panel E ---------------------------------------------------
#---------------------------------------Figure 4 Panel E observation--------------------------------------------
year=2015
Drought_04De_new.data_BStmp<-Drought_04De_new.data_beforscaling
Ecotope_Factors <- select(Drought_04De_new.data_BStmp, SoilFertility, DrySeasonLength,WTD,TreeHeight,SoilSand_content)
names(Ecotope_Factors)<-c('SoilFertility_ori','DrySeasonLength_ori','WTD_ori','TreeHeight_ori','SoilSand_content_ori')
Drought_04De_new.data_BS<-as.data.frame(cbind(Drought_04De_new.data_BStmp,Ecotope_Factors))
#calculate the same process
Scaled.Drought_04De_new.data_BS<-apply(Drought_04De_new.data_BS[,1:12],2, scale)
Drought_04De_new.data_BS <- as.data.frame( cbind(Scaled.Drought_04De_new.data_BS,Drought_04De_new.data_BS[,13:22]))
#threshold_scale=10
#Drought_04De_new.data_BS <- Drought_04De_new.data_BS[which(abs(Drought_04De_new.data_BS$EVI_anomaly) <= 10 &    abs(Drought_04De_new.data_BS$PAR_anomaly)<= threshold_scale &  abs(Drought_04De_new.data_BS$VPD_anomaly)<= threshold_scale & abs(Drought_04De_new.data_BS$Pre_anomaly)<= threshold_scale &
                                                                        # abs(Drought_04De_new.data_BS$MCWD_anomaly)<= threshold_scale &  abs(Drought_04De_new.data_BS$SoilFertility)<= threshold_scale & abs(Drought_04De_new.data_BS$Drought_Length)<= threshold_scale &  abs(Drought_04De_new.data_BS$DrySeasonLength)<= threshold_scale & abs(Drought_04De_new.data_BS$WTD)<= threshold_scale & abs(Drought_04De_new.data_BS$TreeHeight)<= threshold_scale) ,]
#calculate HAND class
Drought_04De_new.data_ha<-Drought_04De_new.data_BS
#Drought_04De_new.data_ha<-Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$WTD_ori <=40 ),]
Drought_04De_new.data_ha$HAND_CLASS=floor(Drought_04De_new.data_ha$WTD_ori)+1
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 1),]$HAND_CLASS=Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 1),]$HAND_CLASS-1
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS=(Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS%/%2+1)*2
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS== -1),]$HAND_CLASS=Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS == -1),]$HAND_CLASS-0
Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS=Drought_04De_new.data_ha[which(Drought_04De_new.data_ha$HAND_CLASS >= 0),]$HAND_CLASS-1
#Calculate mean of region/basin
Drought_04De_new.data_ha_mean<-Drought_04De_new.data_ha
#calculate prediction according to mean climate for Guiana shield
Model_fit_array<-Model_Seg_Prediction(Drought_04De_new.data_ha_mean,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=9, moment = median)
Drought_04De_new.data_ha_mean$EVI_anomaly_fit_09=Model_fit_array$.value
#calculate prediction according to mean climate forSouthern Amazon
Model_fit_array<-Model_Seg_Prediction(Drought_04De_new.data_ha_mean,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=12, moment = median)
Drought_04De_new.data_ha_mean$EVI_anomaly_fit_12=Model_fit_array$.value
#calculate prediction according to mean climate forEverWet Amazon
Model_fit_array<-Model_Seg_Prediction(Drought_04De_new.data_ha_mean,mod.IA_Guiana,x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=21, moment = median)
Drought_04De_new.data_ha_mean$EVI_anomaly_fit_21=Model_fit_array$.value
#Calculate prediction 
Drought_04De_new.data_ha<-add_fitted(Drought_04De_new.data_ha,mod.IA_Guiana)
colnames(Drought_04De_new.data_ha)[23] = 'EVI_anomaly_fit'
#corrected EVI anomaly in Guyana shield
Model_correct_array<-Model_Seg_Correction(Drought_04De_new.data_ha, Drought_04De_new.data_ha_mean, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_09', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=9)
arr_ob_Guyana<-summary_group_full(Model_correct_array,x0.var='EVI_anomaly', x1.var='EVI_anomaly_corrected', flag=1)
arr_ob_Guyana$type="Ob_Guyana"
#corrected EVI anomaly in Southern Amazon
Model_correct_array<-Model_Seg_Correction(Drought_04De_new.data_ha, Drought_04De_new.data_ha_mean, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_12', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=12)
arr_ob_Southern<-summary_group_full(Model_correct_array,x0.var='EVI_anomaly', x1.var='EVI_anomaly_corrected', flag=1)
arr_ob_Southern$type="Ob_SouthernAmazon"
#corrected EVI anomaly in EverWet Amazon
Model_correct_array<-Model_Seg_Correction(Drought_04De_new.data_ha, Drought_04De_new.data_ha_mean, y.var ='EVI_anomaly', y0.other ='EVI_anomaly_fit',y1.other ='EVI_anomaly_fit_21', x0.var = 'WTD_ori', x1.var = 'HAND_CLASS', x4.other='SegGeo_number', Seg=21)
arr_ob_Everwet<-summary_group_full(Model_correct_array,x0.var='EVI_anomaly', x1.var='EVI_anomaly_corrected', flag=1)
arr_ob_Everwet$type="Ob_EverwetAmazon"
#draw the plot
rangeMin=-1.5
range_mm=2.
rangeMax=rangeMin+range_mm
BY_increase=0.5
color_threshold=1
arr_ob_Everwet<-arr_ob_Everwet[which(arr_ob_Everwet$HAND_CLASS >=0 &  arr_ob_Everwet$HAND_CLASS<=21)  ,]
arr_ori<-rbind( arr_ob_Guyana, arr_ob_Southern)
arr_ori2<-rbind( arr_ob_Everwet)
arr_ori<-arr_ori[which(arr_ori$HAND_CLASS >=0 &  arr_ori$HAND_CLASS<=40)  ,]
cbPalette <- c( "#D53E4F", "#238B45","#3288BD")# "#3288BD", "#D53E4F", "#238B45"  blue  red green
Figure4_PanelE_ob<-ggplot(arr_ori, aes(x=HAND_CLASS, y=EVI_anomaly, group=factor(type)))+ 
  geom_hline(yintercept = 0,width=.6 ) + #linetype="dashed"
  geom_errorbar(data=arr_ori2,aes(ymin=EVI_anomaly-ci, ymax=EVI_anomaly+ci),colour=cbPalette[3], size=0.5, width=.3) +
  geom_errorbar(aes(ymin=EVI_anomaly-ci, ymax=EVI_anomaly+ci,colour = factor(type),group=factor(type)), size=0.5, width=.3) +  #,face="bold"
  geom_point( size=3.5, shape=23,aes(group=factor(type),colour=factor(type),fill=factor(type))) + #,fill=factor(Brando_summary_Geo.type)
    scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  geom_point(data=arr_ori2, aes(x=HAND_CLASS, y=EVI_anomaly),colour=cbPalette[3], shape=23, size=3.5, na.rm=TRUE)+
   xlab("HAND (meter)") +
  ylab("EVI anomaly") +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) +
  scale_x_continuous(limits =c(0,40), breaks=seq(0,40,  by=10))+
  theme_bw() +
  theme( legend.position=c(1,0))+
  theme(legend.justification=c(1,0))+
  theme_classic() + 
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=24))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=24))+  #,face="bold"
  theme(axis.title=element_text(size=24))+ #,face="bold"
  theme(plot.title=element_text(size=24))
Figure4_PanelE_ob
#--------------------------------------------------------------------------------------------------------------
##-------------------------------- Figure 4 Panel E prediction continued---------------------------------------
#calculate prediction for each region
#ever-wet Amazon
Model_fit_array<-Model_Seg_Prediction_Region(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=21, moment = median)
Pred_EverWet<-Model_fit_array
Pred_EverWet$type="Prediction_EverWet"
#southern Amazon
Model_fit_array<-Model_Seg_Prediction_Region(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=12, moment = median)
Pred_Southern<-Model_fit_array
Pred_Southern$type="Prediction_Southern"
#Guyana Shield
Model_fit_array<-Model_Seg_Prediction_Region(Drought_04De_new.data_BS,mod.IA_Guiana,WTD_arr, x0.ori = 'WTD_ori',x0.var = 'WTD', x1.var = 'SoilFertility',  x2.var = 'SoilSand_content',  x3.var = 'TreeHeight',  x4.var = 'DrySeasonLength', y.var ='EVI_anomaly', x0.other='PAR_anomaly',x1.other='VPD_anomaly',x2.other='MCWD_anomaly',x3.other='Pre_anomaly', x4.other='SegGeo_number',Seg=09, moment = median)
Pred_Guyana<-Model_fit_array
Pred_Guyana$type="Prediction_Guyana"
#Draw figures
rangeMin=-1.5
range_mm=2.0
rangeMax=rangeMin+range_mm
BY_increase=0.5
Pred_EverWet<-Pred_EverWet[which(Pred_EverWet$WTD_reverse >=0 &  Pred_EverWet$WTD_reverse<=25)  ,]
tgca_combine<-rbind(Pred_EverWet,Pred_Southern,Pred_Guyana) #
tgca_combine<-tgca_combine[which(tgca_combine$WTD_reverse >=0 &  tgca_combine$WTD_reverse<=40)  ,]
cbPalette <- c(   "#3288BD", "#D53E4F", "#238B45") 
Figure4_PanelE_P<-ggplot(tgca_combine, aes(x = WTD_reverse, y =fit, colour = factor(type),group=factor(type))) + 
  theme_bw() +
  scale_y_continuous(limits =c(rangeMin,rangeMax), breaks=seq(rangeMin,rangeMax,  by=BY_increase)) + 
  #  geom_line(linetype = "dashed") +
  geom_hline(yintercept = 0,  width=.6) +
  geom_smooth( data=subset(tgca_combine,type=="Prediction_Southern" |type=="Prediction_EverWet"|type=="Prediction_Guyana"),aes(ymin = lci, ymax = uci,group=factor(type),colour=factor(type),fill=factor(type)), stat = "identity", linetype = "dashed") +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette)  +
  xlab("HAND (meter)") +
  ylab("EVI anomaly") +
  scale_x_continuous(breaks=seq(0,40,10))+ theme_few() %+replace%
  theme(panel.background = element_rect(fill = NA,colour = "black", 
                                        size =1.0))+
  theme(axis.text.x = element_text( color="Black", size=18))+  #,face="bold")
  theme(axis.text.y = element_text( color="Black", size=18))+  #,face="bold"
  theme(axis.title=element_text(size=18))+ #,face="bold"
  theme(plot.title=element_text(size=18,face="bold"))

Figure4_PanelE_P
#merge figures
Figure4_PanelE<-merge_figures_sameunits(Figure4_PanelE_P,Figure4_PanelE_ob)
grid.draw(Figure4_PanelE)
##-------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------functions-------------------------------------------------

##--------------------------------------------------------------------------------------------------------------
