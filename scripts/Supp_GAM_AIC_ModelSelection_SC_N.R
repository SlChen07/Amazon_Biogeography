
##Script may take 3-5 days to run on a machine with 64GB RAM, Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz processor.
## Extend Date Table 1 (a),(b),(c)
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

## Extend Date Table 1 (d)
##------------------------------------Model selection-----------------------------------------------------------
#considering the complexity and non-linearity of fitting, users might choose appropriate parameters of k.  
mod.AmzBiogeograhy<-mgcv::gam(EVI_anomaly ~   s(WTD,k=5,bs='tp') +s(SoilFertility,k=5,bs='tp')+ti(WTD,SoilFertility,k=c(3,3),bs='tp') 
                              +s(TreeHeight,k=5)+ti(WTD,TreeHeight,k=c(3,3),bs='tp') +s(DrySeasonLength,bs='tp')+ti(WTD,DrySeasonLength,k=c(3,3),bs='tp')
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
