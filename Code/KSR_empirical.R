library(phytools)
library(tidyverse)
library(CPR)
library(ape)

setwd("C:/Users/pakno/OneDrive - University of Toronto/BEF_KSR")
data(KSR)
data(KSR_MLtree)
data(KSR_EF)

###ave.biomass, LAI, Damage_effect, Mass.loss.2month
m <- CPR_spaMM(LAI~Real.rich+corrMatrix(1|comp_id),
          data=KSR_EF,
          comm=KSR,
          VCV_sp = vcv(KSR_MLtree),
          method.spaMM="REML")

m$best_model
summary(m$without_comp_model,details=T)
m$best_model_satt
m

m <- CPR_spaMM(litter2012~Real.rich+corrMatrix(1|comp_id),
               data=KSR_EF,
               comm=KSR,
               VCV_sp = vcv(KSR_MLtree),
               method.spaMM="REML",
               control.optim=list(factr=1e15)
               )

m$best_model
summary(m$without_comp_model,details=T)
m$best_model_satt
plot(simulateResiduals(m$optimized_lambda_model))

m <- CPR_spaMM(mean.N.change*1000~Real.rich+corrMatrix(1|comp_id), #all multiply by 1000 to minimize convergence issue
               data=KSR_EF,
               comm=KSR,
               VCV_sp = vcv(KSR_MLtree),
               method.spaMM="REML",
               control.optim=list(factr=1e15)
)

#regardless of multiplications, doesn't change cAIC conclusion (even though optim lambda changed)
m$best_model
summary(m$without_comp_model,details=T)
m$best_model_satt

##need to exclude signle diversity treatment
m <- CPR_spaMM(Damage_effect~Real.rich+corrMatrix(1|comp_id),
               data=subset(KSR_EF,Real.rich != 1),
               comm=KSR[rowSums(KSR)!=1,],
               VCV_sp = vcv(KSR_MLtree),
               method.spaMM="REML",control.optim=list(factr=1e15))

m$best_model
summary(m$without_comp_model,details=T)
m$best_model_satt
plot(simulateResiduals(m$best_model))

library(DHARMa)
## bugs, bug.rich, poll_total, flwr_total???
m2 <- CPR_spaMM(formula=bugs~Real.rich+corrMatrix(1|comp_id),
               data=KSR_EF,
               comm=KSR,
               VCV_sp = vcv(KSR_MLtree),
               method.spaMM="ML",
               family="negbin")

m2$best_model_satt
m2$without_comp_model
plot(simulateResiduals(m$best_model))

###
m3 <- CPR_spaMM(log(flwr_total+1)~Real.rich+corrMatrix(1|comp_id),
                data=KSR_EF,
                comm=KSR,
                VCV_sp = vcv(KSR_MLtree),
                method.spaMM="ML",
                control.optim=list(factr=1e15))

m3$best_model
summary(m3$without_comp_model,details=T)
m3$best_model_satt
plot(simulateResiduals(m3$best_model))

###
library(glmmTMB)
m <- glmmTMB(flwr_total~Real.rich,data=KSR_EF,family="tweedie")
summary(m)
