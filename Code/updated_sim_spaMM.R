library(phytools)
library(tidyverse)
library(CPR)
library(ape)

setwd("C:/Users/pakno/OneDrive - University of Toronto/BEF_KSR")
data(KSR)
data(KSR_MLtree)
data(KSR_EF)

VCV_sp <- vcv(KSR_MLtree)
head(VCV_sp)

#comm <- KSR/rowSums(KSR)
spaMM_formula <-  y~x+corrMatrix(1|comp_id)

###
result_summary <- list()

time1 <- Sys.time()
lambda_true <- c(0,0.5,1)
nspp <- 40
sim <- 500
b1 <- 0.5
  for (i in 1:3) {
    set.seed(999)

    message(i)
    # result <- lapply(1:500,function(x) {message(x)
    #   BEF_simulate(comm = NULL,
    #                VCV_sp = NULL,
    #                scale=1,
    #                nspp=50,
    #                nsite=88,
    #                min_richness=1,
    #                max_richness=15,
    #                spaMM_formula=spaMM_formula,
    #                b1=0,
    #                signals_X="phy_cor",
    #                signals_Y = T,
    #                intercept = 0,
    #                y_mean = 0,
    #                y_sd = 1,
    #                x_mean=0,
    #                x_sd = 1,
    #                noise_mean = 0,
    #                noise_sd = 1,
    #                lambda_true= lambda_true[[i]])
    #   }
    #   )

    result <- list()
    for (j in 1:sim) {
      message(j," ",nspp)
      result[[j]] <- BEF_simulate(comm = NULL,
                                  VCV_sp = NULL,
                                  scale=1,
                                  nspp=nspp,
                                  nsite=88,
                                  min_richness=1,
                                  max_richness= 4,
                                  spaMM_formula=spaMM_formula,
                                  b1=b1,
                                  signals_X="phy_cor",
                                  signals_Y = T,
                                  intercept = 0,
                                  y_mean = 0,
                                  y_sd = 1,
                                  x_mean=0,
                                  x_sd = 1,
                                  noise_mean = 0,
                                  noise_sd = 1,
                                  lambda_true= lambda_true[[i]])

    }

    result <- do.call(rbind,result)
    file_name <- paste0("result_summary",lambda_true[[i]],"_",nspp,"_b1_",b1)
    assign(file_name,result)
  }
time2 <- Sys.time()
time2-time1
###

result_summary <- list()
lambda_true = c(0,0.5,1)
for (i in 1:3) {
  result_summary_2 <- list()
  set.seed(9999)

  message(i)
  result <- lapply(1:1000,function(x) {message(x)
    BEF_simulate(comm = comm,
                 VCV_sp = NULL,
                 scale=1,
                 spaMM_formula,
                 b1=0.5,
                 signals_X="phy_cor",
                 signals_Y = T,
                 intercept = 0,
                 y_mean = 0,
                 y_sd = 1,
                 x_mean=0,
                 x_sd = 1,
                 noise_mean = 0,
                 noise_sd = 1,
                 lambda_true= lambda_true[[i]])
  })
  result_summary[[i]] <- BEF_simulate_eval(result)
  result_summary_2 <- result_summary[[i]]
  save(result_summary_2,file=paste0("Sim_result/spaMM_eval_lambda_",lambda_true[[i]],"_2.Rdata"))
  remove(result_summary_2)

}
