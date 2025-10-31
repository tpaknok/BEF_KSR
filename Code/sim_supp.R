### This script is for reproducing the simulation results for the high-diversity scenario.
### ### Tsang, T. P. N. & Cadotte, M. W. (2025). Species overlap and phylogenetic relatedness result in community statistical non-independence (and what to do about it). Ecology Letters.

library(phytools)
library(tidyverse)
library(EcoCoMix)
library(ape)
library(BBmisc)

### Simulations ----
#### Initial settings ----
set.seed(123)
nspp <- 100 #species pool size
sim <- 500  #number of iterations to be conducted
b1 <- c(0,0.25) #True effect of SR
lambda_true <- runif(sim) #500 true lambdas

count <- 1 #number of iterations conducted (for checking progress)

result_df <- NULL

spaMM_formula <-  y~x1+corrMatrix(1|comp_id)

#### For loop ----
#### a for loop for simulations based on different scenarios.
#### Skip to L77-78 if you want to load simSupp.RData directly

for (k in 1:length(b1)) { #slope
  for (i in 1:length(nspp)) { #species pool size
    for (l in 1:sim) { #number of simulations

      message("sim_",l,"_",nspp[[i]],"_",b1[[k]])
      boo <- 1

      while(boo == 1) {
      result  <- tryCatch(BEF_simulate(comm = NULL,
                                       VCV_sp = NULL,
                                       nspp=nspp[[i]],
                                       nsite=88,
                                       min_richness=10,
                                       max_richness= 50, #local richness = 10-50
                                       spaMM_formula=spaMM_formula,
                                       b1=b1[[k]],
                                       signals_X="sr", #species richness as the predictor
                                       noise_mean = 0,
                                       noise_sd = 0.01, #very small noise
                                       lambda_true= lambda_true[[l]],
                                       conv_fail_drop = T, #drop runs with failed convergence
                                       scale_all=F, #no need to scale the predictor
                                       optim.lambda=T,
                                       init=list(),
                                       int_model=F,
                                       method.spaMM = "REML",
                                       control.optim=list(factr=1e12)),
                          error = function(e) e)
         count <- count+1
         boo <- ifelse(is.error(result),1,0)

         if (boo == 1) {
           message("re-run")
         }
      }

      result_df <- rbind(result,result_df)

      print(result_df %>%
              dplyr::select(b1,nspp,m_optim_sig,m_true_sig,m_original_sig,m_best_sig,m_without_comp_sig) %>%
              pivot_longer(!b1:nspp,names_to="sig") %>%
              group_by(b1,nspp,sig) %>%
              summarize(sig_count = sum(value))) #summarize the results

      print(count)
    }
  }
}

#### Formatting the results (Table S1) ----
####load("./Data/simSupp.RData") #load the results directly

summary_stat <- result_df %>%
  dplyr::select(b1,nspp,m_optim_sig,m_true_sig,m_original_sig,m_best_sig,m_without_comp_sig) %>%
  pivot_longer(!b1:nspp,names_to="Model") %>%
  group_by(b1,nspp,Model) %>%
  summarize(sig_count = sum(value)/sim) #this produces part of Table S1 (typeI error and statistical power)

coef_df_b1 <- result_df %>%
  dplyr::select(b1,nspp,m_optim_slope,m_true_slope,m_original_slope,m_best_slope,m_without_comp_slope) %>%
  pivot_longer(!b1:nspp,names_to="Model") %>%
  filter(Model != "m_best_slope") %>%
  mutate(Model = fct_recode(Model,
                            "True model" = "m_true_slope" ,
                            "Brownian motion" = "m_original_slope" ,
                            "Optimized model" = "m_optim_slope" ,
                            "Linear regression" = "m_without_comp_slope")
  ) %>%
  mutate(Model = fct_relevel(Model,"True model","Optimized model","Brownian motion","Linear regression")) #obtain coefficient estimates in each iteration.

coef_est_df <- coef_df_b1 %>%
  mutate(diff = value-b1) %>%
  group_by(b1,Model) %>%
  summarize(me = mean(diff),
            rmse = sqrt(mean((diff)^2))) #this produces part of Table S1 (RMSE and mean error)

