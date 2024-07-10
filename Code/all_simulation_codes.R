library(phytools)
library(tidyverse)
library(CPR)
data(KSR)
data(KSR_MLtree)
data(KSR_EF)

VCV_sp <- vcv(KSR_MLtree)
head(VCV_sp)

###
sim <- 500
optim_sim <- 100
optim_seq <- sample(c(rep(TRUE, optim_sim), rep(FALSE, sim-optim_sim)), sim ,replace = F)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
GLS_formula <- y~x
###

sim_data_lambda_1 <- BEF_simulate(comm = KSR,
                                  V = VCV_sp,
                                  ef_mean = 0,
                                  sd = 1,
                                  b1 = 0,
                                  signals_X = "sr",
                                  signals_intercept = T,
                                  signals_slope = F,
                                  lambda_true = 1,
                                  sim = sim)


GLS_lambda_1 <- lapply(1:sim, function(x) CPR_GLS(formula = GLS_formula,
                                                  df = sim_data_lambda_1$sim_dat[[x]],
                                                  VCV_sp = VCV_sp,
                                                  optim.lambda = optim_seq[[x]],
                                                  comm = KSR))

GLS_eval_lambda_1 <- BEF_simulate_eval(GLS_lambda_1,type="GLS")
save(GLS_eval_lambda_1,file="Sim_result/GLS_eval_lambda_1.Rdata")
remove(GLS_lambda_1)

INLA_lambda_1 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                               priors=NULL,
                                               df = sim_data_lambda_1$sim_dat[[x]],
                                               VCV_sp = VCV_sp,
                                               comm=KSR,
                                               family= "gaussian",
                                               inla.rerun=3,
                                               optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_1 <- BEF_simulate_eval(INLA_lambda_1,type="INLA")
save(INLA_eval_lambda_1,file="Sim_result/INLA_eval_lambda_1.Rdata")
remove(INLA_lambda_1)

###
sim_data_lambda_0 <- BEF_simulate(comm = KSR,
                                  V = VCV_sp,
                                  ef_mean = 0,
                                  sd = 1,
                                  b1 = 0,
                                  signals_X = "sr",
                                  signals_intercept = T,
                                  signals_slope = F,
                                  lambda_true = 0,
                                  sim = sim)


VCV_sp_lambda0 <- VCV_sp*0
diag(VCV_sp_lambda0) <- diag(VCV_sp)

GLS_lambda_0 <- lapply(1:sim, function(x) CPR_GLS(formula = y~x,
                                                  df = sim_data_lambda_0$sim_dat[[x]],
                                                  VCV_sp = VCV_sp_lambda0,
                                                  optim.lambda = optim_seq[[x]],
                                                  comm = KSR))

GLS_eval_lambda_0 <- BEF_simulate_eval(GLS_lambda_0,type="GLS")
save(GLS_eval_lambda_0,file="Sim_result/GLS_eval_lambda_0.Rdata")
remove(GLS_lambda_0)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
INLA_lambda_0 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                               priors=NULL,
                                               df = sim_data_lambda_0$sim_dat[[x]],
                                               VCV_sp = VCV_sp_lambda0,
                                               comm=KSR,
                                               family= "gaussian",
                                               inla.rerun=3,
                                               optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_0 <- BEF_simulate_eval(INLA_lambda_0,type="INLA")
save(INLA_eval_lambda_0,file="Sim_result/INLA_eval_lambda_0.Rdata")
remove(INLA_lambda_0)

###
sim_data_lambda_0.5 <- BEF_simulate(comm = KSR,
                                    V = VCV_sp,
                                    ef_mean = 0,
                                    sd = 1,
                                    b1 = 0,
                                    signals_X = "sr",
                                    signals_intercept = T,
                                    signals_slope = F,
                                    lambda_true = 0.5,
                                    sim = sim)


VCV_sp_lambda0.5 <- VCV_sp*0.5
diag(VCV_sp_lambda0.5) <- diag(VCV_sp)

GLS_lambda_0.5 <- lapply(1:sim, function(x) CPR_GLS(formula = y~x,
                                                    df = sim_data_lambda_0.5$sim_dat[[x]],
                                                    VCV_sp = VCV_sp_lambda0.5,
                                                    optim.lambda = optim_seq[[x]],
                                                    comm = KSR))

GLS_eval_lambda_0.5 <- BEF_simulate_eval(GLS_lambda_0.5,type="GLS")
save(GLS_eval_lambda_0.5,file="Sim_result/GLS_eval_lambda_05.Rdata")
remove(GLS_lambda_0.5)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
INLA_lambda_0.5 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                                 priors=NULL,
                                                 df = sim_data_lambda_0.5$sim_dat[[x]],
                                                 VCV_sp = VCV_sp_lambda0.5,
                                                 comm=KSR,
                                                 family= "gaussian",
                                                 inla.rerun=3,
                                                 optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_0.5 <- BEF_simulate_eval(INLA_lambda_0.5,type="INLA")
save(INLA_eval_lambda_0.5,file="Sim_result/INLA_eval_lambda_05.Rdata")
remove(INLA_lambda_0.5)


###
sim_data_no_lambda <- BEF_simulate(comm = KSR,
                                   V = VCV_sp,
                                   ef_mean = 0,
                                   sd = 1,
                                   b1 = 0,
                                   signals_X = "sr",
                                   signals_intercept = F,
                                   signals_slope = F,
                                   lambda_true = 0,
                                   sim = sim)


GLS_lambda_no_lambda <- lapply(1:sim, function(x) CPR_GLS(formula = y~x,
                                                          df = sim_data_no_lambda$sim_dat[[x]],
                                                          VCV_sp = VCV_sp,
                                                          optim.lambda = optim_seq[[x]],
                                                          comm = KSR))

GLS_eval_lambda_no_lambda <- BEF_simulate_eval(GLS_lambda_no_lambda,type="GLS")
save(GLS_eval_lambda_no_lambda,file="Sim_result/GLS_eval_lambda_no_lambda.Rdata")
remove(GLS_lambda_no_lambda)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
INLA_lambda_no_lambda <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                                       priors=NULL,
                                                       df = sim_data_no_lambda$sim_dat[[x]],
                                                       VCV_sp = VCV_sp,
                                                       comm=KSR,
                                                       family= "gaussian",
                                                       inla.rerun=3,
                                                       optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_no_lambda <- BEF_simulate_eval(INLA_lambda_no_lambda,type="INLA")
save(INLA_eval_lambda_no_lambda,file="Sim_result/INLA_eval_lambda_no_lambda.Rdata")
remove(INLA_lambda_no_lambda)

##################### Type-II error

sim_data_lambda_1_2 <- BEF_simulate(comm = KSR,
                                  V = VCV_sp,
                                  ef_mean = 0,
                                  sd = 1,
                                  b1 = 0.25,
                                  signals_X = "sr",
                                  signals_intercept = T,
                                  signals_slope = F,
                                  lambda_true = 1,
                                  sim = sim)


GLS_lambda_1_2 <- lapply(1:sim, function(x) CPR_GLS(formula = GLS_formula,
                                                  df = sim_data_lambda_1_2$sim_dat[[x]],
                                                  VCV_sp = VCV_sp,
                                                  optim.lambda = optim_seq[[x]],
                                                  comm = KSR))


GLS_eval_lambda_1_2 <- BEF_simulate_eval(GLS_lambda_1_2,type="GLS")
save(GLS_eval_lambda_1_2,file="Sim_result/GLS_eval_lambda_1_2.Rdata")
remove(GLS_lambda_1_2)

INLA_lambda_1_2 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                               priors=NULL,
                                               df = sim_data_lambda_1_2$sim_dat[[x]],
                                               VCV_sp = VCV_sp,
                                               comm=KSR,
                                               family= "gaussian",
                                               inla.rerun=3,
                                               optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_1_2 <- BEF_simulate_eval(INLA_lambda_1_2,type="INLA")
save(INLA_eval_lambda_1_2,file="Sim_result/INLA_eval_lambda_1_2.Rdata")
remove(INLA_lambda_1_2)

###
sim_data_lambda_0_2 <- BEF_simulate(comm = KSR,
                                  V = VCV_sp,
                                  ef_mean = 0,
                                  sd = 1,
                                  b1 = 0.25,
                                  signals_X = "sr",
                                  signals_intercept = T,
                                  signals_slope = F,
                                  lambda_true = 0,
                                  sim = sim)


VCV_sp_lambda0 <- VCV_sp*0
diag(VCV_sp_lambda0) <- diag(VCV_sp)

GLS_lambda_0_2 <- lapply(1:sim, function(x) CPR_GLS(formula = y~x,
                                                  df = sim_data_lambda_0_2$sim_dat[[x]],
                                                  VCV_sp = VCV_sp_lambda0,
                                                  optim.lambda = optim_seq[[x]],
                                                  comm = KSR))

GLS_eval_lambda_0_2 <- BEF_simulate_eval(GLS_lambda_0_2,type="GLS")
save(GLS_eval_lambda_0_2,file="Sim_result/GLS_eval_lambda_0_2.Rdata")
remove(GLS_lambda_0_2)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
INLA_lambda_0_2 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                               priors=NULL,
                                               df = sim_data_lambda_0_2$sim_dat[[x]],
                                               VCV_sp = VCV_sp_lambda0,
                                               comm=KSR,
                                               family= "gaussian",
                                               inla.rerun=3,
                                               optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_0_2 <- BEF_simulate_eval(INLA_lambda_0_2,type="INLA")
save(INLA_eval_lambda_0_2,file="Sim_result/INLA_eval_lambda_0_2.Rdata")
remove(INLA_lambda_0_2)


###
sim_data_lambda_0.5_2 <- BEF_simulate(comm = KSR,
                                    V = VCV_sp,
                                    ef_mean = 0,
                                    sd = 1,
                                    b1 = 0.25,
                                    signals_X = "sr",
                                    signals_intercept = T,
                                    signals_slope = F,
                                    lambda_true = 0.5,
                                    sim = sim)


VCV_sp_lambda0.5 <- VCV_sp*0.5
diag(VCV_sp_lambda0.5) <- diag(VCV_sp)

GLS_lambda_0.5_2 <- lapply(1:sim, function(x) CPR_GLS(formula = y~x,
                                                    df = sim_data_lambda_0.5_2$sim_dat[[x]],
                                                    VCV_sp = VCV_sp_lambda0.5,
                                                    optim.lambda = optim_seq[[x]],
                                                    comm = KSR))

GLS_eval_lambda_0.5_2 <- BEF_simulate_eval(GLS_lambda_0.5_2,type="GLS")
save(GLS_eval_lambda_0.5_2,file="Sim_result/GLS_eval_lambda_05_2.Rdata")
remove(GLS_lambda_0.5_2)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
INLA_lambda_0.5_2 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                                 priors=NULL,
                                                 df = sim_data_lambda_0.5_2$sim_dat[[x]],
                                                 VCV_sp = VCV_sp_lambda0.5,
                                                 comm=KSR,
                                                 family= "gaussian",
                                                 inla.rerun=3,
                                                 optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_0.5_2 <- BEF_simulate_eval(INLA_lambda_0.5_2,type="INLA")
save(INLA_eval_lambda_0.5_2,file="Sim_result/INLA_eval_lambda_05_2.Rdata")
remove(INLA_lambda_0.5_2)

###
sim_data_no_lambda_2 <- BEF_simulate(comm = KSR,
                                   V = VCV_sp,
                                   ef_mean = 0,
                                   sd = 1,
                                   b1 = 0.25,
                                   signals_X = "sr",
                                   signals_intercept = F,
                                   signals_slope = F,
                                   lambda_true = 0,
                                   sim = sim)


GLS_lambda_no_lambda_2 <- lapply(1:sim, function(x) CPR_GLS(formula = y~x,
                                                          df = sim_data_no_lambda_2$sim_dat[[x]],
                                                          VCV_sp = VCV_sp,
                                                          optim.lambda = optim_seq[[x]],
                                                          comm = KSR))


GLS_eval_lambda_no_lambda_2 <- BEF_simulate_eval(GLS_lambda_no_lambda_2,type="GLS")
save(GLS_eval_lambda_no_lambda_2,file="Sim_result/GLS_eval_lambda_no_lambda_2.Rdata")
remove(GLS_lambda_no_lambda_2)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)
INLA_lambda_no_lambda_2 <- lapply(1:sim, function(x) CPR(formula = INLA_formula,
                                                       priors=NULL,
                                                       df = sim_data_no_lambda_2$sim_dat[[x]],
                                                       VCV_sp = VCV_sp,
                                                       comm=KSR,
                                                       family= "gaussian",
                                                       inla.rerun=3,
                                                       optim.lambda = optim_seq[[x]]))

INLA_eval_lambda_no_lambda_2 <- BEF_simulate_eval(INLA_lambda_no_lambda_2,type="INLA")
save(INLA_eval_lambda_no_lambda_2,file="Sim_result/INLA_eval_lambda_no_lambda_2.Rdata")
remove(INLA_lambda_no_lambda_2)
