comm <- data.frame(A=c(1,0,0,0,0,1,1,0,1,0,0,1,0,0,0),
                   B=c(0,1,0,0,0,1,0,1,0,1,0,0,1,0,0),
                   C=c(0,0,1,0,0,0,1,1,0,0,1,0,0,1,0),
                   D=c(0,0,0,1,0,0,0,0,1,1,1,0,0,0,1),
                   E=c(0,0,0,0,1,0,0,0,0,0,0,1,1,1,1))
tree <- starTree(c("A","B","C","D","E"),1)

#############

comm <- data.frame(A=c(1,0,0,0,1,1,0,1,0,0),
                   B=c(0,1,0,0,1,0,1,0,1,0),
                   C=c(0,0,1,0,0,1,1,0,0,1),
                   D=c(0,0,0,1,0,0,0,1,1,1))
tree <- starTree(c("A","B","C","D"),1)
#################

comm <- data.frame(A=c(1,0,0,1,1,0),
                   B=c(0,1,0,1,0,1),
                   C=c(0,0,1,0,1,1))


library(phytools)
tree <- starTree(c("A","B","C"),1)
#############

comm <- do.call("rbind", replicate(5, comm, simplify = FALSE))

VCV_sp <- vcv(tree)
diag(VCV_sp) <- 1

library(CPR)
#Fig 1c seems that when there are few comms it doesn't work really well...
sim_data <- BEF_simulate(comm,
             VCV_sp,
             0,
             10,
             b1=0,
             signals_X="sr",
             signals_intercept=T,
             lambda_true=0,
             sim=100,
             seed=1000) #note that ef_mean must be zero.

#If ef_mean follow NxPA+NYPA and PA is non-zero, that is "indirectly" contributing to a species richness effect.
#if non-zero constant, that would mean Pa can change across communities (all composition should have the same expected mean).
#e.g. PX = PX+PY, PX=PY, PX = 2PX =, PX = 0 the only solution
# sd assumed to be the same to fulfull normality assumption

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)

results<- lapply(1:100, function(x) CPR(formula = INLA_formula,
                              priors=NULL,
                              df = sim_data$sim_dat[[x]],
                              VCV_sp = VCV_sp,
                              comm=comm,
                              family= "gaussian",
                              inla.rerun=3,
                              optim.lambda = F))

BEF_simulate_eval(results,type="INLA")

df <-  as.data.frame(do.call(rbind,lapply(1:20,function(x) sign(results[[x]]$without_phylo_model[2,3])*sign(results[[x]]$without_phylo_model[2,5]))))
df$true <- unlist(do.call(rbind,lapply(1:20,function(x) sign(results[[x]]$original_VCV_model[2,3])*sign(results[[x]]$original_VCV_model[2,5]))))

df$diff <- df[,1]-df[,2]

GLS_formula <- y~x

GLS <- lapply(1:500, function(x) CPR_GLS(formula = GLS_formula,
                                                  df = sim_data$sim_dat[[x]],
                                                  VCV_sp = VCV_sp,
                                                  optim.lambda = F,
                                                  comm = comm))

BEF_simulate_eval(GLS,type="GLS")

df <-  as.data.frame(do.call(rbind,lapply(1:100,function(x) GLS[[x]]$without_phylo_model[2,4])))
df$true <- unlist(do.call(rbind,lapply(1:100,function(x) GLS[[x]]$original_VCV_model[2,4])))

sum(df$V1 < 0.05)
sum(df$true < 0.05)

df_raw <- sim_data$sim_dat[[12]]
plot(df_raw$x,df_raw$y)

results1 <- CPR(formula = INLA_formula,
                                        priors=NULL,
                                        df = sim_data$sim_dat[[12]],
                                        VCV_sp = VCV_sp,
                                        comm=comm,
                                        family= "gaussian",
                                        inla.rerun=3,
                                        optim.lambda = F)
results1$without_phylo_model
results1$original_VCV_model

#Fig 1b

comm <- data.frame(A=c(1,0,0,1,1,0,1),
                   B=c(0,1,0,1,0,1,1),
                   C=c(0,0,1,0,1,1,1))

sim_data <- BEF_simulate(comm,
                         VCV_sp,
                         0,
                         1,
                         b1=0.75,
                         signals_X="sr",
                         signals_intercept=T,
                         lambda_true=0,
                         sim=100,
                         seed=1000)

INLA_formula <-  y~x+f(comm,model="generic0",Cmatrix=Phylo)

results1c<- lapply(1:100, function(x) CPR(formula = INLA_formula,
                                        priors=NULL,
                                        df = sim_data$sim_dat[[x]],
                                        VCV_sp = VCV_sp,
                                        comm=comm,
                                        family= "gaussian",
                                        inla.rerun=3,
                                        optim.lambda = F))

BEF_simulate_eval(results1c,type="INLA")

df <-  as.data.frame(do.call(rbind,lapply(1:10,function(x) sign(results1c[[x]]$without_phylo_model[2,3])*sign(results1c[[x]]$without_phylo_model[2,5]))))
df$true <- unlist(do.call(rbind,lapply(1:10,function(x) sign(results1c[[x]]$original_VCV_model[2,3])*sign(results1c[[x]]$original_VCV_model[2,5]))))

df$diff <- df[,1]-df[,2]

results1 <- CPR(formula = INLA_formula,
                priors=NULL,
                df = sim_data$sim_dat[[9]],
                VCV_sp = VCV_sp,
                comm=comm,
                family= "gaussian",
                inla.rerun=3,
                optim.lambda = F)

results1$without_phylo_model
results1$original_VCV_model

df_raw <- sim_data$sim_dat[[2]] #9.19,39,41,45
plot(df_raw$x,df_raw$y)
