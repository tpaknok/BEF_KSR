### This script is for reproducing the simulation results described in the main text
### Tsang, T. P. N. & Cadotte, M. W. (2025). Species overlap and phylogenetic relatedness result in community statistical non-independence (and what to do about it). Ecology Letters.

library(phytools)
library(tidyverse)
library(EcoCoMix)
library(ape)

### Simulations - Initial settings
set.seed(123)
nspp <- c(14,28,42,56) #species pool size
sim <- 500 #number of iterations to be conducted
b1 <- c(0,0.25) #True effect of SR
lambda_true <- runif(sim) #500 true lambdas

count <- 1 #number of iterations conducted (for checking progress)

result_df <- NULL

spaMM_formula <-  y~x1+corrMatrix(1|comp_id) #formula used for the regression. See the documentation of ?fitme in EcoCoMix

### a for loop for simulations based on different scenarios.
### This one takes a long time! So an R object containing the simulation results (sim500.Rdata) has been provided.
### You can go to L67 if you don't want to run the simulation.

for (k in 1:length(b1)) { #slope
  for (i in 1:length(nspp)) { #species pool size
      for (l in 1:sim) { #number of simulations

        message("sim_",l,"_",nspp[[i]],"_",b1[[k]])

        result  <- tryCatch(BEF_simulate(comm = NULL,
                                         VCV_sp = NULL,
                                         nspp=nspp[[i]],
                                         nsite=88,
                                         min_richness=1,
                                         max_richness= 4, #local species richness = 1-4
                                         spaMM_formula=spaMM_formula,
                                         b1=b1[[k]],
                                         signals_X="sr", #species richness as the predictor
                                         noise_mean = 0,
                                         noise_sd = 0.01, #very low noise
                                         lambda_true= lambda_true[[l]],
                                         conv_fail_drop = T, #drop runs with failed convergence
                                         scale_all=F, #no need to scale the predictor
                                         optim.lambda=T,
                                         init=list(),
                                         int_model=F,
                                         method.spaMM = "REML",
                                         control.optim=list(factr=1e12)),
                            error = function(e) e)

        result_df <- rbind(result,result_df)

        print(result_df %>%
           dplyr::select(b1,nspp,m_optim_sig,m_true_sig,m_original_sig,m_best_sig,m_without_comp_sig) %>%
           pivot_longer(!b1:nspp,names_to="sig") %>%
           group_by(b1,nspp,sig) %>%
           summarize(sig_count = sum(value))) #summarize the results

        count <- count+1
        print(count)
        }
      }
    }

### you can also load sim500.Rdata and then run the script below
# load("./Data/sim500.RData") loading the result
all_result <- as.data.frame(do.call(rbind,result))

summary_stat <- result_df %>%
  dplyr::select(b1,nspp,m_optim_sig,m_true_sig,m_original_sig,m_best_sig,m_without_comp_sig) %>%
  pivot_longer(!b1:nspp,names_to="Model") %>%
  group_by(b1,nspp,Model) %>%
  summarize(sig_count = sum(value)/sim)

### Making fig 1 - subpanel for type-I error
dir.create("Figure")

typeI_df <- summary_stat %>%
  filter(b1 == 0 & Model != "m_best_sig") %>%
  mutate(Model = fct_recode(Model,
                            "True model" = "m_true_sig" ,
                            "Brownian motion" = "m_original_sig" ,
                            "Optimized model" = "m_optim_sig" ,
                            "Linear regression" = "m_without_comp_sig")
  ) %>%
  mutate(Model = fct_relevel(Model,"True model","Optimized model","Brownian motion","Linear regression"))

p_typeI <- ggplot(typeI_df,aes(y=sig_count*100,x=nspp))+
  geom_hline(yintercept=5)+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Type I error (%)")+
  xlab("Species pool size")+
  scale_x_continuous(breaks=c(14,28,42,56))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(A)",hjust=-0.25,vjust=1.2,size=5.5)+
  ylim(0,85)+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

plot(p_typeI)

### Making fig 1 - subpanel for power
power_df <- summary_stat %>%
  filter(b1 == 0.25 & Model != "m_best_sig") %>%
  mutate(Model = fct_recode(Model,
                            "True model" = "m_true_sig" ,
                            "Brownian motion" = "m_original_sig" ,
                            "Optimized model" = "m_optim_sig" ,
                            "Linear regression" = "m_without_comp_sig")
  ) %>%
  mutate(Model = fct_relevel(Model,"True model","Optimized model","Brownian motion","Linear regression"))

p_power <- ggplot(power_df,aes(y=sig_count*100,x=nspp))+
  geom_hline(yintercept = 80)+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Power (%)")+
  xlab("")+
  scale_x_continuous(breaks=c(14,28,42,56))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(B)",hjust=-0.25,vjust=1.2,size=5.5)+
  ylim(75,105)+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

plot(p_power)

### Making fig 1 - subpanel for rmse and mean error
coef_df_b1_0 <- result_df %>%
  dplyr::select(b1,nspp,m_optim_slope,m_true_slope,m_original_slope,m_best_slope,m_without_comp_slope) %>%
  pivot_longer(!b1:nspp,names_to="Model") %>%
  filter(b1 == 0 & Model != "m_best_slope") %>%
  mutate(Model = fct_recode(Model,
                            "True model" = "m_true_slope" ,
                            "Brownian motion" = "m_original_slope" ,
                            "Optimized model" = "m_optim_slope" ,
                            "Linear regression" = "m_without_comp_slope")
  ) %>%
  mutate(Model = fct_relevel(Model,"True model","Optimized model","Brownian motion","Linear regression"))

library(see)

Error_df_b1_0 <- coef_df_b1_0 %>%
  group_by(nspp,Model) %>%
  summarize(rmse = sqrt(mean((value-b1)^2)),
            me = mean(value-b1))

p_ME_b1_0 <- ggplot(Error_df_b1_0,aes(y=me,x=nspp))+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Bias (Mean Error)")+
  xlab("")+
  scale_x_continuous(breaks=c(14,28,42,56))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(C)",hjust=-0.25,vjust=1.2,size=5.5)+
  ylim(-0.02,0.02)+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

p_RMSE_b1_0 <- ggplot(Error_df_b1_0,aes(y=rmse,x=nspp))+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Accuracy (RMSE)")+
  xlab("")+
  scale_x_continuous(breaks=c(14,28,42,56))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(E)",hjust=-0.25,vjust=1.2,size=5.5)+
  theme_bw()+
  ylim(0,0.15)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

coef_df_b1_0.25 <- result_df %>%
  dplyr::select(b1,nspp,m_optim_slope,m_true_slope,m_original_slope,m_best_slope,m_without_comp_slope) %>%
  pivot_longer(!b1:nspp,names_to="Model") %>%
  filter(b1 == 0.25 & Model != "m_best_slope") %>%
  mutate(Model = fct_recode(Model,
                            "True model" = "m_true_slope" ,
                            "Brownian motion" = "m_original_slope" ,
                            "Optimized model" = "m_optim_slope" ,
                            "Linear regression" = "m_without_comp_slope")
  ) %>%
  mutate(Model = fct_relevel(Model,"True model","Optimized model","Brownian motion","Linear regression"))

Error_df_b1_0.25 <- coef_df_b1_0.25 %>%
  group_by(nspp,Model) %>%
  summarize(rmse = sqrt(mean((value-b1)^2)),
            me = mean(value-b1))

p_ME_b1_0.25 <- ggplot(Error_df_b1_0.25,aes(y=me,x=nspp))+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Bias (Mean Error)")+
  xlab("")+
  scale_x_continuous(breaks=c(14,28,42,56))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(D)",hjust=-0.25,vjust=1.2,size=5.5)+
  theme_bw()+
  ylim(-0.02,0.02)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

p_RMSE_b1_0.25 <- ggplot(Error_df_b1_0.25,aes(y=rmse,x=nspp))+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Accuracy (RMSE)")+
  xlab("")+
  scale_x_continuous(breaks=c(14,28,42,56))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(F)",hjust=-0.25,vjust=1.2,size=5.5)+
  theme_bw()+
  ylim(0,0.15)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

plot(p_RMSE_b1_0.25)

###combining them into one figure
library(ggpubr)

ggarrange(p_typeI,p_power,
          p_ME_b1_0,p_ME_b1_0.25,
          p_RMSE_b1_0,p_RMSE_b1_0.25,
          nrow=3,ncol=2,common.legend=T,legend="bottom")

ggsave(("./Figure/p_sim.tiff"),width=17,height=17,dpi=600,units="cm",compression="lzw",bg="white")

### visualizing coef estimates across simulations
p_coef_b1_0 <- ggplot(coef_df_b1_0,aes(y=value,x=nspp))+
  geom_hline(yintercept = 0)+
  geom_violinhalf(aes(group=interaction(Model,nspp),fill=Model),position=position_dodge(width=4),scale="width",
                  linewidth=0.2)+
  scale_x_continuous(breaks=c(14,28,42,56))+
  scale_fill_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(A)",hjust=-0.25,vjust=1.2,size=5.5)+
  ylab(bquote(β[SR]~estimates))+
  xlab("")+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

plot(p_coef_b1_0)

p_coef_b1_0.25 <- ggplot(coef_df_b1_0.25,aes(y=value,x=nspp))+
  geom_hline(yintercept = 0.25)+
  geom_violinhalf(aes(group=interaction(Model,nspp),fill=Model),position=position_dodge(width=6),scale="width",
                  linewidth=0.2)+
  scale_x_continuous(breaks=c(14,28,42,56))+
  scale_fill_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  annotate("text",x=-Inf,y=Inf,label="(B)",hjust=-0.25,vjust=1.2,size=5.5)+
  ylab(bquote(β[SR]~estimates))+
  xlab("")+
  labs(fill="Model")+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text = element_text(size=12))

plot(p_coef_b1_0.25)

ggarrange(p_coef_b1_0,p_coef_b1_0.25,nrow=1,ncol=2,common.legend=T,legend="bottom")
ggsave(("./Figure/p_coef.tiff"),width=24,height=12,dpi=600,units="cm",compression="lzw",bg="white")

