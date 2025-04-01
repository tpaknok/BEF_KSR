###This is the code for reproducing the results in the manuscript
### For general usage, see the website of EcoCoMix https://tpaknok.github.io/EcoCoMix/articles/Empirical_single.html

library(phytools)
library(tidyverse)
library(EcoCoMix)
library(ape)

data(KSR)
data(KSR_MLtree)
data(KSR_EF)

### a function to get prediction from spaMM models
get_predictions <- function(object) {
  newdata <- data.frame(Real.rich=unique(object$best_model$data$Real.rich))
  predict_with_phylo <- predict(object$best_model,newdata=newdata,re.form=NA,variances=list(fixefVar=TRUE),intervals="fixefVar",binding="Response")
  predict_without_phylo <- predict(object$without_comp_model,newdata=newdata,re.form=NA,variances=list(fixefVar=TRUE),intervals="fixefVar",binding="Response")

  predict_with_phylo <- cbind(predict_with_phylo,
                              attr(predict_with_phylo,"intervals"),
                              data.frame(Name=attr(object$best_model$main_terms_info$Y,"respname")),
                              data.frame(Model="Best model"),
                              data.frame(Sig=ifelse(object$best_model_satt$`Pr(>F)`[[1]] < 0.05,"Sig","Insig"))
                              )

  predict_without_phylo <- cbind(predict_without_phylo,
                                 attr(predict_without_phylo,"intervals"),
                                 data.frame(Name=attr(object$best_model$main_terms_info$Y,"respname")),
                                 data.frame(Model="Linear regression"),
                                 data.frame(Sig=ifelse(summary(object$without_comp_model,details=T,verbose=F)$beta_table[2,"p-value"] < 0.05,"Sig","Insig"))
                                 )

  predict_df <- rbind(predict_with_phylo,predict_without_phylo)
}

### get the phylogenetic vcv matrix
KSR_EF$mpd <- picante::mpd(KSR,cophenetic(KSR_MLtree))
KSR_EF[is.na(KSR_EF$mpd),"mpd"] <- 0

KSR_EF$log_bugs <- log(KSR_EF$bugs+1)
KSR_EF$log_flwr_total <- log(KSR_EF$flwr_total+1)
KSR_EF$log_bug_rich <- log(KSR_EF$bug.rich+1)
KSR_EF$log_poll_total <- log(KSR_EF$poll_total+1)

### a for loop to analyze all functions
resp <- c("litter2012","ave.biomass","LAI","mean.N.change","log_poll_total","log_flwr_total",
                  "Mass.loss.2month","Damage_effect","log_bugs","log_bug_rich")

result_df <- predict_df_sr <- predict_df_mpd <- NULL
for (i in 1:10) {
  message(i)
  y <- KSR_EF[,resp[[i]]]
  m_sr <- EcoCoMix(y~Real.rich+corrMatrix(1|comp_id),
                                     data=KSR_EF,
                                     comm=KSR,
                                     VCV_sp = vcv(KSR_MLtree),
                                     method.spaMM="REML",
                                     init=list())

  assign(paste0("m_",resp[[i]]),m_sr)

  m_lm_sr <- lm(y~Real.rich,data=KSR_EF) #do a linear regression for comparison

  result <- data.frame(MM_p_sr=ifelse(length(ranef(m_sr$best_model)) == 0,
                                      m_sr$best_model_satt[,5], #satterthwaite's method p value in best_model_satt
                                      m_sr$best_model_satt[,6]),
                       without_re_p_sr = m_sr$without_comp_anova[1,5], #without compositional random effect
                       AIC_diff_sr = m_sr$AIC[[1]]-m_sr$AIC[[3]],
                       R2m_sr = ifelse(length(ranef(m_sr$best_model)) == 0, NA, get_R2(m_sr$best_model)[[1]]),
                       R2c_sr = ifelse(length(ranef(m_sr$best_model)) == 0, NA, get_R2(m_sr$best_model)[[2]]),
                       R2_sr_lm = summary(m_lm_sr)$r.sq,
                       optim_lambda_full_sr = m_sr$optimized_lambda, # lambda of full model
                       optim_lambda_int = m_sr$optimized_lambda_int, #lambda of intercept-only model
                       resp=resp[[i]])

  result_df <- rbind(result_df,result)

  predict_df_sr <- rbind(predict_df_sr,cbind(get_predictions(m_sr),resp=resp[[i]]))
}

### Visualize the results
predict_df <- predict_df_sr %>%
  mutate(Name = as.factor(resp)) %>%
  mutate(Name = fct_recode(Name,
                            "Biomass" = "ave.biomass",
                            "Structural complexity" = "LAI" ,
                            "Litter" = "litter2012" ,
                            "Damage reduction" = "Damage_effect",
                            "log(arthropod richness+1)" = "log_bug_rich",
                            "log(arthropod abundance+1)" = "log_bugs",
                            "log(flower production+1)" = "log_flwr_total",
                            "log(pollinator abundance+1)" = "log_poll_total",
                            "Decomposition" = "Mass.loss.2month",
                            "Soil nitrogen (Delta N)" =  "mean.N.change")
  ) %>%
  mutate(Name = fct_relevel(Name,
                            "Biomass",
                            "Structural complexity",
                            "Damage reduction",
                            "Decomposition",
                            "Litter",
                            "log(arthropod richness+1)",
                            "log(arthropod abundance+1)",
                            "log(flower production+1)",
                            "log(pollinator abundance+1)",
                            "Soil nitrogen (Delta N)"))

plot_data <- KSR_EF %>%
  select(Real.rich,litter2012,ave.biomass,LAI,mean.N.change,Mass.loss.2month,Damage_effect,log_bugs,log_flwr_total,log_bug_rich,log_poll_total) %>%
  pivot_longer(cols=litter2012:log_poll_total,names_to="Name",values_to="Value") %>%
  mutate(Name = fct_recode(Name,
                           "Biomass" = "ave.biomass",
                           "Structural complexity" = "LAI" ,
                           "Litter" = "litter2012" ,
                           "Damage reduction" = "Damage_effect",
                           "log(arthropod richness+1)" = "log_bug_rich",
                           "log(arthropod abundance+1)" = "log_bugs",
                           "log(flower production+1)" = "log_flwr_total",
                           "log(pollinator abundance+1)" = "log_poll_total",
                           "Decomposition" = "Mass.loss.2month",
                           "Soil nitrogen (Delta N)" =  "mean.N.change")) %>%
  filter(!(Real.rich == 1 & Name == "Damage reduction"))

predict_df_subset <- subset(predict_df,Name == "Decomposition")
plot_data_subset <- subset(plot_data,Name == "Decomposition")

summary(lm(KSR_EF$Mass.loss.2month~KSR_EF$Real.rich))

label_subset<- data.frame(Real.rich=c(4,4),
                          Response = c(4.5,4.35),
                          labels = c("Best~model:~R[m]^2 ~`=`~0.03~`;`~R[c]^2 ~`=`~0.39","Linear~regression:~R^2~`=`~0.07"))
p_KSR <- ggplot(predict_df_subset,aes(x=Real.rich,y=Response))+
  geom_point(data=plot_data_subset,aes(x=Real.rich,y=Value,group=Real.rich),position=position_dodge2(width=0.1))+
  geom_line(aes(colour=Model,linetype=Sig))+
  geom_ribbon(aes(ymin=fixefVar_0.025,ymax=fixefVar_0.975,fill=Model),alpha=0.25)+
  geom_text(data=label_subset, aes(x=Real.rich,y=Response,label=labels),hjust=0.95,vjust=-3.5,parse=T)+
  xlab("Species Richness")+
  ylab("Mass loss after 2 months (g)")+
  ylim(2.7,5.2)+
  scale_x_continuous(breaks=c(1,2,3,4))+
  scale_fill_manual(values=c("#CC79A7","#0072B2"))+
  scale_linetype_manual(values=c(2,1),guide="none")+
  theme_classic()+
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.position="bottom",
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        strip.text = element_text(size=11))
plot(p_KSR)
ggsave(("Figure/p_KSR.tiff"),width=11,height=11,dpi=600,units="cm",compression="lzw")

p_KSR_all <- ggplot(predict_df,aes(x=Real.rich,y=Response))+
  geom_point(data=plot_data,aes(x=Real.rich,y=Value,group=Real.rich),position=position_dodge2(width=0.1))+
  geom_line(aes(colour=Model,linetype=Sig))+
  geom_ribbon(aes(ymin=fixefVar_0.025,ymax=fixefVar_0.975,fill=Model),alpha=0.25)+
  xlab("Species Richness")+
  ylab("Value")+
  scale_x_continuous(breaks=c(1,2,3,4))+
  scale_fill_manual(values=c("#CC79A7","#0072B2"))+
  scale_linetype_manual(values=c(2,1),guide="none")+
  facet_wrap(~Name,ncol=2,scales="free")+
  theme_bw()+
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.position="bottom",
        legend.text = element_text(size=11),
        legend.title = element_text(size=11),
        strip.text = element_text(size=11))
plot(p_KSR_all)

ggsave(("Figure/p_KSR_all.tiff"),width=18,height=18,dpi=600,units="cm",compression="lzw")

###
c(m_LAI$optimized_lambda_int,
  m_litter2012$optimized_lambda_int,
  m_mean.N.change$optimized_lambda_int,
  m_Damage_effect$optimized_lambda_int,
  m_ave.biomass$optimized_lambda_int,
  m_Mass.loss.2month$optimized_lambda_int,
  m_log_bug_rich$optimized_lambda_int,
  m_log_bugs$optimized_lambda_int,
  m_log_poll_total$optimized_lambda_int,
  m_log_flwr_total$optimized_lambda_int)

m_LAI$AIC[[1]]-m_LAI$AIC[[3]]
m_litter2012$AIC[[1]]-m_litter2012$AIC[[3]]
m_mean.N.change$AIC[[1]]-m_mean.N.change$AIC[[3]]
m_Damage_effect$AIC[[1]]-m_Damage_effect$AIC[[3]]
m_Mass.loss.2month$AIC[[1]]-m_Mass.loss.2month$AIC[[3]]
m_ave.biomass$AIC[[1]]-m_ave.biomass$AIC[[3]]
m_log_bug_rich$AIC[[1]]-m_log_bug_rich$AIC[[3]]
m_log_bugs$AIC[[1]]-m_log_bugs$AIC[[3]]
m_log_poll_total$AIC[[1]]-m_log_poll_total$AIC[[3]]
m_log_flwr_total$AIC[[1]]-m_log_flwr_total$AIC[[3]]

get_R2(m_LAI$best_model)
get_R2(m_Damage_effect$best_model)
get_R2(m_ave.biomass$best_model)
get_R2(m_Mass.loss.2month$best_model)
get_R2(m_log_bug_rich$best_model)
get_R2(m_log_bugs$best_model)
get_R2(m_log_poll_total$best_model)
get_R2(m_log_flwr_total$best_model)

m_LAI$AIC[[7]]-m_LAI$AIC[[6]]
m_litter2012$AIC[[7]]-m_litter2012$AIC[[6]]
m_mean.N.change$AIC[[7]]-m_mean.N.change$AIC[[6]]
m_Damage_effect$AIC[[7]]-m_Damage_effect$AIC[[6]]
m_Mass.loss.2month$AIC[[7]]-m_Mass.loss.2month$AIC[[6]]
m_ave.biomass$AIC[[7]]-m_ave.biomass$AIC[[6]]
m_log_bug_rich$AIC[[7]]-m_log_bug_rich$AIC[[6]]
m_log_bugs$AIC[[7]]-m_log_bugs$AIC[[6]]
m_log_poll_total$AIC[[7]]-m_log_poll_total$AIC[[6]]
m_log_flwr_total$AIC[[7]]-m_log_flwr_total$AIC[[6]]

### Table S1

#As an example, summary(m_ave.biomass$optimized_lambda_model,verbose=F)$beta_table[2,1:2] extract the slope and SE
#as.data.frame(m_ave.biomass$optimized_lambda_model_satt)[5:6] extract the Satterwaite's method p-value and F-value
#Then get the R2 (for lm) or R2m and R2c for mixed model
#Finally, the lambda (based on intercept-only mixed model) and the cAIC of each model.
table_df <- list(c("Biomass",summary(m_ave.biomass$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_ave.biomass$optimized_lambda_model_satt)[5:6],NA,get_R2(m_ave.biomass$optimized_lambda_model)[[1]],get_R2(m_ave.biomass$optimized_lambda_model)[[2]],m_ave.biomass$optimized_lambda_int,m_ave.biomass$AIC[[3]]),
                 c("Biomass",summary(m_ave.biomass$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_ave.biomass$without_comp_anova)[1,4:5],get_R2(m_ave.biomass$without_comp_model)[[1]],NA,NA,NA,m_ave.biomass$AIC[[1]]),
                 c("LAI",summary(m_LAI$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_LAI$optimized_lambda_model_satt)[5:6],NA,get_R2(m_LAI$optimized_lambda_model)[[1]],get_R2(m_LAI$optimized_lambda_model)[[2]],m_LAI$optimized_lambda_int,m_LAI$AIC[[3]]),
                 c("LAI",summary(m_LAI$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_LAI$without_comp_anova)[1,4:5],get_R2(m_LAI$without_comp_model)[[1]],NA,NA,NA,m_LAI$AIC[[1]]),
                 c("Damage_effect",summary(m_Damage_effect$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_Damage_effect$optimized_lambda_model_satt)[5:6],NA,get_R2(m_Damage_effect$optimized_lambda_model)[[1]],get_R2(m_Damage_effect$optimized_lambda_model)[[2]],m_Damage_effect$optimized_lambda_int,m_Damage_effect$AIC[[3]]),
                 c("Damage_effect",summary(m_Damage_effect$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_Damage_effect$without_comp_anova)[1,4:5],get_R2(m_Damage_effect$without_comp_model)[[1]],NA,NA,NA,m_Damage_effect$AIC[[1]]),
                 c("Mass.loss.2month",summary(m_Mass.loss.2month$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_Mass.loss.2month$optimized_lambda_model_satt)[5:6],NA,get_R2(m_Mass.loss.2month$optimized_lambda_model)[[1]],get_R2(m_Mass.loss.2month$optimized_lambda_model)[[2]],m_Mass.loss.2month$optimized_lambda_int,m_Mass.loss.2month$AIC[[3]]),
                 c("Mass.loss.2month",summary(m_Mass.loss.2month$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_Mass.loss.2month$without_comp_anova)[1,4:5],get_R2(m_Mass.loss.2month$without_comp_model)[[1]],NA,NA,NA,m_Mass.loss.2month$AIC[[1]]),
                 c("log_bug_rich",summary(m_log_bug_rich$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_bug_rich$optimized_lambda_model_satt)[5:6],NA,get_R2(m_log_bug_rich$optimized_lambda_model)[[1]],get_R2(m_log_bug_rich$optimized_lambda_model)[[2]],m_log_bug_rich$optimized_lambda_int,m_log_bug_rich$AIC[[3]]),
                 c("log_bug_rich",summary(m_log_bug_rich$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_bug_rich$without_comp_anova)[1,4:5],get_R2(m_log_bug_rich$without_comp_model)[[1]],NA,NA,NA,m_log_bug_rich$AIC[[1]]),
                 c("log_bugs",summary(m_log_bugs$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_bugs$optimized_lambda_model_satt)[5:6],NA,get_R2(m_log_bugs$optimized_lambda_model)[[1]],get_R2(m_log_bugs$optimized_lambda_model)[[2]],m_log_bugs$optimized_lambda_int,m_log_bugs$AIC[[3]]),
                 c("log_bugs",summary(m_log_bugs$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_bugs$without_comp_anova)[1,4:5],get_R2(m_log_bugs$without_comp_model)[[1]],NA,NA,NA,m_log_bugs$AIC[[1]]),
                 c("log_poll_total",summary(m_log_poll_total$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_poll_total$optimized_lambda_model_satt)[5:6],NA,get_R2(m_log_poll_total$optimized_lambda_model)[[1]],get_R2(m_log_poll_total$optimized_lambda_model)[[2]],m_log_poll_total$optimized_lambda_int,m_log_poll_total$AIC[[3]]),
                 c("log_poll_total",summary(m_log_poll_total$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_poll_total$without_comp_anova)[1,4:5],get_R2(m_log_poll_total$without_comp_model)[[1]],NA,NA,NA,m_log_poll_total$AIC[[1]]),
                 c("log_flwr_total",summary(m_log_flwr_total$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_flwr_total$optimized_lambda_model_satt)[5:6],NA,get_R2(m_log_flwr_total$optimized_lambda_model)[[1]],get_R2(m_log_flwr_total$optimized_lambda_model)[[2]],m_log_flwr_total$optimized_lambda_int,m_log_flwr_total$AIC[[3]]),
                 c("log_flwr_total",summary(m_log_flwr_total$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_log_flwr_total$without_comp_anova)[1,4:5],get_R2(m_log_flwr_total$without_comp_model)[[1]],NA,NA,NA,m_log_flwr_total$AIC[[1]]),
                 c("litter2012",summary(m_litter2012$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_litter2012$optimized_lambda_model_satt)[5:6],NA,get_R2(m_litter2012$optimized_lambda_model)[[1]],get_R2(m_litter2012$optimized_lambda_model)[[2]],m_litter2012$optimized_lambda_int,m_litter2012$AIC[[3]]),
                 c("litter2012",summary(m_litter2012$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_litter2012$without_comp_anova)[1,4:5],get_R2(m_litter2012$without_comp_model)[[1]],NA,NA,NA,m_litter2012$AIC[[1]]),
                 c("mean.N.change",summary(m_mean.N.change$optimized_lambda_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_mean.N.change$optimized_lambda_model_satt)[5:6],NA,get_R2(m_mean.N.change$optimized_lambda_model)[[1]],get_R2(m_mean.N.change$optimized_lambda_model)[[2]],m_mean.N.change$optimized_lambda_int,m_mean.N.change$AIC[[3]]),
                 c("mean.N.change",summary(m_mean.N.change$without_comp_model,verbose=F)$beta_table[2,1:2],as.data.frame(m_mean.N.change$without_comp_anova)[1,4:5],get_R2(m_mean.N.change$without_comp_model)[[1]],NA,NA,NA,m_mean.N.change$AIC[[1]])
)

table_df <- do.call(rbind,lapply(table_df,unlist))

write.csv(table_df,"Table/table_df.csv")

###


