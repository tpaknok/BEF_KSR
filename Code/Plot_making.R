library(tidyverse)
library(ggplot2)

typeI_df <- rbind(result_summary0_20,result_summary0.5_20,result_summary1_20,
                  result_summary0_40,result_summary0.5_40,result_summary1_40,
                  result_summary0_60,result_summary0.5_60,result_summary1_60)
typeI_df <- as.data.frame(typeI_df)
typeI_df$b1 <- 0

plot_typeI <- typeI_df %>%
  pivot_longer(cols=ends_with("_sig"),names_to="Model",values_to="Sig_count") %>%
  group_by(nspp,true_lambda,Model) %>%
  summarize(typeI = sum(Sig_count)/500) %>%
  filter(Model != "m_best_sig") %>%
  mutate(ID = interaction(Model,true_lambda)) %>%
  filter(ID != "m_original_sig.1") %>%
  mutate(Model = fct_recode(Model,
                            "True VCV" = "m_true_sig" ,
                            "Brownian motion" = "m_original_sig" ,
                            "Optimized VCV" = "m_optim_sig" ,
                            "Linear regression" = "m_without_comp_sig")
  ) %>%
  mutate(Model = fct_relevel(Model,"True VCV","Optimized VCV","Brownian motion","Linear regression")) %>%
  mutate(true_lambda = as.factor(true_lambda)) %>%
  mutate(true_lambda = fct_recode(true_lambda,
                                  "λ~`=`~0" = "0",
                                  "λ~`=`~0.5" = "0.5",
                                  "λ~`=`~1" = "1"
  )
  )

p_typeI <- ggplot(plot_typeI,aes(y=typeI*100,x=nspp))+
  geom_hline(yintercept=5)+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Type I error (%)")+
  xlab("Species pool size")+
  scale_x_continuous(breaks=c(20,40,60))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  ylim(0,60)+
  facet_grid(~true_lambda,labeller = label_parsed)+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

plot(p_typeI)

ggsave(("Figure/p_typeI.tiff"),width=18,height=6,dpi=600,compression="lzw")
### Type II error
library(tidyverse)
typeII_df <- rbind(result_summary0_20_b1_0.5,result_summary0_20_b1_1,
                   result_summary0_40_b1_0.5,result_summary0_40_b1_1,
                   result_summary0_60_b1_0.5,result_summary0_60_b1_1,
                   result_summary0.5_20_b1_0.5,result_summary0.5_20_b1_1,
                   result_summary0.5_40_b1_0.5,result_summary0.5_40_b1_1,
                   result_summary0.5_60_b1_0.5,result_summary0.5_60_b1_1,
                   result_summary1_20_b1_0.5,result_summary1_20_b1_1,
                   result_summary1_40_b1_0.5,result_summary1_40_b1_1,
                   result_summary1_60_b1_0.5,result_summary1_60_b1_1
                   )

typeII_df <- as.data.frame(typeII_df)
typeII_df$b1 <- rep(c(rep(0.5,500),rep(1,500)),9)

plot_typeII <- typeII_df %>%
  pivot_longer(cols=ends_with("_sig"),names_to="Model",values_to="Sig_count") %>%
  group_by(nspp,true_lambda,b1,Model) %>%
  summarize(typeII = sum(Sig_count)/500) %>%
  filter(Model != "m_best_sig") %>%
  mutate(ID = interaction(Model,true_lambda)) %>%
  filter(ID != "m_original_sig.1") %>%
  mutate(Model = fct_recode(Model,
                            "True VCV" = "m_true_sig" ,
                            "Brownian motion" = "m_original_sig" ,
                            "Optimized VCV" = "m_optim_sig" ,
                            "Linear regression" = "m_without_comp_sig"
  )
  ) %>%
  mutate(true_lambda = as.factor(true_lambda),
         b1 = as.factor(b1)) %>%
  mutate(Model = fct_relevel(Model,"True VCV","Optimized VCV","Brownian motion","Linear regression")) %>%
  mutate(true_lambda = fct_recode(true_lambda,
                             "λ ~`=`~0" = "0",
                             "λ~`=`~0.5" = "0.5",
                             "λ~`=`~1" = "1"),
         b1 = fct_recode(b1,
                         "b[1]~`=`~0.5" = "0.5",
                         "b[1]~`=`~1" = "1"),
  )

p_typeII <- ggplot(plot_typeII, aes(y=100-typeII*100,x=nspp))+
  geom_line(aes(group=Model,colour=Model))+
  ylab("Type II error (%)")+
  xlab("Species pool size")+
  scale_x_continuous(breaks=c(20,40,60))+
  geom_point(aes(group=Model,colour=Model))+
  scale_colour_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"))+
  facet_grid(b1~true_lambda,labeller = label_parsed)+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


plot(p_typeII)

ggsave(("Figure/p_typeII.tiff"),width=18,height=18,dpi=600,units="cm",compression="lzw")

### slope estimate

slope_df <- rbind(typeI_df,typeII_df)
plot_slope <- slope_df %>%
  pivot_longer(cols=ends_with("_slope"),names_to="Model",values_to="Slope") %>%
  select("Slope","true_lambda","nspp","Model","b1") %>%
  filter(Model != "m_best_slope") %>%
  mutate(ID = interaction(Model,true_lambda)) %>%
  filter(ID != "m_original_slope.1") %>%
  mutate(Model = fct_recode(Model,
                            "True VCV" = "m_true_slope" ,
                            "Brownian motion" = "m_original_slope" ,
                            "Optimized VCV" = "m_optim_slope" ,
                            "Linear regression" = "m_without_comp_slope"
  )
  ) %>%
  mutate(true_lambda = as.factor(true_lambda),
         b1_factor = as.factor(b1)) %>%
  mutate(Model = fct_relevel(Model,"True VCV","Optimized VCV","Brownian motion","Linear regression")) %>%
  mutate(true_lambda = fct_recode(true_lambda,
                                  "λ ~`=`~0" = "0",
                                  "λ~`=`~0.5" = "0.5",
                                  "λ~`=`~1" = "1"),
         b1_factor = fct_recode(b1_factor,
                         "b[1]~`=`~0" = "0",
                         "b[1]~`=`~0.5" = "0.5",
                         "b[1]~`=`~1" = "1"),
  )

p_slope <- ggplot(plot_slope, aes(x=as.factor(nspp),y=Slope))+
  geom_violin(aes(group=interaction(Model,nspp),fill=interaction(Model)),position="dodge")+
  geom_hline(aes(yintercept=b1))+
  ylab(bquote(b[1]~mean~estimate))+
  xlab("Species pool size")+
  facet_grid(b1_factor~true_lambda,labeller = label_parsed,scales="free")+
  scale_fill_manual(values=c("#009E73","#CC79A7","#E69F00","#0072B2"),name="Model")+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


plot(p_slope)
ggsave(("Figure/slope.tiff"),width=18,height=18,dpi=600,units="cm",compression="lzw")

### lambda signal
lambda_df <- rbind(typeI_df,typeII_df)
plot_lambda <- slope_df %>%
  select("V11","true_lambda","nspp","b1") %>%
  mutate(b1_factor = as.factor(b1),
         true_lambda_factor = as.factor(true_lambda)) %>%
  mutate(true_lambda_factor = fct_recode(true_lambda_factor,
                                  "λ ~`=`~0" = "0",
                                  "λ~`=`~0.5" = "0.5",
                                  "λ~`=`~1" = "1"),
         b1_factor = fct_recode(b1_factor,
                                "b[1]~`=`~0" = "0",
                                "b[1]~`=`~0.5" = "0.5",
                                "b[1]~`=`~1" = "1"),
  )


p_lambda <- ggplot(plot_lambda, aes(x=as.factor(nspp),y=V11))+
  geom_violin()+
  geom_hline(aes(yintercept=true_lambda))+
  ylab("λ estimates")+
  xlab("Species pool size")+
  facet_grid(b1_factor~true_lambda_factor,labeller = label_parsed,scales="free")+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position="bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


plot(p_lambda)
