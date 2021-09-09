################## model summary
################## 06/13/2020
###### summary for posterior inference of lysozyme dataset
############################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("foreach", "doParallel", 
              "boot", "parallel","rjags", 
              "ggplot2", "dplyr", "invgamma", "R2jags", "MCMCvis", "R2WinBUGS",
              "ggmcmc", "ggpubr", "xtable", "scales","AICcmodavg")
ipak(packages)
###############################################
experimental.conditions <- data.frame(NaCl=c(250,125,125,
                                             75,50,300,
                                             150,75,300,
                                             150,
                                             100,
                                             100,
                                             50,200), 
                                      pH=c(6.9,6.9,4.7,4.7,4.7,6.9,6.9,6.9,4.7,
                                           4.7,
                                           4.7,
                                           6.9,6.9,6.9)) # NaCl: mM
# # ################################
###############################
#### 06/22/2020
#### dimer, monomer, M12, Mass
###############################
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior.rds") 
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior.rds")
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_dimer <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_dimer.rds") 
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_Mass <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_Mass.rds") 
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_M12 <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_M12.rds") 

# pull out the DIC
AICcmodavg::DIC(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior)
# AICcmodavg::DIC(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior)
AICcmodavg::DIC(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_dimer)
AICcmodavg::DIC(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_Mass)
AICcmodavg::DIC(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_M12)

##########################
# setwd("/Users/fan/Documents/UCI/Research/light_scattering_new/lysozyme/")
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior.rds")
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior.rds")
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment <- readRDS("./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment.rds")

# lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior.rds
############# monomer 
# monomer, loose
use.experimental.conditions <- (experimental.conditions[1:14,])[c(1,10,11,12,13,14, 2:9),]
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior = as.mcmc(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior)
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_A2_df <- data.frame(A2 = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior[[1]][,1:14]),
                                                                                                  pH = rep(use.experimental.conditions$pH, rep(400,14)),
                                                                                                  NaCl = rep(use.experimental.conditions$NaCl, rep(400,14)))

lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_A2_df  %>% ggplot(aes(x=factor(NaCl),y=A2)) + geom_boxplot(outlier.size = 0) + facet_wrap(~pH) + xlab("NaCl (mM)") + ylab("A2")

############
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_A2_df %>% group_by(pH, NaCl) %>% summarise(m_A2 = scientific(mean(A2), digits=3),
                                                                                                                         sd_A2 = scientific(sd(A2),digits=3),
                                                                                                                         q025_A2 = scientific(quantile(A2,0.025),digits=3),
                                                                                                                         q0975 = scientific(quantile(A2, 0.975), digits=3)) %>% xtable()

lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_A2_df %>% group_by(pH, NaCl) %>% summarise(m_A2 = mean(A2) * 10^5,
                                                                                                                         sd_A2 = sd(A2) * 10^5,
                                                                                                                         q025_A2 = quantile(A2,0.025) * 10^5,
                                                                                                                         q0975 = quantile(A2, 0.975)*10^5,
                                                                                                                         prob_positive = mean(A2>0)) %>% xtable(digits=3)
# lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_full_loose_prior[[1]][,"s2n"] %>% summary()
# monomer, diffuse
use.experimental.conditions <- (experimental.conditions[1:14,])[c(1,10,11,12,13,14, 2:9),]
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior <- as.mcmc(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior)


lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior_A2_df <- data.frame(A2 = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior[[1]][,1:14]),
                                                                                                       pH = rep(use.experimental.conditions$pH, rep(400,14)),
                                                                                                       NaCl = rep(use.experimental.conditions$NaCl, rep(400,14)))

lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior_A2_df  %>% ggplot(aes(x=factor(NaCl),y=A2)) + geom_boxplot(outlier.size = 0) + facet_wrap(~pH) + xlab("NaCl (mM)") + ylab("A2")

# monomer, no adjustment 
use.experimental.conditions <- (experimental.conditions[1:14,])[c(1,10,11,12,13,14, 2:9),]
lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment <- as.mcmc(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment)


lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment_A2_df <- data.frame(A2 = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment[[1]][,1:14]),
                                                                                                    pH = rep(use.experimental.conditions$pH, rep(400,14)),
                                                                                                    NaCl = rep(use.experimental.conditions$NaCl, rep(400,14)))


lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment_A2_df %>% group_by(pH, NaCl) %>% summarise(m_A2 = mean(A2) * 10^5,
                                                                                                                         sd_A2 = sd(A2) * 10^5,
                                                                                                                         q025_A2 = quantile(A2,0.025) * 10^5,
                                                                                                                         q0975 = quantile(A2, 0.975)*10^5,
                                                                                                                         prob_positive = mean(A2>0)) %>% xtable(digits=3)

# lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_full_diffuse_prior[[1]][,"s2n"] %>% summary()
# combine them together
lysozyme_complete_pool_constant_dndc_noconmeas_all_A2_df <- data.frame(rbind(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_A2_df,
                                                                             lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior_A2_df,
                                                                             lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment_A2_df),
                                                                       model = rep(c("M1", "M1-loose-prior", "No adjustment"),
                                                                                   rep(400*14,3)) )
#####  
# https://stackoverflow.com/questions/45618181/change-font-size-of-titles-from-facet-wrap/45618363
lysozyme_complete_pool_constant_dndc_noconmeas_all_A2_df %>% ggplot(aes(x=factor(NaCl),
                                                                        y=A2, fill=model)) + geom_boxplot(outlier.size = 0) + facet_wrap(~pH, labeller = label_both) + xlab("NaCl (mM)") + ylab(expression(A[~"2"])) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                              axis.text.x = element_text(angle = 0, hjust=1),
                                                              strip.text = element_text(size = 20)) 



##### dndc
lysozyme_complete_pool_constant_dndc_noconmeas_all_dndc_df <- data.frame(dndc = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior[[1]][,"dndc"],
                                                                                  lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior[[1]][,"dndc"],
                                                                                  lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment[[1]][,"dndc"]),
                                                                         model = rep(c("M1", "M1-loose-prior", "No adjustment"),
                                                                                     rep(400,3)))

p_dndc <- lysozyme_complete_pool_constant_dndc_noconmeas_all_dndc_df %>% ggplot(aes(x=model, y=dndc)) + geom_boxplot(outlier.size = 0) + xlab("") + ylab("") + ggtitle("dndc") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                                axis.text.x = element_text(angle = 30, hjust=1, size=20),
                                                                                                                                                                                                                                                                                                                                                                strip.text = element_text(size = 20)) 


p_dndc_density <- lysozyme_complete_pool_constant_dndc_noconmeas_all_dndc_df %>% ggplot(aes(x=dndc, fill=model)) + geom_density() + xlab("") + ylab("") + ggtitle("dndc") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                          axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                                                                                                                                                                                          strip.text = element_text(size = 20)) 


p_dndc_density
#### s2R
lysozyme_complete_pool_constant_dndc_noconmeas_all_s2R_df <- data.frame(s2R = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior[[1]][,"s2R"],
                                                                                  lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior[[1]][,"s2R"],
                                                                                lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment[[1]][,"s2R"]),
                                                                         model = rep(c("M1", "M1-loose-prior", "No adjustment"),
                                                                                     rep(400,3)))


p_s2R <- lysozyme_complete_pool_constant_dndc_noconmeas_all_s2R_df %>% ggplot(aes(x=model,y=sqrt(s2R))) + geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("") + ggtitle(expression(sigma[~"R"])) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                       axis.text.x = element_text(angle = 30, hjust=1,size=20),
                                                                                                                                                                                                                       strip.text = element_text(size = 20)) 
p_s2R
# expression in R
# https://stackoverflow.com/questions/10156417/subscripts-in-plots-in-r
#### s2n
lysozyme_complete_pool_constant_dndc_noconmeas_all_s2n_df <- data.frame(s2n = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior[[1]][,"s2n"],
                                                                                lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior[[1]][,"s2n"],
                                                                                lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment[[1]][,"s2n"]),
                                                                        model = rep(c("M1", "M1-loose-prior", "No adjustment"),
                                                                                    rep(400,3)))


p_s2n <- lysozyme_complete_pool_constant_dndc_noconmeas_all_s2n_df %>% ggplot(aes(x=model,y=sqrt(s2n))) + geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("") + ggtitle(expression(sigma[Delta~"n"])) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                       axis.text.x = element_text(angle = 30, hjust=1, size=20),
                                                                                                                                                                                                                       strip.text = element_text(size = 20)) 
p_s2n

#### s2c
lysozyme_complete_pool_constant_dndc_noconmeas_all_s2c_df <- data.frame(s2c = c(lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior[[1]][,"s2c"],
                                                                                lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_diffuse_prior[[1]][,"s2c"]
                                                                                ),
                                                                        model = rep(c("M1", "M1-loose-prior"),
                                                                                    rep(400,2)))


p_s2c <- lysozyme_complete_pool_constant_dndc_noconmeas_all_s2c_df %>% ggplot(aes(x = model, y=sqrt(s2c))) + geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("") + ggtitle(expression(sigma[~"u"])) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                            axis.text.x = element_text(angle = 30, hjust=1,size=20),
                                                                                                                                                                                                                                                            strip.text = element_text(size = 20)) 
p_s2c

#####
# par(mfrow=c(2,2))
figure <- ggarrange(p_dndc, p_s2R, p_s2n, p_s2c,
                    ncol = 2, nrow = 2)
figure




