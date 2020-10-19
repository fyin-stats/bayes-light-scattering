#Analysis script for the three-detector light scattering data (lysozyme experiments)
#Last modified 7/18/15 by CTB
#####################################################
# library(parallel)
# library(ggplot2)
# library(dplyr)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("foreach", "doParallel", 
              "boot", "parallel","rjags", "ggplot2", "dplyr", "invgamma", "R2jags",
              "R2WinBUGS")
ipak(packages)
mc.cores<-25
####################################
draws<-100000
thin<- 250
burn<-100000
nadapt<-100000
nchains <- 5
# data preprocess, remove the comments if you have not preprocessed the data
# source("./lysozyme/lysozyme.data.preprocess.R")
########################################
lysozyme_model <- function(){
  #mu_dndc ~ dnorm(prior.dndc.mu.mean, 1/pow(prior.dndc.mu.sd,2))
  #inv.s2dndc ~ dgamma(a_prior.dndc.s2, b_prior.dndc.s2)
  #s2dndc <- 1/inv.s2dndc
  ## variance parameters (gamma in jags, shape rate)
  inv.s2R ~ dgamma(a_R, b_R)
  inv.s2n ~ dgamma(a_n, b_n)
  s2R <- 1/inv.s2R
  s2n <- 1/inv.s2n
  
  # inv.s2cd ~ dgamma(a_cd, b_cd) 
  inv.s2c ~ dgamma(a_c, b_c)
  # s2cd <- 1/inv.s2cd
  s2c <- 1/inv.s2c
  ### 
  for(l in 1:L){
    A2[l] ~ dnorm(0, 1/pow(prior.A2.sd,2) )
  }
  
  ## stage 2
  #for(l in 1:L){
  #   dndc[l] ~ dnorm(mu_dndc, inv.s2dndc)T(0,)
  # }
  
  dndc ~ dnorm(0.197, 1/pow(0.005,2));T(0,)
  
  for(i in 1:n){
    for(l in 1:L){
      # con[i,l] ~ dlnorm( log(con.meas[i,l]), inv.s2c )
      con[i,l] <- con.meas[i,l]
    }
  }
  
  # observational model
  #for(i in 1:n){
  #  for(l in 1:L){
  #   con.meas[i,l] ~ dlnorm( log(con[i,l]), inv.s2c )
  #  }
  #}
  
  # measured refractive index of solution
  for(i in 1:nri){
    for(l in 1:L){
      ri.meas[i,l] ~ dnorm( mu_dRI[i,l], inv.s2n );T(0,)
      mu_dRI[i,l] <- dndc*con[i,l]
    }
  }
  
  # measured light scattering signal
  for(i in 1:n){
    for(j in 1:nang){
      for(l in 1:L)
      {
        scat.meas[i,j,l] ~ dnorm(mu_R[i,j,l], inv.s2R) # Wyatt, 1993; eq. 8
        mu_R[i,j,l] <- f*con[i,l]*pow(n0[l]*dndc,2)*(N+1)*monomer.mass*(1-2*A2[l]*con[i,l]*(N+1)*monomer.mass)
      }
    } 
  }
  # derived quantities
  # Assess model fit using a sum-of-squares-type discrepancy
  # measured refractive index of solution
  for(i in 1:nri){
    for(l in 1:L){
      predicted_dRI[i,l] <- mu_dRI[i,l]
      residual_dRI[i,l] <- ri.meas[i,l] - predicted_dRI[i,l]
      sq_dRI[i,l] <- pow(residual_dRI[i,l], 2)
      
      # Generate replicate data and compute fit statistics for them
      ri.meas.new[i,l] ~ dnorm(mu_dRI[i,l], inv.s2n)
      residual_dRI.new[i,l] <- ri.meas.new[i,l]-predicted_dRI[i,l]
      sq_dRI.new[i,l] <- pow(residual_dRI.new[i,l], 2)
      # ri.meas[i,l] ~ dnorm( mu_dRI[i,l], inv.s2n );T(0,)
      # mu_dRI[i,l] <- dndc*con[i,l]
    }
  }
  #####
  for(i in 1:n){
    for(j in 1:nang){
      for(l in 1:L){
        predicted_R[i,j,l] <- mu_R[i,j,l]
        residual_R[i,j,l] <- scat.meas[i,j,l] - predicted_R[i,j,l] # residuals for observed data
        sq_R[i,j,l] <- pow(residual_R[i,j,l], 2) # squared residuals
        # Generate replicate data and compute fit statistics for them
        scat.meas.new[i,j,l] ~ dnorm(mu_R[i,j,l], inv.s2R) # one new data set at each MCMC iteration
        residual_R.new[i,j,l] <- scat.meas.new[i,j,l]-predicted_R[i,j,l]
        sq_R.new[i,j,l] <- pow(residual_R.new[i,j,l], 2) # squared residuals for new data
      } 
    }
  }
  # Sum of squared residuals for actual data set and new data set
  for(l in 1:L){
    fit_R[l] <- sum(sq_R[,1,l])
    fit_R.new[l] <- sum(sq_R.new[,1,l])
    test_R[l] <- step(fit_R.new[l] - fit_R[l]) # Test whether new data more extreme
    bpval_R[l] <- mean(test_R[l])
    
    # 
    fit_dRI[l] <- sum(sq_dRI[,l])
    fit_dRI.new[l] <- sum(sq_dRI.new[,l])
    test_dRI[l] <- step(fit_dRI.new[l] - fit_dRI[l])
    bpval_dRI[l] <- mean(test_dRI[l])
  }
  
  # fit_R <- sum(sq_R[,,])
  # fit_R.new <- sum(sq_R.new[,,])
  # test_R <- step(fit_R.new - fit_R) # Test whether new data more extreme
  # bpval_R <- mean(test_R)
  # 
  # # 
  # fit_dRI <- sum(sq_dRI[,])
  # fit_dRI.new <- sum(sq_dRI.new[,])
  # test_dRI <- step(fit_dRI.new - fit_dRI)
  # bpval_dRI <- mean(test_dRI)
}
#
write.model(lysozyme_model, "./lysozyme/lysozyme_model_log_normal_no_adjustment.bug")
############################
time1 <- Sys.time()
dat_jags <- list(scat.meas = datscat_array,
                 n0.meas=n0.meas,
                 n0=solvent.base.ri.meas,
                 ri.meas=datri_matrix,
                 con.meas=con_matrix,
                 con.design=con_design_vector,
                 prior.A2.sd=1,
                 prior.dndc.mu.sd=prior.dndc.sd,
                 prior.dndc.mu.mean=prior.dndc.mean,
                 max.oligomer.size=3,
                 N=0,
                 L=length(lysozyme.df.list),
                 nri=dim(datri_matrix)[1],
                 n=dim(datscat_array)[1],
                 npred = 100,
                 f=f,
                 monomer.mass=monomer.mass,
                 nang=dim(datscat_array)[2],
                 a_prior.dndc.s2=2/2,
                 b_prior.dndc.s2=(2/2)*(0.0002^2),
                 a_R = 2/2,
                 b_R = (2/2)*(10^(-10)),
                 a_n = 2/2,
                 b_n = (2/2)*(10^(-8)),
                 a_cd = 2/2,
                 b_cd = 2/2 * prior.con.lsd^2,
                 a_c = 2/2,
                 b_c = 2/2*(log(1.05)/1.96)^2 )
# > pinvgamma(0.001^2, shape=1, rate = 0.0002^2)
# [1] 0.9607894
# a_R = ls.prsdwgt/2,
# b_R = (ls.prsdwgt/2)*ls.prsdest^2, 
# a_n = ri.prsdwgt/2,
# b_n = (ri.prsdwgt/2)*ri.prsdest^2,
# we need s2dndc <= 0.002^2 with very very high probability, thus we have inv gamma with shape 2, rate 0.0001
jags.inits <- function(){
  list("dndc"=0.197)
}
####
out <- do.call(jags.parallel, 
               list(data=dat_jags, inits=jags.inits, 
                    parameters.to.save = c("A2","dndc","s2cd","s2R", "s2c", "s2n","con", 
                                           "bpval_dRI", "bpval_R", "residual_R.new", "residual_dRI.new",
                                           "scat.meas.new", "ri.meas.new"),
                    n.chains = nchains,
                    n.iter = burn+nadapt+draws,
                    n.burnin = burn+nadapt,
                    n.thin = thin, 
                    digits = 16, 
                    model.file="./lysozyme/lysozyme_model_log_normal_no_adjustment.bug"))
# model <- jags.model(textConnection(model_string), 
#                     data = dat_jags, 
#                     inits = list(dndc = prior.dndc.mean),
#                     n.chains = nchains, 
#                     n.adapt = nadapt)
# inits = list(con=con_matrix),
# update(model, burn, progress.bar="none") # burn-in for 200000 samples
# samp <- coda.samples(model, 
#                      variable.names=c("A2","dndc","s2cd","s2R", "s2c", "s2n","con"), 
#                      n.iter=draws, thin=thin)
# saveRDS(samp, paste0("samp_jags_", i, ".rds"))
# samp
time2 <- Sys.time()
####
saveRDS(out, "./lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_no_adjustment.rds")
out_mcmc <- as.mcmc(out)
# ################# posterior predictive
# out <- readRDS("lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior.rds")
# out_mcmc <- as.mcmc(out)
# # residuals of posterior predictive draws vs measured
# R_df <- data.frame(residuals = NA,
#                    prednew = NA,
#                    observed = NA,
#                    con_level = NA,
#                    pH = NA,
#                    NaCl = NA)
# dRI_df <- data.frame(residuals=NA,
#                      prednew = NA,
#                      observed = NA, 
#                      con_level = NA,
#                      pH = NA,
#                      NaCl = NA)
# con_df <- data.frame(con = NA, 
#                      con.meas = NA,
#                      con_level = NA,
#                      pH = NA,
#                      NaCl = NA)
# #################
# L = dim(dat_jags$scat.meas)[3]
# IR = dim(dat_jags$scat.meas)[1]
# use.experimental.conditions <- experimental.conditions[1:17,][-c(11,12,14),]
# ##################
# for(l in 1:L){
#   for(i in 1:IR){
#     R_residuals_colname <- paste0("residual_R.new","[",i,",","1",",",l,"]")
#     R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
#     # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
#     R_df <- rbind(R_df, data.frame(residuals = c(out_mcmc[[1]][,R_residuals_colname]),
#                                    prednew = c(out_mcmc[[1]][,R_pred_colname]),
#                                    observed = dat_jags$scat.meas[i,1,l],
#                                    con_level = i,
#                                    pH = use.experimental.conditions$pH[l],
#                                    NaCl = use.experimental.conditions$NaCl[l]))
#   }
# }
# ####################
# for(l in 1:L){
#   for(i in 1:IR){
#     con_colname <- paste0("con","[",i,",",l,"]")
#     # R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
#     # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
#     con_df <- rbind(con_df, data.frame(con = c(out_mcmc[[1]][,con_colname]),
#                                        con.meas = dat_jags$con.meas[i,l],
#                                        con_level = i,
#                                        pH = use.experimental.conditions$pH[l],
#                                        NaCl = use.experimental.conditions$NaCl[l]))
#   }
# }
# ####################
# IdRI <- dim(dat_jags$ri.meas)[1]
# for(l in 1:L){
#   for(i in 1:IdRI){
#     dRI_residuals_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
#     dRI_pred_colname <- paste0("ri.meas.new","[",i,",",l,"]")
#     dRI_df <- rbind(dRI_df, data.frame(residuals = c(out_mcmc[[1]][,dRI_residuals_colname]),
#                                        prednew = c(out_mcmc[[1]][,dRI_pred_colname]),
#                                        observed = dat_jags$ri.meas[i,l],
#                                        con_level = i,
#                                        pH = use.experimental.conditions$pH[l],
#                                        NaCl = use.experimental.conditions$NaCl[l]))
#   }
# }
# ########################
# 
# # remove the first row
# # R_df <- R_df[-1,]
# # dRI_df <- dRI_df[-1,]
# 
# # l indexes the experiment
# # i index the concentration level
# # R_colnames <- paste0("residual_R.new","[",i,",","1",",",l,"]")
# 
# # R_df <- readRDS("/Users/fan/Documents/UCI/Research/light_scattering_new/lysozyme/R_df.rds")
# # dRI_df <- readRDS("/Users/fan/Documents/UCI/Research/light_scattering_new/lysozyme/dRI_df.rds")
# # con_df <- readRDS("/Users/fan/Documents/UCI/Research/light_scattering_new/lysozyme/con_df.rds")
# # ########################
# # ########################
# # pdf("lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_shorter_residuals_LS.pdf")
# # p_R <- ggplot(R_df[-1,], aes(x=factor(con_level), y=(prednew-observed)/observed)) + geom_boxplot(outlier.size = 0) + facet_wrap(~pH+NaCl,labeller = label_wrap_gen(multi_line=FALSE)) + geom_hline(yintercept=0, color = "red", size=1) + 
# #   xlab("Concentration level") + ylab("Relative difference in predicted and observed Rayleigh ratio") + ggtitle("") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=12),
# #                                                                                                                                                                                                                                                                                                         axis.text.x = element_text(angle = 90, hjust=1),
# #                                                                                                                                                                                                                                                                                                         axis.text=element_text(size=8))
# # 
# # print(p_R)
# # dev.off()
# # # ##########################
# # pdf("lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_shorter_residuals_dRI.pdf")
# # p_dRI <- ggplot(dRI_df[-1,], aes(x=factor(con_level), y=(prednew-observed)/observed)) + geom_boxplot(outlier.size = 0) + facet_wrap(~pH+NaCl, labeller = label_wrap_gen(multi_line=FALSE)) + geom_hline(yintercept=0, color = "red", size=1) +
# #   xlab("Concentration level") + ylab(expression(~"Relative difference in predicted and observed"~Delta~"n") )+ ggtitle("") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=12),
# #                                                                                                                                                                                                                                                                                                                 axis.text.x = element_text(angle = 90, hjust=1),
# #                                                                                                                                                                                                                                                                                                                 axis.text=element_text(size=8))
# # 
# # 
# # print(p_dRI)
# # dev.off()
# # #########################
# # pdf("lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_shorter_residuals_con.pdf")
# # p_con <- ggplot(con_df[-1,], aes(x=factor(con_level), y=(con/con.meas))) + geom_boxplot(outlier.size = 0) + facet_wrap(~pH+NaCl, labeller = label_wrap_gen(multi_line=FALSE)) + geom_hline(yintercept=1, color = "red", size=1) +
# #   xlab("Concentration level") + ylab("Inferred concentration / measured concentration")+ ggtitle("") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=12),
# #                                                                                                                                                                                                                                                                                           axis.text.x = element_text(angle = 90, hjust=1),
# #                                                                                                                                                                                                                                                                                           axis.text=element_text(size=8))
# # 
# # 
# # 
# # print(p_con)
# # dev.off()
# # ########################
# # #########################
# # # prior vs posterior
# # #########################
# # out <- readRDS("/Users/fan/Documents/UCI/Research/light_scattering_new/lysozyme/lysozyme_complete_pool_constant_dndc_noconmeas_log_normal_loose_prior_shorter_post_check.rds")
# # out_mcmc <- as.mcmc(out)
# # # A2
# # A2_prior_samples <- rnorm(draws/thin, mean = 0, sd = 1)
# # A2_post_samples <- c(out_mcmc[[1]][,c(paste0("A2[",1:14,"]"))])
# # A2_prior_post_df <- data.frame(A2 = c(A2_prior_samples, A2_post_samples),
# #                                type = rep(c("prior",c(paste0("A2[",1:14,"]"))), rep(400,15)) )
# # ggplot(A2_prior_post_df, aes(x=factor(type), y=A2)) + geom_boxplot(outlier.size = 0)
# # #########################
# # # dndc
# # dndc_prior_samples <- rnorm(draws/thin, mean = 0.197, sd = 0.005)
# # dndc_post_samples <- c(out_mcmc[[1]][,"dndc"])
# # dndc_prior_post_df <- data.frame(dndc = c(dndc_prior_samples, dndc_post_samples),
# #                                  type = rep(c("prior", "post"), c(400,400)) )
# # ggplot(dndc_prior_post_df, aes(x=factor(type), y=dndc)) + geom_boxplot(outlier.size = 0)
# # #########################
# # # s2R 
# # s2R_prior_samples <- invgamma::rinvgamma(n=draws/thin, shape=1, rate=(2/2)*(10^(-10)) )
# # s2R_post_samples <- c(out_mcmc[[1]][,"s2R"])
# # s2R_prior_post_df <- data.frame(s2R = c(s2R_prior_samples, s2R_post_samples),
# #                                 type = rep(c("prior", "post"), c(400,400)) )
# # ggplot(s2R_prior_post_df, aes(x=sqrt(s2R), fill=type)) + geom_density()
# # ggplot(s2R_prior_post_df, aes(x=sqrt(s2R), fill=type)) + geom_histogram()
# # ggplot(s2R_prior_post_df, aes(x=factor(type),y=sqrt(s2R), fill=type)) + geom_boxplot(outlier.size=0)
# # ###########################
# # # s2c
# # s2c_prior_samples <- invgamma::rinvgamma(n=draws/thin, shape=1, rate=(2/2)*(10^(-10)) )
# # s2c_post_samples <- c(out_mcmc[[1]][,"s2c"])
# # s2c_prior_post_df <- data.frame(s2R = c(s2R_prior_samples, s2R_post_samples),
# #                                 type = rep(c("prior", "post"), c(400,400)) )
# # ggplot(s2R_prior_post_df, aes(x=sqrt(s2R), fill=type)) + geom_density()