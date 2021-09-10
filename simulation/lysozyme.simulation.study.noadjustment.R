######################################
######################################
### lysozyme simulation study
######################################
######################################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("scales",
              "foreach", "doParallel", 
              "dplyr", "xtable",
              "coda", "mclust", "mcclust", 
              "ggplot2", "sna", "intergraph", "truncnorm",
              "boot", "parallel","rjags", "invgamma", "R2jags",
              "R2WinBUGS")
ipak(packages)
detectCores()
#######################################
registerDoParallel(cores = 34)
#######################################
LS_data_simulator <- function(A2, sn=sqrt(2*10^(-9)), sR=sqrt(1*10^(-11)), dndc=0.20, sc, nang = 1,
                              con.levels=c(2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
                              num_experiments=1, wl = 658e-7, n0=1.33, monomer.mass=14307, N=0, seed=12345){
  ### generate concentration data
  # total concentration levels
  set.seed(seed)
  con.meas = matrix(rep(con.levels, rep(num_experiments, length(con.levels))),ncol=1)
  I = length(con.meas)
  L = 1 # only one experimental condition
  # actual concentration
  con <- matrix(0, nrow=I, ncol=L)
  for(i in 1:I){
    con[i,1] <- rlnorm(1, meanlog = log(con.meas[i,1]), sdlog = sc)
  }
  # dndc measurements
  nri <- sum(con.meas<=0.02)
  ri.meas <- matrix(0,nrow=nri,ncol=L)
  for(i in 1:nri){
    ri.meas[i,1] <- rtruncnorm(1, a=0, b=Inf, mean = dndc*con[i,1], sd=sn)
  }
  # LS measurements
  # dn.meas <- rtruncnorm()
  # simulate Light scattering measurement
  f <- (4*pi^2)/(wl^4*6.022e23)
  K <- 4*pi^2*n0^2*dndc^2 / (6.022*10^23 * wl^4)
  # scat.meas <- matrix(NA, nrow = n, ncol=nang)
  scat.meas <- array(0, dim=c(I,1,L))
  for(i in 1:I){
    scat.meas[i,1,1] <- rnorm(1,
                              mean=f*con[i,1]*(n0*dndc)^2*(N+1)*monomer.mass*(1-2*A2*con[i,1]*(N+1)*monomer.mass), 
                              sd=sR)
  }
  # f*con[i,l]*pow(n0[l]*dndc,2)*(N+1)*monomer.mass*(1-2*A2[l]*con[i,l]*(N+1)*monomer.mass)
  return(list(scat.meas = scat.meas, ri.meas = ri.meas, con.meas = con.meas, n = I, nri = nri, nang=1, L=L,
              monomer.mass=monomer.mass, f=(4*pi^2)/(wl^4*6.022e23) ) )
}
##############################################
##############################################
#### simulation settings
A2_vec <- c(10^(-2), 10^(-3), 10^(-4), 10^(-5), -10^(-2), -10^(-3), -10^(-4), -10^(-5)) 
num_experiments_vec <- c(1,2,5,10)
# num_experiments_vec <- c(10)
sc_vec <- c(log(1.01)/1.96, log(1.05)/1.96, log(1.10)/1.96, log(1.20)/1.96)
# priors_vec <- c("informative", "somewhat informative", "weakly informative")
priors_vec <- c("no adjustment")
num_replicates <- 100
### MCMC settings
draws<-300000
thin<- 250
burn<-100000
nadapt<-100000
# three priors: inverse-chi squared (), inverse gamma (shape=1, but with rate = true value), inverse gamma (shape=1, but with mean = log(1.25)/1.96)
##############################################
#### how do you store the results?
##############################################
simulation_settings_df <- expand.grid(A2=A2_vec, num_experiments = num_experiments_vec, 
                                      sc = sc_vec, priors = priors_vec)
##############################################
# model specification
lysozyme_model <- function(){
  #mu_dndc ~ dnorm(prior.dndc.mu.mean, 1/pow(prior.dndc.mu.sd,2))
  #inv.s2dndc ~ dgamma(a_prior.dndc.s2, b_prior.dndc.s2)
  #s2dndc <- 1/inv.s2dndc
  ## variance parameters (gamma in jags, shape rate)
  inv.s2R ~ dgamma(a_R, b_R)
  inv.s2n ~ dgamma(a_n, b_n)
  s2R <- 1/inv.s2R
  s2n <- 1/inv.s2n
  
  inv.s2c ~ dgamma(a_c, b_c)
  s2c <- 1/inv.s2c
  ### 
  for(l in 1:L){
    A2[l] ~ dnorm(0, 1/pow(prior.A2.sd,2) )
  }
  
  dndc ~ dnorm(0.195, 1/pow(0.005,2));T(0,)
  
  for(i in 1:n){
    for(l in 1:L){
      # con[i,l] ~ dlnorm( log(con.meas[i,l]), inv.s2c )
      con[i,l] <- con.meas[i,l]
    }
  }
  
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
}
#
write.model(lysozyme_model, "./simulation/lysozyme_sim_model_no_adjustment.bug")
##
simulation_settings_df$dndc <- 0.20
simulation_settings_df$s2R <- 1*10^(-11)
simulation_settings_df$s2n <- 2*10^(-9)
simulation_settings_df$s2c <- simulation_settings_df$sc^2
simulation_summary <- simulation_settings_df
###
simulation_summary$A2_mean_post_mean <- NA
simulation_summary$dndc_mean_post_mean <- NA
simulation_summary$s2R_mean_post_mean <- NA
simulation_summary$s2n_mean_post_mean <- NA
simulation_summary$s2c_mean_post_mean <- NA
# coverage prob
simulation_summary$A2_coverage_prob <- NA
simulation_summary$dndc_coverage_prob <- NA
simulation_summary$s2R_coverage_prob <- NA
simulation_summary$s2n_coverage_prob <- NA
simulation_summary$s2c_coverage_prob <- NA
# average width of the 95% interval
simulation_summary$A2_width <- NA
simulation_summary$dndc_width <- NA
simulation_summary$s2R_width <- NA
simulation_summary$s2n_width <- NA
simulation_summary$s2c_width <- NA
# time
simulation_summary$time <- NA
##############################################
for(i in 1:nrow(simulation_settings_df)){
  # figure out the setting
  # if(simulation_settings_df$priors[i] == "informative"){
  #   b_c = 0.5
  #   a_c = 1/(2*simulation_settings_df$sc[i]^2) + 1
  # }
  # else if(simulation_settings_df$priors[i] == "somewhat informative"){
  #   a_c = 1
  #   b_c = simulation_settings_df$sc[i]^2
  # }
  # else {
  #   a_c = 1
  #   b_c = log(1.25)/1.96
  # }
  #
  a_c = 1
  b_c = 1
  cat("This is the", i, "-th iteration")
  cat("\n")
  time1 <- Sys.time()
  rlt_list <- foreach(j = 1:num_replicates)%dopar%{
    # generate the data
    dat = LS_data_simulator(A2 = simulation_settings_df$A2[i],
                            sc = simulation_settings_df$sc[i],
                            num_experiments = simulation_settings_df$num_experiments[i], 
                            seed=12345*j)
    dat_jags <- list(scat.meas = dat$scat.meas,
                     n0=1.33,
                     ri.meas=dat$ri.meas,
                     con.meas=dat$con.meas,
                     prior.A2.sd=1,
                     N=0,
                     L=dat$L,
                     nri=dat$nri,
                     n=dat$n,
                     f=dat$f,
                     monomer.mass=dat$monomer.mass,
                     nang=1,
                     a_R = 2/2,
                     b_R = (2/2)*(10^(-10)),
                     a_n = 2/2,
                     b_n = (2/2)*(10^(-8)),
                     a_c = a_c,
                     b_c = b_c )
    # estimation
    jags.inits <- function(){
      list("dndc"=0.197)
    }
    ####
    out <-  jags(data=dat_jags, inits=jags.inits, 
                 parameters.to.save = c("A2","s2R", "s2c", "s2n", "dndc"),
                 n.chains = 1,
                 n.iter = burn+nadapt+draws,
                 n.burnin = burn+nadapt,
                 n.thin = thin, 
                 digits = 16, 
                 model.file="./simulation/lysozyme_sim_model_no_adjustment.bug")
    out_mcmc <- as.mcmc(out)
    out_mcmc
  }
  time2 <- Sys.time()
  cat("The", i, "-th iteration is done")
  cat("\n")
  # generate the summary statistics
  # mean of posterior means
  #browser()
  simulation_summary$A2_mean_post_mean[i] <- mean(do.call("c", lapply(rlt_list, function(x) mean(x[[1]][,"A2"]))))
  simulation_summary$dndc_mean_post_mean[i] <- mean(do.call("c", lapply(rlt_list, function(x) mean(x[[1]][,"dndc"]))))
  simulation_summary$s2R_mean_post_mean[i] <- mean(do.call("c", lapply(rlt_list, function(x) mean(x[[1]][,"s2R"]))))
  simulation_summary$s2n_mean_post_mean[i] <- mean(do.call("c", lapply(rlt_list, function(x) mean(x[[1]][,"s2n"]))))
  simulation_summary$s2c_mean_post_mean[i] <- mean(do.call("c", lapply(rlt_list, function(x) mean(x[[1]][,"s2c"]))))
  # coverage prob
  simulation_summary$A2_coverage_prob[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    (quantile(x[[1]][,"A2"],0.025) <= simulation_settings_df$A2[i]) * (quantile(x[[1]][,"A2"],0.975) >= simulation_settings_df$A2[i])   )  )) 
  simulation_summary$dndc_coverage_prob[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    (quantile(x[[1]][,"dndc"],0.025) <= simulation_settings_df$dndc[i]) * (quantile(x[[1]][,"dndc"],0.975) >= simulation_settings_df$dndc[i])   )  )) 
  
  simulation_summary$s2R_coverage_prob[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    (quantile(x[[1]][,"s2R"],0.025) <= simulation_settings_df$s2R[i]) * (quantile(x[[1]][,"s2R"],0.975) >= simulation_settings_df$s2R[i])   )  )) 
  
  simulation_summary$s2n_coverage_prob[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    (quantile(x[[1]][,"s2n"],0.025) <= simulation_settings_df$s2n[i]) * (quantile(x[[1]][,"s2n"],0.975) >= simulation_settings_df$s2n[i])   )  )) 
  
  simulation_summary$s2c_coverage_prob[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    (quantile(x[[1]][,"s2c"],0.025) <= simulation_settings_df$s2c[i]) * (quantile(x[[1]][,"s2c"],0.975) >= simulation_settings_df$s2c[i])   )  )) 
  # average width of the 95% interval
  simulation_summary$A2_width[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    abs(quantile(x[[1]][,"A2"],0.025) - quantile(x[[1]][,"A2"],0.975))   )))
  simulation_summary$dndc_width[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    abs(quantile(x[[1]][,"dndc"],0.025) - quantile(x[[1]][,"dndc"],0.975))   )))
  simulation_summary$s2R_width[i] <-mean(do.call("c", lapply(rlt_list, function(x) 
    abs(quantile(x[[1]][,"s2R"],0.025) - quantile(x[[1]][,"s2R"],0.975))   )))
  simulation_summary$s2n_width[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    abs(quantile(x[[1]][,"s2n"],0.025) - quantile(x[[1]][,"s2n"],0.975))   )))
  simulation_summary$s2c_width[i] <- mean(do.call("c", lapply(rlt_list, function(x) 
    abs(quantile(x[[1]][,"s2c"],0.025) - quantile(x[[1]][,"s2c"],0.975))   )))
  # time
  simulation_summary$time[i] <- time2-time1
  sink("./simulation/lysozyme_simulation_summary.txt", append = TRUE)
  print(simulation_summary[i,])
  sink()
}
#############################
saveRDS(simulation_summary, 
        "./simulation/simulation_summary_log_normal_no_adjustment.rds")
