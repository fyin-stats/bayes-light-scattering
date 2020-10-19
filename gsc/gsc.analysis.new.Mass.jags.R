### figure out x-mer for gsc
#Analysis script for gamma-S crystallin light scattering data
#Meta-notes:
#  - All data is entered in reverse temporal order (most recent to oldest)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("foreach", "doParallel", 
              "boot", "parallel","rjags", "dplyr", "coda", 
              "ggplot2","R2jags", "R2WinBUGS")
ipak(packages)
#Load required libraries and such
mc.cores<-25
#Load the source file
source("./gsc/gsc.data.preprocess.R")
##############################
gsc_model <- function(){
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
  A2 ~ dnorm(0, 1/pow(prior.A2.sd,2))
  dndc ~ dnorm(prior.dndc.mean, 1/pow(prior.dndc.sd,2));T(0,)
  # N ~ dbern(0.5)
  Mass ~ dunif(monomer.mass, 20*monomer.mass)
  for(i in 1:n){
    for(l in 1:L){
      con[i,l] ~ dlnorm( log(con.meas[i,l]), inv.s2c )
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
        mu_R[i,j,l] <- f*con[i,l]*pow(n0[l]*dndc,2)*Mass*(1-2*A2*con[i,l]*Mass)
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
write.model(gsc_model, "./gsc/gsc_model_log_normal_Mass.bug")
# mle<-virialEstBoot(adat[[i]],solvent.base.RI=solvent.base.ri.meas[i],monomer.mass=monomer.mass,aggregate.size=1:3,
#                    RI.dat=adat[-i][fit.use[fit.use!=i]],
#                    method="MLE",use.RI.correction=FALSE, bootstrap.replicates=1,max.con.for.RI=RI.max.con)[[1]]$t0
# print(mle["M"])
jags.inits <- function(){
  list("dndc"=0.1985)
}
####
#Various MCMC settings
draws<-100000
thin<-250
burn<-100000
nadapt<-100000
nchains<-5
# read the preprocessed data
# gsc_dat_jags <- readRDS("./gsc/gsc_dat_jags.rds")
gsc_dat_jags <- list(scat.meas = datscat_array,
                     n0.meas=n0.meas,
                     n0=rep(mean(solvent.base.ri.meas[1:3]),3),
                     ri.meas=datri_matrix,
                     con.meas=con_matrix,
                     con.design=con_design_vector,
                     prior.A2.sd=1,
                     prior.dndc.sd=prior.dndc.sd,
                     prior.dndc.mean=prior.dndc.mean,
                     N=N_vec[i],
                     L=dim(datscat_array)[3],
                     nri=dim(datri_matrix)[1],
                     n=dim(datscat_array)[1],
                     npred = 100,
                     f=f,
                     monomer.mass=monomer.mass,
                     nang=dim(datscat_array)[2],
                     a_R = 2/2,
                     b_R = (2/2)*(10^(-10)),
                     a_n = 2/2,
                     b_n = (2/2)*(10^(-8)),
                     a_c = 2/2,
                     b_c = (2/2) * (log(1.05)/1.96)^2)
#
out <- do.call(jags.parallel, 
               list(data=gsc_dat_jags, inits=jags.inits, 
                    parameters.to.save = c("A2","dndc","s2R", "s2c", "s2n","con", "Mass",
                                           "bpval_dRI", "bpval_R", "residual_R.new", "residual_dRI.new",
                                           "scat.meas.new", "ri.meas.new"),
                    n.chains = nchains,
                    n.iter = burn+nadapt+draws,
                    n.burnin = burn+nadapt,
                    n.thin = thin, 
                    digits = 16, 
                    model.file="./gsc/gsc_model_log_normal_Mass.bug"))
########
time2 <- Sys.time()
####
saveRDS(out, "./gsc/gsc_log_normal_loose_Mass_post_check.rds")
out_mcmc <- as.mcmc(out)

####
# out <- readRDS("gsc_log_normal_loose_Mass_post_check.rds")
# out
# out <- readRDS("gsc_log_normal_loose_post_check.rds")
# MCMCtrace(out)

# gsc_log_normal_loose_Mass_post_check <- readRDS("gsc_log_normal_loose_Mass_post_check.rds")
# MCMCtrace(gsc_log_normal_loose_Mass_post_check)
# 
# 
# #
# out_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_post_check)
# data.frame(out_mcmc[[1]][,c("A2", "Mass")]) %>% ggplot(aes(x=A2, y=Mass)) + geom_density2d()