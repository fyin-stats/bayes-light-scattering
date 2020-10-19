################ gsc post summary
################
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("foreach", "doParallel", 
              "boot", "parallel","rjags", 
              "ggplot2", "dplyr", "invgamma", "R2jags", "MCMCvis", 
              "R2WinBUGS",
              "ggmcmc", "ggpubr", "xtable", "scales","AICcmodavg")
ipak(packages)
###################
monomer.mass <- 20959.80
### xmer: DIC
gsc_xmer <- readRDS("./gsc/gsc_log_normal_loose_xmer_post_check.rds")
### Model: 21, 22
gsc_log_normal_loose_Mass_post_check <- readRDS("./gsc/gsc_log_normal_loose_Mass_post_check.rds")
gsc_log_normal_loose_Mass_normal_post_check <- readRDS("./gsc/gsc_log_normal_loose_Mass_normal_post_check.rds")###
###
gsc_log_normal_loose_Mass_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_post_check)
gsc_log_normal_loose_Mass_normal_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_normal_post_check)
###
xmer_DIC <- data.frame(x = 1:20,
                       DIC = do.call("c", lapply(gsc_xmer, function(x) AICcmodavg::DIC(x))))
pdf("./gsc/DICxmer.pdf")
par(mgp=c(2,1,0)) 
plot(DIC~x, data=xmer_DIC, type="b",cex.lab=2, cex.axis = 1.5, lwd=3)
points(x=mean(do.call("c", lapply(gsc_log_normal_loose_Mass_post_check_mcmc, function(x) x[,c("Mass")])))/monomer.mass,
       y=AICcmodavg::DIC(gsc_log_normal_loose_Mass_post_check), col="blue", pch=1, cex=3)
points(x=mean(do.call("c", lapply(gsc_log_normal_loose_Mass_normal_post_check_mcmc, function(x) x[,c("Mass")])))/monomer.mass,
       y=AICcmodavg::DIC(gsc_log_normal_loose_Mass_normal_post_check), col="red", pch=2, cex=3)
# abline(h=AICcmodavg::DIC(gsc_log_normal_loose_Mass_post_check), col="blue", lty=2, lwd=3) # uniform mass prior
# abline(h=AICcmodavg::DIC(gsc_log_normal_loose_Mass_normal_post_check), col="red", lwd=3) # normal mass prior
# legend("topright", legend=c("M21", "M22"),
#        col=c("blue","red"), lty=1:2, cex=2, lwd=3)
legend("topright", legend=c("M21", "M22"),
       col=c("blue","red"), pch=1:2, cex=2)
dev.off()
#####################
#####################
gsc_12mer <- gsc_xmer[[12]]
# readRDS("./gsc/gsc_12mer.rds")
######################
# check convergence
# MCMCtrace(gsc_log_normal_loose_Mass_normal_post_check) # converged
# MCMCtrace(gsc_log_normal_loose_Mass_post_check)
# MCMCtrace(gsc_12mer,n.eff = T)
######################
gsc_log_normal_loose_Mass_normal_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_normal_post_check)
gsc_log_normal_loose_Mass_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_post_check)
gsc_12mer_mcmc <- as.mcmc(gsc_12mer)
######################
######################
# dndc
gsc_dndc_df <- data.frame(dndc = c(gsc_log_normal_loose_Mass_post_check_mcmc[[1]][,"dndc"],
                               gsc_log_normal_loose_Mass_normal_post_check_mcmc[[1]][,"dndc"],
                               gsc_12mer_mcmc[[1]][,"dndc"]),
                        model = rep(c("M21", "M22", "M12"), rep(400,3)))
p_gsc_dndc_density <- gsc_dndc_df %>% ggplot(aes(x=model,y=dndc)) + geom_boxplot() + xlab("") + ylab("") + ggtitle("dndc") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                                                   axis.text.x = element_text(angle = 0, hjust=1, size=20),
                                                                                                                                                                                                                                                                                                                                                                                   strip.text = element_text(size = 20)) 

#######################
######################
# density of Mw
gsc_Mw_df <- data.frame(Mw = c(gsc_log_normal_loose_Mass_post_check_mcmc[[1]][,"Mass"],
                               gsc_log_normal_loose_Mass_normal_post_check_mcmc[[1]][,"Mass"]),
                        model = rep(c("M21", "M22"), rep(400,2)))
p_gsc_Mw_density <- gsc_Mw_df %>% ggplot(aes(x=model,y=Mw)) + geom_boxplot() + geom_hline(yintercept=monomer.mass*12,colour="red",size=1.5) + xlab("") + ylab("") + ggtitle("Mw") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                                                           axis.text.x = element_text(angle = 0, hjust=1, size=20),
                                                                                                                                                                                                                                                                                           strip.text = element_text(size = 20)) 

#######################
# posterior density of A2
########################
gsc_A2_df <- data.frame(A2 = c(gsc_log_normal_loose_Mass_post_check_mcmc[[1]][,"A2"],
                               gsc_log_normal_loose_Mass_normal_post_check_mcmc[[1]][,"A2"],
                               gsc_12mer_mcmc[[1]][,"A2"]),
                        model = rep(c("M21", "M22", "M12"), rep(400,3)))
gsc_A2_df %>% group_by(model) %>% summarise(A2_median = median(A2)*10^5,
                                            p_positive = mean(A2>0))

p_gsc_A2_density <- gsc_A2_df %>% ggplot(aes(x=model,y=A2)) + geom_boxplot() + xlab("") + ylab("") + ggtitle(expression(A[~"2"])) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                                                                                  axis.text.x = element_text(angle = 0, hjust=1, size=20),
                                                                                                                                                                                                                                                                                                                strip.text = element_text(size = 20)) 

#######################
# s2u
gsc_s2u_df <- data.frame(s2u = c(gsc_log_normal_loose_Mass_post_check_mcmc[[1]][,"s2c"],
                                 gsc_log_normal_loose_Mass_normal_post_check_mcmc[[1]][,"s2c"],
                                 gsc_12mer_mcmc[[1]][,"s2c"]),
                          model = rep(c("M21", "M22", "M12"), rep(400,3)))
p_gsc_su_density <- gsc_s2u_df %>% ggplot(aes(x=model,y=sqrt(s2u))) + geom_boxplot() + xlab("") + ylab("") + ggtitle(expression(sigma[~"u"])) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=20),
                                                                                                                                                                                                                                                                                                                      axis.text.x = element_text(angle = 0, hjust=1, size=20),
                                                                                                                                                                                                                                                                                                                      strip.text = element_text(size = 20)) 

#######################

#######################


figure <- ggarrange(p_gsc_A2_density,p_gsc_Mw_density,
                    p_gsc_dndc_density, p_gsc_su_density,
                    ncol = 2, nrow = 2)
figure
#######################
pdf("./gsc/gsc_A2_Mw.pdf")
figure
dev.off()


#########################
# posterior predictive checks
##########################
gsc_log_normal_loose_Mass_post_check <- readRDS("./gsc/gsc_log_normal_loose_Mass_post_check.rds")
gsc_log_normal_loose_Mass_normal_post_check <- readRDS("./gsc/gsc_log_normal_loose_Mass_normal_post_check.rds")
gsc_log_normal_loose_12mer_post_check <- gsc_xmer[[12]]
gsc_dat_jags <- readRDS("./gsc/gsc_dat_jags.rds")
dat_jags <- gsc_dat_jags
###########################
gsc_log_normal_loose_Mass_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_post_check)
gsc_log_normal_loose_Mass_normal_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_Mass_normal_post_check)
gsc_log_normal_loose_12mer_post_check_mcmc <- as.mcmc(gsc_log_normal_loose_12mer_post_check)
############################
############################
out_mcmc <- gsc_log_normal_loose_Mass_post_check_mcmc
############################
R_df_Mass <- data.frame(residuals = NA,
                   prednew = NA,
                   observed = NA,
                   con_level = NA,
                   pH = NA,
                   NaCl = NA,
                   replicate=NA)
dRI_df_Mass <- data.frame(residuals=NA,
                     prednew = NA,
                     observed = NA, 
                     con_level = NA,
                     pH = NA,
                     NaCl = NA,
                     replicate=NA)
con_df_Mass <- data.frame(con = NA, 
                     con.meas = NA,
                     con_level = NA,
                     pH = NA,
                     NaCl = NA,
                     replicate=l)
#################
L = dim(dat_jags$scat.meas)[3]
IR = dim(dat_jags$scat.meas)[1]
# use.experimental.conditions <- experimental.conditions[1:17,][-c(11,12,14),]
##################
for(l in 1:L){
  for(i in 1:IR){
    R_residuals_colname <- paste0("residual_R.new","[",i,",","1",",",l,"]")
    R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
    # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    R_df_Mass <- rbind(R_df_Mass, data.frame(residuals = c(out_mcmc[[1]][,R_residuals_colname]),
                                   prednew = c(out_mcmc[[1]][,R_pred_colname]),
                                   observed = dat_jags$scat.meas[i,1,l],
                                   con_level = i,
                                   pH = 6.9,
                                   NaCl = 100,
                                   replicate=l))
  }
}
####################
for(l in 1:L){
  for(i in 1:IR){
    con_colname <- paste0("con","[",i,",",l,"]")
    # R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
    # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    con_df_Mass <- rbind(con_df_Mass, data.frame(con = c(out_mcmc[[1]][,con_colname]),
                                       con.meas = dat_jags$con.meas[i,l],
                                       con_level = i,
                                       pH = 6.9,
                                       NaCl = 100,
                                       replicate = l))
  }
}
####################
IdRI <- dim(dat_jags$ri.meas)[1]
for(l in 1:L){
  for(i in 1:IdRI){
    dRI_residuals_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    dRI_pred_colname <- paste0("ri.meas.new","[",i,",",l,"]")
    dRI_df_Mass <- rbind(dRI_df_Mass, data.frame(residuals = c(out_mcmc[[1]][,dRI_residuals_colname]),
                                       prednew = c(out_mcmc[[1]][,dRI_pred_colname]),
                                       observed = dat_jags$ri.meas[i,l],
                                       con_level = i,
                                       pH = 6.9,
                                       NaCl = 100,
                                       replicate=l))
  }
}
##################################################
######
out_mcmc <- gsc_log_normal_loose_Mass_normal_post_check_mcmc
############################
R_df_Mass_normal <- data.frame(residuals = NA,
                        prednew = NA,
                        observed = NA,
                        con_level = NA,
                        pH = NA,
                        NaCl = NA,
                        replicate=NA)
dRI_df_Mass_normal <- data.frame(residuals=NA,
                          prednew = NA,
                          observed = NA, 
                          con_level = NA,
                          pH = NA,
                          NaCl = NA,
                          replicate=NA)
con_df_Mass_normal <- data.frame(con = NA, 
                          con.meas = NA,
                          con_level = NA,
                          pH = NA,
                          NaCl = NA,
                          replicate=l)
#################
L = dim(dat_jags$scat.meas)[3]
IR = dim(dat_jags$scat.meas)[1]
# use.experimental.conditions <- experimental.conditions[1:17,][-c(11,12,14),]
##################
for(l in 1:L){
  for(i in 1:IR){
    R_residuals_colname <- paste0("residual_R.new","[",i,",","1",",",l,"]")
    R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
    # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    R_df_Mass_normal <- rbind(R_df_Mass_normal, data.frame(residuals = c(out_mcmc[[1]][,R_residuals_colname]),
                                             prednew = c(out_mcmc[[1]][,R_pred_colname]),
                                             observed = dat_jags$scat.meas[i,1,l],
                                             con_level = i,
                                             pH = 6.9,
                                             NaCl = 100,
                                             replicate=l))
  }
}
####################
for(l in 1:L){
  for(i in 1:IR){
    con_colname <- paste0("con","[",i,",",l,"]")
    # R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
    # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    con_df_Mass_normal <- rbind(con_df_Mass_normal, data.frame(con = c(out_mcmc[[1]][,con_colname]),
                                                 con.meas = dat_jags$con.meas[i,l],
                                                 con_level = i,
                                                 pH = 6.9,
                                                 NaCl = 100,
                                                 replicate = l))
  }
}
####################
IdRI <- dim(dat_jags$ri.meas)[1]
for(l in 1:L){
  for(i in 1:IdRI){
    dRI_residuals_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    dRI_pred_colname <- paste0("ri.meas.new","[",i,",",l,"]")
    dRI_df_Mass_normal <- rbind(dRI_df_Mass_normal, data.frame(residuals = c(out_mcmc[[1]][,dRI_residuals_colname]),
                                                 prednew = c(out_mcmc[[1]][,dRI_pred_colname]),
                                                 observed = dat_jags$ri.meas[i,l],
                                                 con_level = i,
                                                 pH = 6.9,
                                                 NaCl = 100,
                                                 replicate=l))
  }
}

#######################
out_mcmc <- gsc_log_normal_loose_12mer_post_check_mcmc
# gsc_log_normal_loose_12mer_post_check_mcmc
############################
R_df_12 <- data.frame(residuals = NA,
                        prednew = NA,
                        observed = NA,
                        con_level = NA,
                        pH = NA,
                        NaCl = NA,
                        replicate=NA)
dRI_df_12 <- data.frame(residuals=NA,
                          prednew = NA,
                          observed = NA, 
                          con_level = NA,
                          pH = NA,
                          NaCl = NA,
                          replicate=NA)
con_df_12 <- data.frame(con = NA, 
                          con.meas = NA,
                          con_level = NA,
                          pH = NA,
                          NaCl = NA,
                          replicate=NA)
#################
L = dim(dat_jags$scat.meas)[3]
IR = dim(dat_jags$scat.meas)[1]
# use.experimental.conditions <- experimental.conditions[1:17,][-c(11,12,14),]
##################
for(l in 1:L){
  for(i in 1:IR){
    R_residuals_colname <- paste0("residual_R.new","[",i,",","1",",",l,"]")
    R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
    # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    R_df_12 <- rbind(R_df_12, data.frame(residuals = c(out_mcmc[[1]][,R_residuals_colname]),
                                             prednew = c(out_mcmc[[1]][,R_pred_colname]),
                                             observed = dat_jags$scat.meas[i,1,l],
                                             con_level = i,
                                             pH = 6.9,
                                             NaCl = 100,
                                             replicate=l))
  }
}
####################
for(l in 1:L){
  for(i in 1:IR){
    con_colname <- paste0("con","[",i,",",l,"]")
    # R_pred_colname <- paste0("scat.meas.new","[",i,",","1",",",l,"]")
    # dRI_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    con_df_12 <- rbind(con_df_12, data.frame(con = c(out_mcmc[[1]][,con_colname]),
                                                 con.meas = dat_jags$con.meas[i,l],
                                                 con_level = i,
                                                 pH = 6.9,
                                                 NaCl = 100,
                                                 replicate = l))
  }
}
####################
IdRI <- dim(dat_jags$ri.meas)[1]
for(l in 1:L){
  for(i in 1:IdRI){
    dRI_residuals_colname <- paste0("residual_dRI.new","[",i,",",l,"]")
    dRI_pred_colname <- paste0("ri.meas.new","[",i,",",l,"]")
    dRI_df_12 <- rbind(dRI_df_12, data.frame(residuals = c(out_mcmc[[1]][,dRI_residuals_colname]),
                                                 prednew = c(out_mcmc[[1]][,dRI_pred_colname]),
                                                 observed = dat_jags$ri.meas[i,l],
                                                 con_level = i,
                                                 pH = 6.9,
                                                 NaCl = 100,
                                                 replicate=l))
  }
}
#####################
# Mass
R_df_Mass$model = "M21"
dRI_df_Mass$model = "M21"
con_df_Mass$model = "M21"
# Mass normal
R_df_Mass_normal$model = "M22"
dRI_df_Mass_normal$model = "M22"
con_df_Mass_normal$model = "M22"
# M12
R_df_12$model = "M12"
dRI_df_12$model = "M12"
con_df_Mass_normal$model = "M12"
#######################
R_post_pred_Mass <- R_df_Mass[-1,] %>% group_by(con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
R_post_pred_Mass_normal <- R_df_Mass_normal[-1,] %>% group_by(con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
R_post_pred_12 <- R_df_12[-1,] %>% group_by(con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
xtable(cbind(R_post_pred_12[,2], R_post_pred_Mass[,2], R_post_pred_Mass_normal[,2]), digits=3)
# dRI_df_Mass[-1,] %>% group_by(replicate, con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
dRI_post_pred_Mass <- dRI_df_Mass[-1,] %>% group_by(con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
dRI_post_pred_Mass_normal <- dRI_df_Mass_normal[-1,] %>% group_by(con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
dRI_post_pred_12 <- dRI_df_12[-1,] %>% group_by(con_level) %>% summarise(post_pred_pval = mean(prednew>observed))
xtable(cbind(dRI_post_pred_12[,2], dRI_post_pred_Mass[,2], dRI_post_pred_Mass_normal[,2]), digits=3)
