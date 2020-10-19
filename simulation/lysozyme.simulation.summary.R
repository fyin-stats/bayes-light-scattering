#######################################
###### lysozyme simulation studies summary
##### Last modified: 07/01/2020
#######################################
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
#######################################
#######################################
#######################################
simulation_summary <- readRDS("./simulation/simulation_summary_log_normal.rds")
simulation_summary_log_normal_no_adjustment <- readRDS("./simulation/simulation_summary_log_normal_no_adjustment.rds")
########################################
simulation_summary_all <- rbind(simulation_summary,
                                simulation_summary_log_normal_no_adjustment)
#  simulation_summary_10,
########################################
simulation_summary_all <- mutate(simulation_summary_all, con_error = exp(sc*1.96)-1)
##########
# setwd("/Users/fan/Documents/trunk/")
# point estimates
p_A2_001 <- simulation_summary_all %>% filter(abs(A2)==0.01) %>% 
  ggplot(aes(x=factor(num_experiments),y=(A2_mean_post_mean-A2)/A2, color=priors, shape=priors )) + 
  geom_point(size=3) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Relative Bias of A2") + xlab("Number of experiments") + geom_hline(yintercept=0,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                                                             axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                                                             strip.text = element_text(size = 10)) 
p_A2_001

pdf("./simulation/simulation_A2_10-2.pdf")
p_A2_001
dev.off()
###########
p_A2_0001 <- simulation_summary_all %>% filter(abs(A2)==0.001) %>% 
  ggplot(aes(x=factor(num_experiments),y=(A2_mean_post_mean-A2)/A2, color=priors, shape=priors )) + 
  geom_point(size=3) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Relative Bias of A2") + xlab("Number of experiments") + geom_hline(yintercept=0,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
p_A2_0001
pdf("./simulation/simulation_A2_10-3.pdf")
p_A2_0001
dev.off()
############
p_A2_00001 <- simulation_summary_all %>% filter(abs(A2)==0.0001) %>% 
  ggplot(aes(x=factor(num_experiments),y=(A2_mean_post_mean-A2)/A2, color=priors, shape=priors )) + 
  geom_point(size=2) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Relative Bias of A2") + xlab("Number of experiments") + geom_hline(yintercept=0,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
p_A2_00001
pdf("./simulation/simulation_A2_10-4.pdf")
p_A2_00001
dev.off()
#############
p_A2_000001 <- simulation_summary_all %>% filter(abs(A2)==0.00001) %>% 
  ggplot(aes(x=factor(num_experiments),y=(A2_mean_post_mean-A2)/A2, color=priors, shape=priors )) + 
  geom_point(size=2) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Relative Bias of A2") + xlab("Number of experiments") + geom_hline(yintercept=0,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
pdf("./simulation/simulation_A2_10-5.pdf")
p_A2_000001
dev.off()

#############
### coverage probability
#############
p_A2_001_coverage <- simulation_summary_all %>% filter(abs(A2) == 0.01) %>% 
  ggplot(aes(x=factor(num_experiments),y=A2_coverage_prob, color=priors, shape=priors )) + 
  geom_point(size=2) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Coverage probability of A2, 95% credible interval") + xlab("Number of experiments") + geom_hline(yintercept=0.95,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
p_A2_001_coverage

pdf("./simulation/simulation_A2_10-2_coverage.pdf")
p_A2_001_coverage
dev.off()
################
p_A2_0001_coverage <- simulation_summary_all %>% filter(abs(A2)==0.001) %>% 
  ggplot(aes(x=factor(num_experiments),y=A2_coverage_prob, color=priors, shape=priors )) + 
  geom_point(size=2) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Coverage probability of A2, 95% credible interval") + xlab("Number of experiments") + geom_hline(yintercept=0.95,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
p_A2_001_coverage

pdf("./simulation/simulation_A2_10-3_coverage.pdf")
p_A2_0001_coverage
dev.off()
##################
p_A2_00001_coverage <- simulation_summary_all %>% filter(abs(A2)==0.0001) %>% 
  ggplot(aes(x=factor(num_experiments),y=A2_coverage_prob, color=priors, shape=priors )) + 
  geom_point(size=2) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Coverage probability of A2, 95% credible interval") + xlab("Number of experiments") + geom_hline(yintercept=0.95,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
p_A2_00001_coverage

pdf("./simulation/simulation_A2_10-4_coverage.pdf")
p_A2_00001_coverage
dev.off()
####################
p_A2_000001_coverage <- simulation_summary_all %>% filter(abs(A2)==0.00001) %>% 
  ggplot(aes(x=factor(num_experiments),y=A2_coverage_prob, color=priors, shape=priors )) + 
  geom_point(size=2) + facet_wrap(~A2+con_error,nrow=2,label=label_both) + ylab("Coverage probability of A2, 95% credible interval") + xlab("Number of experiments") + geom_hline(yintercept=0.95,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=10),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 10)) 
p_A2_000001_coverage

pdf("./simulation/simulation_A2_10-5_coverage.pdf")
p_A2_000001_coverage
dev.off()
#############
# dndc
#############
pdf("./simulation/dndc_relative_bias.pdf")
simulation_summary_all %>% 
  ggplot(aes(x=factor(num_experiments),y=(dndc_mean_post_mean-dndc)/dndc, color=priors, shape=factor(con_error))) + 
  geom_point(size=3) +  facet_wrap(~A2,label=label_both) + ylab("Relative Bias of dndc") + xlab("Number of experiments") + geom_hline(yintercept=0,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=15),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 15)) 
dev.off()



pdf("./simulation/dndc_coverage.pdf")
simulation_summary_all %>% 
  ggplot(aes(x=factor(num_experiments),y=dndc_coverage_prob, color=priors, shape=factor(con_error))) + 
  geom_point(size=3) + facet_wrap(~A2,label=label_both) + ylab("Coverage probability of dndc, 95% credible interval") + xlab("Number of experiments") + geom_hline(yintercept=0.95,colour="red") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_rect(colour = "black")) + theme(text = element_text(size=15),
                                                                                                                                                                                     axis.text.x = element_text(angle = 0, hjust=1),
                                                                                                                                                                                     strip.text = element_text(size = 15)) 
dev.off()


##############
# Table
##############
# A2
A2_summary_table <- mutate(simulation_summary_all,
                           Absolute_bias = A2_mean_post_mean-A2,
                           Relative_bias = (A2_mean_post_mean-A2)/A2,
                           Coverage_prob = A2_coverage_prob,
                           Avg_width = A2_width,
                           Coverage_Width = A2_coverage_prob/A2_width)







##############
# Table for A2
simulation_summary_table_relative_bias <- matrix(NA, nrow = 1, ncol = 3+16)
# 4 parts, priors
unique_A2 <- unique(simulation_summary_all$A2)
unique_con_error <- unique(simulation_summary_all$con_error)
unique_prior <- unique(simulation_summary_all$priors)
unique_num_experiments <- unique(simulation_summary_all$num_experiments)


# relative bias
for(i in 1:length(unique_A2)){
  temp_table <- matrix(NA, nrow = length(unique_con_error), ncol = 3+16)
  temp_table[,1] <- abs(unique_A2[i])
  temp_table[,2] <- unique_A2[i]
  for(j in 1:length(unique_con_error)){
    # first column absolute value of A2
    # simulation_summary_table_relative_bias[,1] <- 
    temp_table[j,3] <- unique_con_error[j]
    # informative
    temp_table[j,4] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,5] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,6] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,7] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    # somewhat informative
    temp_table[j,8] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,9] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,10] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,11] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    # weakly informative
    temp_table[j,12] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,13] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,14] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,15] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    # no adjustment
    temp_table[j,16] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,17] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,18] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
    temp_table[j,19] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_bias = (A2_mean_post_mean-A2)/A2))$relative_bias
  }
  simulation_summary_table_relative_bias <- rbind(simulation_summary_table_relative_bias, temp_table)
}
#
simulation_summary_table_relative_bias <- simulation_summary_table_relative_bias[complete.cases(simulation_summary_table_relative_bias),]
simulation_summary_table_relative_bias <- data.frame(simulation_summary_table_relative_bias)
colnames(simulation_summary_table_relative_bias) <- c("abs_A2",
                                                      "A2",
                                                      "true_errors",
                                                      rep(c(1,2,5,10),4))

# sort by absolute value of A2, from large to small
simulation_summary_table_relative_bias <- simulation_summary_table_relative_bias %>% data.frame() %>% arrange(desc(abs_A2))
# convert to percentage
simulation_summary_table_relative_bias[,3] <- 100*simulation_summary_table_relative_bias[,3]
  # paste0(100 * c(matrix(simulation_summary_table_relative_bias[,4:19])), "%")
xtable(simulation_summary_table_relative_bias)


#### frequentist coverage 
simulation_summary_table_coverage <- matrix(NA, nrow = 1, ncol = 3+16)
####
for(i in 1:length(unique_A2)){
  temp_table <- matrix(NA, nrow = length(unique_con_error), ncol = 3+16)
  temp_table[,1] <- abs(unique_A2[i])
  temp_table[,2] <- unique_A2[i]
  for(j in 1:length(unique_con_error)){
    # first column absolute value of A2
    # simulation_summary_table_relative_bias[,1] <- 
    temp_table[j,3] <- unique_con_error[j]
    # informative
    temp_table[j,4] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob 
    temp_table[j,5] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,6] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,7] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    # somewhat informative
    temp_table[j,8] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,9] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,10] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,11] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    # weakly informative
    temp_table[j,12] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,13] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,14] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,15] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    # no adjustment
    temp_table[j,16] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,17] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,18] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
    temp_table[j,19] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]))$A2_coverage_prob
  }
  simulation_summary_table_coverage <- rbind(simulation_summary_table_coverage, temp_table)
}
######
simulation_summary_table_coverage <- simulation_summary_table_coverage[complete.cases(simulation_summary_table_coverage),]
simulation_summary_table_coverage <- data.frame(simulation_summary_table_coverage)
colnames(simulation_summary_table_coverage) <- c("abs_A2",
                                                      "A2",
                                                      "true_errors",
                                                      rep(c(1,2,5,10),4))

# sort by absolute value of A2, from large to small
simulation_summary_table_coverage <- simulation_summary_table_coverage %>% data.frame() %>% arrange(desc(abs_A2))

simulation_summary_table_coverage[,3] <- 100 * simulation_summary_table_coverage[,3]
xtable(simulation_summary_table_coverage)




##############################
### average width / true value
simulation_summary_table_width <- matrix(NA, nrow = 1, ncol = 3+16)
####
for(i in 1:length(unique_A2)){
  temp_table <- matrix(NA, nrow = length(unique_con_error), ncol = 3+16)
  temp_table[,1] <- abs(unique_A2[i])
  temp_table[,2] <- unique_A2[i]
  for(j in 1:length(unique_con_error)){
    # first column absolute value of A2
    # simulation_summary_table_relative_bias[,1] <- 
    temp_table[j,3] <- unique_con_error[j]
    # informative
    temp_table[j,4] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,5] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,6] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,7] <- (filter(simulation_summary_all, priors=="informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    # somewhat informative
    temp_table[j,8] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,9] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,10] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,11] <- (filter(simulation_summary_all, priors=="somewhat informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    # weakly informative
    temp_table[j,12] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,13] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,14] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,15] <- (filter(simulation_summary_all, priors=="weakly informative" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    # no adjustment
    temp_table[j,16] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 1 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,17] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 2 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,18] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 5 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
    temp_table[j,19] <- (filter(simulation_summary_all, priors=="no adjustment" & num_experiments == 10 & con_error == unique_con_error[j] & A2 == unique_A2[i]) %>% mutate(relative_width = (A2_width)/A2))$relative_width %>% abs()
  }
  simulation_summary_table_width <- rbind(simulation_summary_table_width, temp_table)
}

######
simulation_summary_table_width <- simulation_summary_table_width[complete.cases(simulation_summary_table_width),]
simulation_summary_table_width <- data.frame(simulation_summary_table_width)
colnames(simulation_summary_table_width) <- c("abs_A2",
                                                 "A2",
                                                 "true_errors",
                                                 rep(c(1,2,5,10),4))

# sort by absolute value of A2, from large to small
simulation_summary_table_width <- simulation_summary_table_width %>% data.frame() %>% arrange(desc(abs_A2))

simulation_summary_table_width[,3] <- 100 * simulation_summary_table_width[,3]
xtable(simulation_summary_table_width)



