# Readme file for the Bayes-light-scattering codes
# Please do not forget to set up the path as the local path of this folder before running the following scripts
# Feel free to contact Fan Yin via email, yinf2@uci.edu, if you have any questions.


##### 
# simulation study 
source("./simulation/lysozyme.simulation.study.lognormal.R")
source("./simulation/lysozyme.simulation.study.noadjustment.R")

# summarize the results of simulation study
source("./simulation/lysozyme.simulation.summary.R")


#####
# case study: lysozyme
source("./lysozyme/lysozyme.main.R")


#####
# case study: gamma S crystallin
source("./gsc/gsc.main.R")
