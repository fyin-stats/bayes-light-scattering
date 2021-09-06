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
# setwd("./gsc/")
source("./master.scripts/photoFinish.R")
#Define the design concentrations, in g/mL
con.design<-list(
  c(0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30)*1e-3,
  c(0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30)*1e-3,
  c(0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3
) 
con.design.vector.list<-list(
  c(0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30)*1e-3,
  c(0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30)*1e-3,
  c(0.5,1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30)*1e-3
) 
#Monomer mass, in Da
monomer.mass<-20959.80
#Solvent base RI (n0); stored on a per-experiment basis (but a reference value is all
#that was available in some cases)
solvent.base.ri.mean<-rep(1.3244,3)
solvent.base.ri.meas<-c(1.3265,1.3286,1.3248)
solvent.base.ri.meas.sd<-sd(solvent.base.ri.meas)
# solvent.base.ri.meas.sd<-solvent.base.ri.meas*rep(0.001,5)/1.96
solvent.base.ri.lsd<-log(1.05)/1.96
#RI increment prior
prior.dndc.mean<-0.1985
prior.dndc.sd<-0.005 # how do we dertermine these two values?
#Prior precision of the concentration measurements (and other stuff related)
pathlen<-1       #Path length in cm
extcoef<-NA      #Extinction coefficient for gSC in mL/(mg*cm) (NA b/c I don't know it)
specacc<-0.005   #Spectrophotometer accuracy (+/- absorbance interval)
prior.con.meas.lsd<-log(1+specacc)/1.96  #Effective standard deviation
prior.con.meas.lsd.wgt<-7      #Enough weight to keep SD w/in factor of 2 of prior mean
prior.con.lsd<-log(1.05)/1.96       #Usually within apx 10% of intended concentration
prior.con.lsd.wgt<-2           #We don't put too much prior weight on this
#Various MCMC settings
draws<-100000
thin<-250
burn<-100000
nadapt<-100000
nchains<-5
#
wl <- 658e-7
f<-(4*pi^2)/(wl^4*6.022e23)  
# Data files and conditions  # a total of five conditions
# measurement data and corresponding concentration data
# what is the difference between them?
tsdat<-c(
  "./gsc/raw.data/16_02_04_gSWT_100mMNaCl_pH6.9_MALS.csv",
  "./gsc/raw.data/16_01_26_gSWT_100mMNaCl_pH6.9_MALS.csv",
  "./gsc/raw.data/16_01_21_gSWT_100mMNaCl_pH6.9_MALS.csv"
)
cdat<-c(
  "./gsc/raw.data/16_02_04_gSWT_100mMNaCl_pH6.9_conc.csv",
  "./gsc/raw.data/16_01_26_gSWT_100mMNaCl_pH6.9_conc.csv",
  "./gsc/raw.data/16_01_21_gSWT_100mMNaCl_pH6.9_conc.csv"
)
NaClmM<-c(100,100,100)
pH<-c(6.9,6.9,6.9)
replicate<-c(1,1,1)
dset<-rep(c("UCIHC"),3) # 
tokeep<-list(     #Columns to keep (weeding out bad detectors and low angles)
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34)
)
trim.tol<-c(1,1,1)  #Trimming tolerance
trim.rismoo<-c(T,T,T) #Local smoothing for RI stability during trimming
trim.lssmoo<-c(T,T,T) #Local smoothing for RI stability during trimming
fit.use<-which(dset=="UCIHC")      #When fitting, use these data sets for RI data
RI.max.con<-0.02   #Discard RI for concentrations >this; avoids faulty readings
##########################################################
# Clean data
clean.dat<-vector(mode="list",length=length(tsdat))
for(i in 1:length(tsdat)){
  cat("Working on NaCl concentration",NaClmM[i],"pH",pH[i],"\n")
  cat("Cleaning data...\n")
  dat<-read.csv(tsdat[i],header=T,as.is=T,fileEncoding="Latin1")[,tokeep[[i]]]
  dup<-which(diff(dat[,1])==0)  #deal with bizarre LS time duplication issue
  if(any(dup)){
    nondup<-dat[-dup,1:(NCOL(dat)-2)]
    dat[1:NROW(nondup),1:(NCOL(dat)-2)]<-nondup
    dat[(NROW(nondup)+1):NROW(dat),1:(NCOL(dat)-2)]<-NA
    dat<-dat[!apply(is.na(dat),1,all),]
  }
  clean.dat[[i]]<-photoFinish(dat,cdat[i],trim.tol=trim.tol[i],trim.locsmoo=c(trim.lssmoo[i],trim.rismoo[i]))
  print(clean.dat[[i]])
  save(con.design,monomer.mass,clean.dat,tsdat,cdat,NaClmM,dset,replicate, pH,i,file="./gsc/gsc.analysis.cleaned.Rdata")
}
##########################################################
#### read the cleaned data
load("./gsc/gsc.analysis.cleaned.Rdata")
#Adjust concentrations using secondary measurements; taking the geometric mean
#and using to refine estimate of concentration reliability
adat<-clean.dat
temp<-vector()
for(i in which(dset=="UCIHC")){
  dat<-read.csv(cdat[i],header=T,as.is=T,fileEncoding="Latin1")
  # print(dat)
  dat <- dat[complete.cases(dat),]
  sel<-1:(NROW(adat[[i]]$results)-1)
  temp<-c(temp,log(adat[[i]]$results[-1,1])-log(dat[sel,2])) #It's in the second col
  adat[[i]]$results[,1]<-(adat[[i]]$results[,1]*c(0,dat[sel,2]))^0.5
}
#prior.con.meas.lsd<-sd(temp)/2 #s_cm = 0.205
#prior.con.meas.lsd.wgt<-length(temp)  # m_cm = 68

#Use combined RI and scattering measurements to set prior SDs
ls.cin<-vector()      #LS CI Ns
ls.ciwd<-vector()     #LS CI widths
ls.dsdn<-vector()     #LS detector SD Ns
ls.dsd<-vector()      #LS detector SDs
ri.cin<-vector()
ri.ciwd<-vector()
for(i in fit.use){
  scat.ang<-as.numeric(adat[[i]]$scat.ang)
  lscols<-6+((1:length(scat.ang))-1)*4 #LS measurement cols
  ricol<-max(lscols)+6                 #RI measurement cols
  risel<-adat[[i]]$results[,1]<RI.max.con
  ls.cin<-c(ls.cin,as.vector(adat[[i]]$results[,lscols-1]))
  ri.cin<-c(ri.cin,as.vector(adat[[i]]$results[risel,ricol-1]))
  ls.ciwd<-c(ls.ciwd,as.vector(adat[[i]]$results[,lscols+2]-adat[[i]]$results[,lscols+1])) # width of ci
  ri.ciwd<-c(ri.ciwd,adat[[i]]$results[risel,ricol+2]-adat[[i]]$results[risel,ricol+1])
  ls.dsd<-c(ls.dsd,apply(adat[[i]]$results[,lscols],1,sd))
  ls.dsdn<-c(ls.dsdn,rep(length(lscols),NROW(adat[[i]]$results)))
}
ls.prsdest<-((sum((ls.ciwd*ls.cin)/sum(ls.cin))/(1.96*2))^2 + (sum(ls.dsd*ls.dsdn)/sum(ls.dsdn))^2)^0.5  # S_R
ri.prsdest<-(sum((ri.ciwd*ri.cin)/sum(ri.cin))/(1.96*2)) # S_n
ls.prsdwgt<-length(ls.cin) # m_R
ri.prsdwgt<-length(ri.cin) # m_n
# combine data from multiple experiments together
#################################################
datcon <- vector()
datscat <- matrix(NA, ncol = length(as.numeric(adat[[i]]$scat.ang)))
n0.meas <- solvent.base.ri.meas # solvent base refractive index measurements
datri <- vector()
eda.df <- data.frame(con=NA, A2=NA, dndc=NA, NaCl=NA, pH=NA, con.design=NA)
#################################################
# fit.use
fit.use <- c(1,2,3)
gsc.df.list <- vector(mode="list", length = length(fit.use))
for(j in 1:length(fit.use)){
  i = fit.use[j]
  cat("Working on NaCl concentration",NaClmM[i],"pH",pH[i],"\n")
  cat("\tSetting up data...\n")
  scat.ang<-as.numeric(adat[[i]]$scat.ang) # angles
  lscols<-6+((1:length(scat.ang))-1)*4 #LS measurement cols, LS means light scattering
  ricol<-max(lscols)+6                 #RI measurement cols, RI means refractive index
  # browser()
  # scattering data
  temp_datscat<-sweep(adat[[i]]$results[,lscols],2,adat[[i]]$results[1,lscols],"-") # subtract reference row from each row
  temp_datscat<-temp_datscat[-1,] # remove the first row, which is the reference row, light scattering intensity data
  datscat <- rbind(datscat, temp_datscat)
  #print(datscat)
  # concentration measurements
  temp_datcon <- adat[[i]]$results[-1,1]
  datcon <- c(datcon, temp_datcon)
  # Approximate A2
  # ri increment measurements
  temp_datri<-adat[[i]]$results[,ricol] # refractive index data
  temp_datri<-(temp_datri-temp_datri[1])[-1]
  temp_dndc <- temp_datri/temp_datcon
  # temp_datri[adat[[i]]$results[-1,1]>RI.max.con]<-NA  #Remove untrustworthy observations, according to the concentration
  datri <- c(datri, temp_datri)
  cat("\n")
  cat("\n")
  # remove measurements collected for the case in which the design concentration is zero
  gsc.df.list[[j]]$datcon <- temp_datcon[1:14]
  gsc.df.list[[j]]$datscat <- temp_datscat[1:14,]
  gsc.df.list[[j]]$datri <- temp_datri[1:12]
  gsc.df.list[[j]]$n0.meas <- solvent.base.ri.meas[i]
  gsc.df.list[[j]]$con.design <- con.design[[i]][1:14]
  id.ri <- which(!is.na(gsc.df.list[[j]]$datri))
  id.scat <- which(!is.na(gsc.df.list[[j]]$datscat))
  gsc.df.list[[j]]$ri.meas <- gsc.df.list[[j]]$datri[id.ri]
  print(i)
} # experiment 3 and 4 is more consitent in terms of estimating
# remove the NA row
datscat <- datscat[complete.cases(datscat),]
id.scat <- which(!is.na(datscat[,1]))
id.ri <- which(!is.na(datri))
ri.meas <- datri[id.ri]
scat.meas <- datscat[id.scat,]
con.design.vector <- do.call("c", con.design.vector.list[fit.use])
##############################################
##############################################
####################### same design concentration 
# data 
# relative scattering intensity
# i,j,l (3D array)
datscat_array <- array(NA, dim = c(nrow(gsc.df.list[[1]]$datscat),
                                   1,
                                   length(gsc.df.list))) # ncol(gsc.df.list[[1]]$datscat)
for(l in 1:length(gsc.df.list)){
  datscat_array[,1,l] <- gsc.df.list[[l]]$datscat[,4]
}
# lapply(lysozyme.df.list, function(x) x$datcon[8])
datri_matrix <- matrix(NA, 
                       nrow = sum(gsc.df.list[[1]]$con.design<=0.02), 
                       ncol= length(gsc.df.list))
for(l in 1:length(gsc.df.list)){
  datri_matrix[,l] <- gsc.df.list[[l]]$datri[1:12]
}
# concentration c^m
# i,l (2D matrix)
con_matrix <- matrix(NA, nrow = length(gsc.df.list[[1]]$con.design), 
                     ncol = length(gsc.df.list))
for(l in 1:length(gsc.df.list)){
  con_matrix[,l] <- gsc.df.list[[l]]$datcon
}
# con_design : a vector (because all experiments have the same design concentrations)
con_design_vector <-gsc.df.list[[1]]$con.design

## save the data
############################################
# saveRDS(gsc_dat_jags, "./gsc/gsc_dat_jags.rds")
