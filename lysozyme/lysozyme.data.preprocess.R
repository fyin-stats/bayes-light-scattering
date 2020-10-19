#Analysis script for the three-detector light scattering data (lysozyme experiments)
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  try(sapply(pkg, require, character.only = TRUE), silent = TRUE)
}
packages <- c("foreach", "doParallel", 
              "boot", "parallel","rjags", "ggplot2", 
              "dplyr", "invgamma", "R2jags",
              "R2WinBUGS","xtable")
ipak(packages)
mc.cores<-25
######################################################
#Load the source file
source("./master.scripts/photoFinish.R")
#setwd("E:/UCI/Research/light_scattering/")
#Define the design concentrations, in g/mL
con.design<-list(
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  # c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  # c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3,
  c(0,2.5e-3,5e-3,7.5e-3,1e-2,1.25e-2,1.5e-2,1.75e-2,2e-2,2.5e-2,3e-2,3.5e-2,4e-2,4.5e-2,5e-2),
  # c(0,2.5e-3,5e-3,7.5e-3,1e-2,1.25e-2,1.5e-2,1.75e-2,2e-2,2.5e-2,3e-2,3.5e-2,4e-2,4.5e-2,5e-2),
  c(0,2.5e-3,5e-3,7.5e-3,1e-2,1.25e-2,1.5e-2,1.75e-2,2e-2,2.5e-2,3e-2,3.5e-2,4e-2,4.5e-2,5e-2),
  c(0,2.5e-3,5e-3,7.5e-3,1e-2,1.25e-2,1.5e-2,1.75e-2,2e-2,2.5e-2,3e-2,3.5e-2,4e-2,4.5e-2,5e-2),
  c(0,2.5e-3,5e-3,7.5e-3,1e-2,1.25e-2,1.5e-2,1.75e-2,2e-2,2.5e-2,3e-2,3.5e-2,4e-2,4.5e-2,5e-2)
)
#Monomer mass, in Da 
monomer.mass<-14307
#Solvent base RI (n0); stored on a per-experiment basis (but a reference value is all
#that was available in some cases)
solvent.base.ri.mean<-rep(1.3244,14)
solvent.base.ri.meas<-c(1.3331, 1.3289, 1.3279, 
                        1.3272, 1.3272, 1.3310, 
                        1.3309, 1.3295, 1.3296, 
                        1.3276, 1.3285, 1.3273, 
                        1.3291, 1.3312)
# solvent.base.ri.meas.sd<-solvent.base.ri.meas*c(rep(0.001,17),rep(1.01,9))/1.96
solvent.base.ri.meas.sd <- sd(solvent.base.ri.meas)
solvent.base.ri.lsd<-log(1.05)/1.96
#RI increment prior
prior.dndc.mean<- 0.197 # 0.197, from paper Protein refractive index increment is determined by conformation as well
#as composition
prior.dndc.sd<- 0.002
#Prior precision of the concentration measurements (and other stuff related)
pathlen<-1       #Path length in cm
extcoef<-2.64    #Extinction coefficient for lysozyme in mL/(mg*cm)
specacc<-0.005   #Spectrophotometer accuracy (+/- absorbance interval)
prior.con.meas.lsd<-log(1+specacc)/1.96  #Effective standard deviation
prior.con.meas.lsd.wgt<-7      #Enough weight to keep SD w/in factor of 2 of prior mean
prior.con.lsd<-log(1.05)/1.96      #Usually within apx 10% of intended concentration
prior.con.lsd.wgt<-2           #We don't put too much prior weight on this
#Priors for the radius of gyration.  For a small globular protein, this should be
#extremely small; for lysozyme, at least one source estimates it at approximately
#1.3nm (Voets et al.; J Phys Chem B; 2010).  We set the prior mean at 1nm, with a prior
#that is relatively diffuse.
prior.Rg.lmean<-log(1e-9)
prior.Rg.lsd<-log(10)/1.96
#Various MCMC settings
# similar results with cheaper MCMC settings?
draws<-100000
thin<- 250
burn<-100000
nadapt<-100000
nchains <- 5
#########################
wl <- 658e-7
f<-(4*pi^2)/(wl^4*6.022e23)  
#Data files and conditions
# conditions
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
#
tsdat<-c(
  "./lysozyme/16.angle/high.conc/15_9_2_lysozyme_250mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_8_31_lysozyme_125mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_8_25_lysozyme_125mMNaCl_pH4.7_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_27_lysozyme_75mMNaCl_pH4.7_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_21_lysozyme_50mMNaCl_pH4.7_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_14_lysozyme_300mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_09_lysozyme_150mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_08_lysozyme_75mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_06_lysozyme_300mMNaCl_pH4.7_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_7_01_lysozyme_150mMNaCl_pH4.7_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_5_22_lysozyme_100mMNaCl_pH4.7_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_5_19_lysozyme_100mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_5_15_lysozyme_50mMNaCl_pH6.9_MALS+RI.csv",
  "./lysozyme/16.angle/high.conc/15_5_14_lysozyme_200mMNaCl_pH6.9_MALS+RI.csv")
cdat<-c(
  "./lysozyme/16.angle/high.conc/15_9_2_lysozyme_250mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_8_31_lysozyme_125mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_8_25_lysozyme_125mMNaCl_pH4.7_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_27_lysozyme_75mMNaCl_pH4.7_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_21_lysozyme_50mMNaCl_pH4.7_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_14_lysozyme_300mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_09_lysozyme_150mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_08_lysozyme_75mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_06_lysozyme_300mMNaCl_pH4.7_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_7_01_lysozyme_150mMNaCl_pH4.7_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_5_22_lysozyme_100mMNaCl_pH4.7_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_5_19_lysozyme_100mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_5_15_lysozyme_50mMNaCl_pH6.9_dndc.csv",
  "./lysozyme/16.angle/high.conc/15_5_14_lysozyme_200mMNaCl_pH6.9_dndc.csv")
#
suppcon<-read.csv("./lysozyme/master_concentration_list.csv", #Load supplemental concentration data 
                  header=TRUE,as.is=TRUE)                         #for the UCIHC set
### 
suppcon[,-(1:3)]<-suppcon[,-(1:3)]*1e-3 #Fix the units
suppcon.cols<-18:31                     #Cols w/the secondary concentration measurements
NaClmM<-c(250,125,125,75,50,300,150,75,300,150,100,100,50,200)
pH<-c(6.9,6.9,4.7,4.7,4.7,6.9,6.9,6.9,4.7,4.7,4.7,6.9,6.9,6.9)
replicate <- rep(1, 14)
#
dset <- rep("UCIHC", 14)
tokeep<-list(     #Columns to keep (weeding out bad detectors and low angles)
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34),
  c(11:18,21:28,33:34)
)
###
trim.tol<-c(0,0.5,0,
            1,0.5,0.5,
            1,0,1,
            1,
            1,

            0,0,1)  #Trimming tolerance
trim.rismoo<-c(T,T,F,F,F,T,F,T,
               T,T,
               T,
               T,T,F) #Local smoothing for RI stability during trimming
trim.lssmoo<-c(F,T,F,F,
               F,F,F,F,
               F,
               F,
               F,
               F,F,F) #Local smoothing for RI stability during trimming
#
fit.use<-which(dset=="UCIHC")      #When fitting, use these data sets for RI data
RI.max.con<-0.02   #Discard RI for concentrations >this; avoids faulty readings
#Adjust concentrations using secondary measurements; taking the geometric mean
#and using to refine estimate of concentration reliability
#Clean data
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
  #Suckful special-case patches due to one data set w/f'd up terminal sequence
  if((pH[i]==6.9)&&(NaClmM[i]==125)&&(replicate[i]==1)){
    dat[5872:6325,2*(1:7)]<-NA  #These obs weren't being used, but mess up clustering...
  }
  clean.dat[[i]]<-photoFinish(dat,cdat[i],trim.tol=trim.tol[i],trim.locsmoo=c(trim.lssmoo[i],trim.rismoo[i]))
  save(con.design,monomer.mass,clean.dat,tsdat,cdat,NaClmM, pH,i,file="./lysozyme/lysozyme.3d.analysis.cleaned.Rdata")
}
###
load("./lysozyme/lysozyme.3d.analysis.cleaned.Rdata")
adat<-clean.dat
temp<-vector()
for(i in which(dset=="UCIHC")){
  sel<-which((suppcon$Salt==NaClmM[i])&(suppcon$pH==pH[i])& (suppcon$Replicate==replicate[i]))
  if(length(sel)>0){
    temp<-c(temp,log(adat[[i]]$results[-1,1])-log(unlist(suppcon[sel,suppcon.cols])))
    adat[[i]]$results[,1]<-(adat[[i]]$results[,1]*c(0,unlist(suppcon[sel,suppcon.cols])))^0.5
  }
}
prior.con.meas.lsd<-sd(temp)/2
prior.con.meas.lsd.wgt<-length(temp)
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
  ls.ciwd<-c(ls.ciwd,as.vector(adat[[i]]$results[,lscols+2]-adat[[i]]$results[,lscols+1]))
  ri.ciwd<-c(ri.ciwd,adat[[i]]$results[risel,ricol+2]-adat[[i]]$results[risel,ricol+1])
  ls.dsd<-c(ls.dsd,apply(adat[[i]]$results[,lscols],1,sd))
  ls.dsdn<-c(ls.dsdn,rep(length(lscols),NROW(adat[[i]]$results)))
}
## 
# sigma_c ~ inv-gamma(shape=m_c/2, rate=m_c/2 * S_c^2)
# sd means the standard deviation of measurements
ls.prsdest<-((sum((ls.ciwd*ls.cin)/sum(ls.cin))/(1.96*2))^2 + (sum(ls.dsd*ls.dsdn)/sum(ls.dsdn))^2)^0.5 # estimate of standard deviation
# square root of variance,
ri.prsdest<-sum((ri.ciwd*ri.cin)/sum(ri.cin))/(1.96*2)
ls.prsdwgt<-length(ls.cin) # weight, how many data points there are
ri.prsdwgt<-length(ri.cin)
###############################################
lysozyme.df.list <- vector(mode="list", length = length(fit.use))
datcon <- vector()
datscat <- matrix(NA, ncol = length(as.numeric(adat[[i]]$scat.ang)))
n0.meas <- solvent.base.ri.meas # solvent base refractive index measurements
datri <- vector()
#
eda.df <- data.frame(con=NA, A2=NA, dndc=NA, NaCl=NA, pH=NA, con.design=NA)
#################################################
for(i in fit.use){
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
  #browser()
  temp_datri<-adat[[i]]$results[,ricol] # refractive index data
  temp_datri<-(temp_datri-temp_datri[1])[-1]
  temp_dndc <- temp_datri/temp_datcon
  # temp_datri[adat[[i]]$results[-1,1]>RI.max.con]<-NA  #Remove untrustworthy observations, according to the concentration
  datri <- c(datri, temp_datri)
  #print(datri)
  #eda.df <- rbind(eda.df, data.frame(con=temp_datcon, A2=temp_A2,dndc=temp_dndc,NaCl = NaClmM[i], pH = pH[i], 
  #                                   con.design=con.design[[i]][1:length(temp_datcon)]))
  cat("\n")
  cat("\n")
  # browser()
  # remove measurements collected for the case in which the design concentration is zero
  lysozyme.df.list[[i]]$datcon <- temp_datcon
  lysozyme.df.list[[i]]$datscat <- temp_datscat
  lysozyme.df.list[[i]]$datri <- temp_datri
  lysozyme.df.list[[i]]$n0.meas <- solvent.base.ri.meas[i]
  lysozyme.df.list[[i]]$con.design <- con.design[[i]][-1]
  id.ri <- which(!is.na(lysozyme.df.list[[i]]$datri))
  id.scat <- which(!is.na(lysozyme.df.list[[i]]$datscat))
  lysozyme.df.list[[i]]$ri.meas <- lysozyme.df.list[[i]]$datri[id.ri]
  
  print(i)
  # print(lysozyme.df.list[[i]]$ri.meas/prior.dndc.mean)
  temp_c1 <- lysozyme.df.list[[i]]$ri.meas/prior.dndc.mean * monomer.mass
  temp_A2 <- (1-(temp_datscat[,4]/(f*1.32^2*0.21^2*temp_c1) ) )/(2*temp_c1) 
  print(temp_A2)
  #cat("\n")
} # experiment 3 and 4 is more consitent in terms of estimating
# remove the NA row
datscat <- datscat[complete.cases(datscat),]
id.scat <- which(!is.na(datscat[,1]))
#id.ri <- which(!is.na(datri))
ri.meas <- datri[id.ri]
scat.meas <- datscat[id.scat,]
eda.df <- eda.df[complete.cases(eda.df),]
#######################################
lysozyme.df.list <- lysozyme.df.list
# [-c(11,12,14)]
#######################################
##### same design concentration 
#######################################
# data 
# relative scattering intensity
# i,j,l (3D array) # ncol(lysozyme.df.list[[1]]$datscat)
datscat_array <- array(NA, dim = c(nrow(lysozyme.df.list[[1]]$datscat),
                                   1,
                                   length(lysozyme.df.list)))
for(l in 1:length(lysozyme.df.list)){
  datscat_array[,,l] <- lysozyme.df.list[[l]]$datscat[,4]
}
# delta n
# nri,l (2D matrix)
# lapply(lysozyme.df.list, function(x) x$datcon[8])
datri_matrix <- matrix(NA, 
                       nrow = sum(lysozyme.df.list[[1]]$con.design<=0.02), 
                       ncol= length(lysozyme.df.list))
for(l in 1:length(lysozyme.df.list)){
  datri_matrix[,l] <- lysozyme.df.list[[l]]$datri[1:8]
}
# concentration c^m
# i,l (2D matrix)
con_matrix <- matrix(NA, nrow = length(lysozyme.df.list[[1]]$con.design), 
                     ncol = length(lysozyme.df.list))
for(l in 1:length(lysozyme.df.list)){
  con_matrix[,l] <- lysozyme.df.list[[l]]$datcon
}
# con_design : a vector (because all experiments have the same design concentrations)
con_design_vector <-lysozyme.df.list[[1]]$con.design

### gather the data
# dat_jags <- list(scat.meas = datscat_array,
#                  n0.meas=n0.meas,
#                  n0=solvent.base.ri.meas[1:17][-c(11,12,14)],
#                  ri.meas=datri_matrix,
#                  con.meas=con_matrix,
#                  con.design=con_design_vector,
#                  prior.A2.sd=1,
#                  prior.dndc.mu.sd=prior.dndc.sd,
#                  prior.dndc.mu.mean=prior.dndc.mean,
#                  max.oligomer.size=3,
#                  N=0,
#                  L=length(lysozyme.df.list),
#                  nri=dim(datri_matrix)[1],
#                  n=dim(datscat_array)[1],
#                  npred = 100,
#                  f=f,
#                  monomer.mass=monomer.mass,
#                  nang=dim(datscat_array)[2],
#                  a_prior.dndc.s2=2/2,
#                  b_prior.dndc.s2=(2/2)*(0.0002^2),
#                  a_R = 2/2,
#                  b_R = (2/2)*(10^(-10)),
#                  a_n = 2/2,
#                  b_n = (2/2)*(10^(-8)),
#                  a_c = 2/2,
#                  b_c = 2/2*(log(1.05)/1.96)^2 )
#### 