#Validation experiments for the estimation code
#CTB 9/16/15

library(parallel)

#Experimental conditions
solvent.base.ri<-1.33
dndc<-0.188
con.design<-c(0,2.5,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50)*1e-3
A2<-c(-10,5,1,0,1,5,10)*1e-4
M<-c(1,2,3)
monomer.mass<-14307
con.lsd<-log(1.1)
con.meas.lsd<-log(1.1)
scat.sd<-1e-6
ri.sd<-1e-4
scat.ang<-c(60.02,69.37,79.73,90.00,110.63,121.17,132.07,142.37)
cores<-40
draws<-250
thin<-5e3
burn<-1e4

#Priors
prior.scat.sd<-1e-6
prior.scat.sd.wgt<-10
prior.ri.sd<-ri.sd
prior.ri.sd.wgt<-10
prior.con.lsd<-con.lsd
prior.con.lsd.wgt<-10
prior.con.meas.lsd<-con.meas.lsd
prior.con.meas.lsd.wgt<-10
prior.dndc.mean<-dndc
prior.dndc.sd<-0.2
prior.n0.mean<-solvent.base.ri
prior.n0.lsd<-log(1.5)
prior.n0.meas.sd<-ri.sd

valsim<-mclapply(unlist(apply(cbind(rep(M,length(A2)),rep(A2,each=length(M))),1,list),recursive=FALSE),function(z){
  #Simulate data
  dat<-simulate.virialEst(con.design=con.design, solvent.base.ri=solvent.base.ri, monomer.mass=monomer.mass, scat.ang=scat.ang, A2=z[2], aggregate.size=z[1], wl=658e-7, dndc=dndc, con.lsd=con.lsd, con.meas.lsd=con.meas.lsd, scat.sd=scat.sd, ri.sd=ri.sd)
  #Take posterior draws
  fit<-virialEst(con.meas=dat$concentration.measured[-1], scat.meas=dat$scattering.measured[-1,], con.design=dat$concentration.design[-1], ri.meas=dat$RI.measured[-1], n0.meas=dat$base.RI.measured, monomer.mass=monomer.mass, scat.ang=scat.ang, burn=burn, thin=thin, draws=draws, prior.scat.sd=prior.scat.sd, prior.scat.sd.wgt=prior.scat.sd.wgt, prior.ri.sd=prior.ri.sd, prior.ri.sd.wgt=prior.ri.sd.wgt, prior.con.lsd=prior.con.lsd, prior.con.lsd.wgt=prior.con.lsd.wgt, prior.con.meas.lsd=prior.con.meas.lsd, prior.con.meas.lsd.wgt=prior.con.meas.lsd.wgt, prior.dndc.mean=prior.dndc.mean, prior.n0.mean=prior.n0.mean, prior.n0.lsd=prior.n0.lsd, prior.n0.meas.sd=prior.n0.meas.sd, prior.dndc.sd=prior.dndc.sd,max.oligomer.size=4)
  #Return the results
  list(dat=dat,fit=fit)
},mc.cores=cores)

#BELOW: Old stuff from earlier versions of the code/data, left for historical reasons.

#Example of a fit that didn't work well....
simdat2<-simulate.virialEst(c(0,2.5e-4,5e-4,7.5e-4,1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3,(6:10)*1e-3), 1.33411, monomer.mass=14307, scat.ang=c(40,90,138), A2=3e-4)

simfit2<-virialEst(con.meas=simdat2$concentration.measured,scat.meas=simdat2$scattering.measured,con.design=simdat2$concentration.design,solvent.base.ri=simdat2$solvent.base.RI,ri.meas=simdat2$RI.measured,monomer.mass=simdat2$monomer.mass,scat.ang=simdat2$scat.ang,burn=2e5,thin=5e3,draws=250,prior.scat.sd=1e-6,prior.scat.sd.wgt=2,prior.ri.sd=1e-7,prior.ri.sd.wgt=2,prior.con.lsd=log(1.1),prior.con.lsd.wgt=2,prior.con.meas.lsd=log(1.005),prior.con.meas.lsd.wgt=2,prior.dndc.mean=c(0,0.18),prior.dndc.sd=0.3,prior.Rg.lmean=log(1e-9),prior.Rg.lsd=log(10))


#Let's see what's needed to get signal.  Answer: much higher concentrations!
dR<-vector()
a2<-seq(from=-5e-4,to=5e-4,length=11)
for(i in 1:length(a2))
  dR<-rbind(dR,rowMeans(simulate.virialEst(c(0,2.5e-4,5e-4,7.5e-4,1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3,(6:10)*1e-3,(2:5)*1e-2), 1.33411, monomer.mass=14307, scat.ang=c(40,90,138), A2=a2[i])$scattering.measured))
plot(0,0,type="n",xlim=c(0,5e-2),ylim=range(dR),xlab="Concentration (g/mL)",ylab=expression(Delta*R(theta)))
for(i in 1:length(a2))
  points(c(0,2.5e-4,5e-4,7.5e-4,1e-3,1.5e-3,2e-3,3e-3,4e-3,5e-3,(6:10)*1e-3,(2:5)*1e-2),dR[i,], col=hsv((i-1)/(length(a2)-1)*0.6),type="b")
legend("topleft",lty=1,col=hsv(((1:length(a2))-1)/(length(a2)-1)*0.6),title=expression(A[2]),legend=round(a2,5))

#Let's see what happens now...
simdat<-simulate.virialEst(c(0,1e-3,5e-3,1e-2,2e-2,3e-2,4e-2,5e-2), 1.33411, monomer.mass=14307, scat.ang=c(40,90,138), A2=1e-4)
simdat2<-simulate.virialEst(c(0,1e-3,5e-3,1e-2,2e-2,3e-2,4e-2,5e-2), 1.33411, monomer.mass=14307, scat.ang=c(40,90,138), A2=-1e-4)

simfit2<-virialEst(con.meas=simdat2$concentration.measured,scat.meas=simdat2$scattering.measured,con.design=simdat2$concentration.design,solvent.base.ri=simdat2$solvent.base.RI,ri.meas=simdat2$RI.measured,monomer.mass=simdat2$monomer.mass,scat.ang=simdat2$scat.ang,burn=1e5,thin=1e3,draws=250,prior.scat.sd=1e-6,prior.scat.sd.wgt=2,prior.ri.sd=1e-7,prior.ri.sd.wgt=2,prior.con.lsd=log(1.1),prior.con.lsd.wgt=2,prior.con.meas.lsd=log(1.005),prior.con.meas.lsd.wgt=2,prior.dndc.mean=c(0,0.18),prior.dndc.sd=0.3,prior.Rg.lmean=log(1e-9),prior.Rg.lsd=log(10), con.prop.lsd=log(1.001))

a2<-seq(from=-5e-4,to=5e-4,length=5)
a2.est<-vector()
dndc.est<-vector()
m.est<-vector()
for(i in 1:length(a2)){
  cat("a2:",a2[i],"\n")
  simdat<-simulate.virialEst(c(0,1e-3,5e-3,1e-2,1.5e-2,2e-2,2.5e-2,3e-2,3.5e-2, 4e-2,4.5e-2, 5e-2), 1.33411, monomer.mass=14307, scat.ang=c(40,90,138), A2=a2[i], con.lsd=log(1.025), con.meas.lsd=log(1.0025), scat.sd=1e-9, ri.sd=1e-7)
  simfit<-virialEst(con.meas=simdat$concentration.measured,scat.meas=simdat$scattering.measured,con.design=simdat$concentration.design, solvent.base.ri=simdat$solvent.base.RI, ri.meas=simdat$RI.measured,monomer.mass=simdat$monomer.mass, scat.ang=simdat$scat.ang,burn=5e5, thin=5e3,draws=250,prior.scat.sd=1e-9, prior.scat.sd.wgt=10, prior.ri.sd=1e-7,prior.ri.sd.wgt=10, prior.con.lsd=log(1.025),prior.con.lsd.wgt=2, prior.con.meas.lsd=log(1.0025), prior.con.meas.lsd.wgt=10, prior.dndc.mean=c(0,0.18), prior.dndc.sd=0.3, prior.Rg.lmean=log(1e-9), prior.Rg.lsd=log(10), con.prop.lsd=log(1.005),prior.agg.lambda=0.0001,max.oligomer.size=1)
  a2.est<-rbind(a2.est,quantile(simfit$A2,c(0.025,0.25,0.5,0.75,0.975)))
  dndc.est<-rbind(dndc.est,quantile(simfit$RI.reg.param[,2],c(0.025,0.25,0.5,0.75,0.975)))
  m.est<-rbind(m.est,quantile(simfit$aggregate.size,c(0.025,0.25,0.5,0.75,0.975)))
}

#Function to infer A2 directly in the ideal case...
#Recall: dR = K/(1/Mc + 2A)
#        (K/dR-1/Mc)/2 = A
a2fun<-function(n0,dndc,dR,wl,sa,Rg,m,W,c){
  ((2*pi*n0*dndc)^2/((dR*wl^4*6.022e23)*(1+(4*pi*n0*sin(sa/360*pi)/wl)^2*Rg^2/3))-1/(m*W*c))*0.5
}

#Examination of detector-level errors in the lysozyeme data
load("lysozyme/lysozyme.3d.analysis.cleaned.Rdata")
detectDetect<-function(x,logscale=FALSE,suppress.zero=FALSE,plot.ci=TRUE,...){
  sa<-as.numeric(x$scat.ang)
  nang<-length(sa)
  lscol<-6+((1:nang)-1)*4
  if(logscale){
    lsdat<-log(x$results[,lscol])
    lsub<-log(x$results[,lscol+1])
    lslb<-log(x$results[,lscol+2])
  }else{
    lsdat<-x$results[,lscol]
    lsub<-x$results[,lscol+1]
    lslb<-x$results[,lscol+2]
  }
  if(suppress.zero){
    lsdat<-lsdat[-1,]
    lsub<-lsub[-1,]
    lslb<-lslb[-1,]
  }
  lsmean<-rowMeans(lsdat)
  lserr<-sweep(lsdat,1,lsmean,"-")
  lsub<-sweep(lsub,1,lsmean,"-")
  lslb<-sweep(lslb,1,lsmean,"-")
  if(logscale)
    ylab<-"log(R(theta,c))-mean(log(R(theta,c)))"
  else
    ylab<-"R(theta,c)-mean(R(theta,c))"
  plot(90,0,type="n",xlab="Detector Angle",ylab=ylab,xlim=range(sa),ylim=range(lserr),...)
  abline(h=0,lty=3)
  pal<-hsv(((1:NROW(lsdat))-1)/(NROW(lsdat)-1)*0.6)
  for(i in 1:NROW(lsdat)){
    lines(sa,lserr[i,],col=pal[i],type="b")
    if(plot.ci)
      segments(sa,lslb[i,],sa,lsub[i,],col=pal[i])
  }
  conc<-round(x$results[,1],4)
  if(suppress.zero)
    conc<-conc[-1]
  legend("bottomright",lty=1,col=pal,pch=1,legend=paste(conc,x$conc.units))
}
pdf("lysozyme.meas.err.pdf",12,12)
par(mfrow=c(3,3))
for(i in 1:length(clean.dat))
  detectDetect(clean.dat[[i]],main=tsdat[i]) #See lysozyme.3d.analysis.R
for(i in 1:length(clean.dat))
  detectDetect(clean.dat[[i]],logscale=TRUE,suppress.zero=TRUE,main=tsdat[i])
dev.off()


#Examination of naive predictions in the lysozyme data, to get a sense of detection issues
crudeA2<-function(x,conc.corrected=TRUE,n0=1.33411,dndc=0.181,wl=658e-7,Rg=1.3e-7,W=14307,legend.pos="right", ...){
  sa<-as.numeric(x$scat.ang)
  nang<-length(sa)
  lscol<-6+((1:nang)-1)*4
  lsdat<-x$results[,lscol]
  lsdat<-sweep(lsdat,2,lsdat[1,],"-")[-1,]
  conc<-x$results[-1,1+conc.corrected]
  a2pred<-vector()
  for(i in 1:nang)
    a2pred<-rbind(a2pred,a2fun(n0=n0,dndc=dndc,dR=lsdat[,i],wl=wl,sa=sa[i],Rg=Rg,m=1,W=W,c=conc))
  plot(1,1,type="n",xlim=range(conc),ylim=range(a2pred),xlab=paste("Concentration (",x$conc.units,")",sep=""),ylab=expression(hat(A)[2]),...)
  abline(h=0,lty=3)
  pal<-hsv(((1:nang)-1)/(nang-1)*0.6)
  for(i in 1:nang)
    lines(conc,a2pred[i,],type="b",col=pal[i])
  legend(legend.pos,lty=1,pch=1,col=pal,legend=x$scat.ang)
  invisible(list(conc=conc,A2=a2pred))
}


