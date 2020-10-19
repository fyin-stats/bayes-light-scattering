#photoFinish - A simple tool for processing light scattering data
#Carter T. Butts, 2/18/16

#A simple tool for processing light scattering data.  photoFinish takes as input two
#respective files with time series data on scattering and refractive index and
#concentration.  This information is merged to identify time points whose concentration
#and scattering/RI measurements are stable and well-characterized; these are clustered
#to produce a quality-trimmed sample in each concentration regime.  Median measurements
#are computed for each measurement in each regime, along with associated statistics 
#(e.g., confidence intervals and the like).  The resulting object (of class
#photoFinish) has a basic print method, as well as a plot method that can produce
#various summary and/or diagnostic displays.
#  The primary virtue of photoFinish is that it removes experimental artifacts 
#related to changes in concentration and to the occasional presence of gas bubbles
#(which can wreak havoc on the results).  For the most part, this is done via the
#above-mentioned quality trimming; however, the use of robust central tendency
#estimation is also important for minimizing the leverage of occasional "bubble
#observations" that make it past the trimming process.  The plot method can be used
#to more closely examine the points that were or were not removed.  In cases for which
#refractive index and scattering data have been taken from the same sample, the latter
#can also be used to correct the former (which may improve results for some analyses).
#  Note that the clustering method used to process the data is fairly slow, and can
#take a while.  The solution does have good properties, however, so bear with it.  Also,
#note that concentration corrections are computed whether or not the RI and scattering
#data was taken on the same sample (since the program doesn't know that); these should
#be ignored where not applicable.
#
#  Arguments:
#    tsfile - the name of the file containing the time series data.
#    confile - the name of the file containing the concentration data.
#    trim.tol - the maximum number of IQRs above the median that the absolute
#      derivative of an observation is allowed to be before being trimmed; by
#      default, any point more than 1 IQR above the median is trimmed.  Smaller
#      values will produce more stable results, at the cost of excluding more data.
#    trim.precombine - logical; precombine series before assessing whether points
#      should be trimmed?  (Defaults to FALSE)
#    trim.locsmooth - logical vector; should derivative estimates (for light scattering
#      and RI data series, respectively) be smoothed prior to assessing whether points
#      should be trimmed?  This can be useful in removing point pairs or triplets
#      in unstable regions of the time series that would otherwise produce spurious
#      data clusters.  (Defaults to c(FALSE,FALSE))
#    time.units - the units in which temporal measurements were given (defaults to
#      minutes).
#    conc.units - the units in which concentration measurements were given (defaults
#      to g/mL).
#    scat.ang - the angles at which scattering was measured (defaults to the first
#      string in each applicable column header after a space).
#    scat.cols - data columns containing scattering information (defaults to every
#      second column up to the last two).
#    intermeasure.flush - logical; was the sample flushed out between each measurement?
#      (Defaults to FALSE)
#    verbose - logical; should progress messages be printed?  (Defaults to TRUE)
#
#  Return value:
#    An object of class "photoFinish."  Print and plot methods are currently
#    implemented.  Important elements of the returned object include:
#      results - table containing summary results
#      ols.coefs - table containing OLS coefficients for each measurement median vs
#        concentration
#      param.est - table containing light scattering parameter estimates based on the
#        model 1/R = f/M 1/c + 2 f A2 (where M is the total mass of the scattering
#        particle, c is the concentration, R is the scattering measurement,
#        A2 is the second virial coefficient, and f is a constant (depending on
#        the experimental conditions and protein).  Estimation is performed using
#        a linear model on R'=1/R and c'=1/c, which is not very stable; robust
#        fits and fits using corrected concentrations are provided, which may help.
#      include.LS, include.RI - vectors indicating untrimmed points
#      clust.LS, clust.RI - vectors indicating point cluster membership
#      time.LS, time.RI - vectors containing raw time points
#      vals.LS - time by angle matrix of raw scattering observations
#      vals.RI - time by dRI vector of raw refractive index observations
#      instability.LS, instability.RI - instability indices (used in trimming)
#      conc.corrected - dRI-corrected concentration measurements
#      conc.correction.model - model used to produce the corrections
#
photoFinish<-function(tsfile,confile, trim.tol=1,trim.precombine=FALSE, trim.locsmooth=c(FALSE,FALSE), time.units="min", conc.units="g/mL", scat.ang=NULL, scat.cols=NULL, intermeasure.flush=FALSE, verbose=TRUE){
  require(robustbase)
  #Read the data
  if((is.character(tsfile)||is.character(confile))&&verbose)
      cat("Reading data files...\n")
  if(is.character(tsfile))
    tdat<-read.csv(tsfile,header=TRUE,as.is=TRUE,fileEncoding="latin1")
  else
    tdat<-tsfile
  if(is.character(confile))
    cdat<-read.csv(confile,header=TRUE,as.is=TRUE,fileEncoding="latin1")
  else
    cdat<-confile
  #Extract the time points for the time series data, and other useful info
  if(verbose)
    cat("Identifying and trimming unstable observations...\n")
  stimes<-tdat[,1]        #Scattering time points
  stimes<-stimes[!is.na(stimes)]
  if(is.null(scat.cols)){
    nang<-NCOL(tdat)/2-1  #Number of angles
    angcol<-(1:nang)*2
    ricol<-max(angcol)+2
  }else{
    nang<-length(scat.cols)
    angcol<-scat.cols
    ricol<-max(angcol)+2
  }
  if(is.null(scat.ang))   #Names of the angles (for plots and such)
    scat.ang<-sapply(strsplit(colnames(tdat)[angcol],"[.][.]"),"[[",2)
  sa<-tdat[!is.na(tdat[,1]),angcol]       #Scattering at each angle
  sa<-apply(sa,c(1,2),as.numeric)
  rtimes<-tdat[,ricol-1]  #Refractive index times
  dri<-tdat[,ricol]       #dRI values
  dri<-dri[!is.na(rtimes)]
  rtimes<-rtimes[!is.na(rtimes)]
  #Determine which points should be trimmed due to instability
  #  First, estimate the absolute time derivatives for each measurement
  #  (We do this in a slightly non-traditional way, given our purpose and the
  #  nature of the data.)
  if(trim.precombine){
    d<-apply(sa,1,median)
    d<-abs(diff(d))/diff(stimes)
    dsdt<-rowMeans(cbind(c(NA,d),c(d,NA)),na.rm=TRUE)
  }else{
    d<-rowMeans(abs(apply(sa,2,diff)),na.rm=TRUE)/diff(stimes)
    dsdt<-rowMeans(cbind(c(NA,d),c(d,NA)),na.rm=TRUE)
  }
  d<-abs(diff(dri))/diff(rtimes)
  drdt<-rowMeans(cbind(c(NA,d),c(d,NA)),na.rm=TRUE)
  if(trim.locsmooth[1])
    dsdt<-rowMeans(cbind(c(NA,dsdt[-length(dsdt)]),dsdt,c(dsdt[-1],NA)),na.rm=TRUE)
  if(trim.locsmooth[2])
    drdt<-rowMeans(cbind(c(NA,drdt[-length(drdt)]),drdt,c(drdt[-1],NA)),na.rm=TRUE)
  #  Now, keep observations no more than trim.tol IQRs above the median
  keeps<-(dsdt-median(dsdt,na.rm=TRUE))/IQR(dsdt,na.rm=TRUE)<trim.tol
  keeps[is.na(keeps)]<-FALSE
  keepr<-(drdt-median(drdt,na.rm=TRUE))/IQR(drdt,na.rm=TRUE)<trim.tol
  keepr[is.na(keepr)]<-FALSE
  #Pull the concentration data, and find out how many clusters we should have
  conc<-cdat[,1]
  conc<-conc[!is.na(conc)]         #Remove any excess lines
  if(intermeasure.flush)
    nc<-2*length(conc)+1           #One cluster/concentration plus interstitial/pre/post
  else
    nc<-length(conc)+2             #One cluster/concentration value plus pre/post
  #Perform a constrained agglomerative clustering to find the intervals associated with
  #each concentration value (plus the initial/final baselines)
  if(verbose)
    cat("Clustering remaining observations (may take a while)...\n")
  kst<-stimes[keeps]                    #Kept times for scattering data
  krt<-rtimes[keepr]                    #Kept times for refraction data
  ksval<-sa[keeps,]                     #Kept scattering vals 
  krval<-dri[keepr]                     #Kept refraction vals 
  scmemb<-constrained.clust(ksval,nc)   #Cluster the scattering data
  rcmemb<-constrained.clust(krval,nc)   #Cluster the refraction data
  scord<-(1:nc)[order(sapply(1:nc,function(cl){median(kst[scmemb==cl],na.rm=TRUE)}))]
  rcord<-(1:nc)[order(sapply(1:nc,function(cl){median(krt[rcmemb==cl],na.rm=TRUE)}))]
  if(intermeasure.flush)                #"Measurement" clusters (1st is baseline)
    mcl<-c(1,2*(1:length(conc)))
  else
    mcl<-1:(1+length(conc))
  conc<-c(0,conc)                       #Baseline concentration
  #Collect various statistics on the clustered points
  if(verbose)
    cat("Collecting summary statistics...\n")
  nam<-c(paste("Conc (",conc.units,")",sep=""), paste("LS TMin (",time.units,")",sep=""), paste("LS TMax (",time.units,")",sep=""), 
    sapply(scat.ang,function(a){
      c(
        paste("N LS ",a,"deg",sep=""),
        paste("Med LS ",a,"deg",sep=""),
        paste("Med CI 2.5% ",a,"deg",sep=""),
        paste("Med CI 97.5% ",a,"deg",sep="")
      )
    }),
    paste("dRITMin (",time.units,")",sep=""), paste("dRITMax (",time.units,")",sep=""),"N dRI","Med dRI","Med CI 2.5%","Med CI 97.5%")
  #print(scord)
  #print(scmemb)
  tab<-t(sapply(1:length(mcl),function(i){
    #Tabulate the scattering data
    sel<-scmemb==scord[mcl[i]]
    if(any(sel)){
      tr<-c(conc[i],min(kst[sel],na.rm=TRUE),max(kst[sel],na.rm=TRUE))
      for(j in 1:nang){
        v<-ksval[sel,j]
        tr<-c(tr,sum(!is.na(v)))
        v<-v[!is.na(v)]
        lv<-length(v)
        if(lv>1){
          tr<-c(tr,median(v),sort(v)[pmax(1,qbinom(c(0.025,0.975),lv,0.5))])
        }else{
           warning("Cluster associated with scattering measurement",i,"had one element - you should alter the thinning level.\n")
          tr<-c(tr,rep(median(v),3))
        }
      }
    }else{
      warning("Cluster associated with scattering measurement",i,"was empty - you should alter the thinning level.\n")
      tr<-c(conc[i],NA,NA,rep(c(0,NA,NA,NA),nang))
    }
    #Tabulate the refraction data
    sel<-rcmemb==rcord[mcl[i]]
    if(any(sel)){
      tr<-c(tr,min(krt[sel],na.rm=TRUE),max(krt[sel],na.rm=TRUE))
      v<-krval[sel]
      tr<-c(tr,sum(!is.na(v)))
      v<-v[!is.na(v)]
      lv<-length(v)
      if(lv>1){
        tr<-c(tr,median(v),sort(v)[pmax(1,qbinom(c(0.025,0.975),lv,0.5))])
      }else{
        warning("Cluster associated with refraction measurement",i,"had one element - you should alter the thinning level.\n")
        tr<-c(tr,rep(median(v),3))
      }
    }else{
      warning("Cluster associated with refraction measurement",i,"was empty - you should alter the thinning level.\n")
      tr<-c(tr,NA,NA,0,NA,NA,NA)
    }
    #cat("i=",i,"mcl[i]=",mcl[i],"\n")
    #print(tr)
    tr
  }))
  #print(tab)
  #print(dim(tab))
  #print(length(nam))
  #print(mcl)
  colnames(tab)<-nam
  rownames(tab)<-(1:NROW(tab))-1
  #Correct the concentrations using the refractive index data
  corrmod<-lm(tab[,1]~tab[,7+nang*4])  #Regress concentration on dRI
  conccorr<-predict(corrmod)           #dRI corrected concentrations
  conccorr<-conccorr-conccorr[1]       #Correct the corrections (at 0)
  tab<-cbind(tab,conccorr)
  colnames(tab)[NCOL(tab)]<-paste("Corr.Conc (",conc.units,")",sep="")
  tab<-tab[,c(1,NCOL(tab),2:(NCOL(tab)-1))]
  #Calculate parameters for regressions of scattering on concentration
  tab2<-t(sapply(6+((1:nang)-1)*4,function(i){
    c(coef(lm(tab[,i]~tab[,1])),coef(lm(tab[,i]~tab[,2])))
  }))
  tab2<-rbind(tab2, c(coef(lm(tab[,NCOL(tab)-2]~tab[,1])), coef(lm(tab[,NCOL(tab)-2]~tab[,2]))))
  rownames(tab2)<-c(paste("LS",scat.ang,"deg"),"dRI")
  colnames(tab2)<-c("Intercept",paste("Conc (",conc.units,")",sep=""), "Intercept",paste("Conc Corr. (",conc.units,")",sep=""))
  tab3<-t(sapply(6+((1:nang)-1)*4,function(i){  #1/R = f/M 1/c + 2fA2
    c(coef(lm(I(1/tab[-1,i])~I(1/tab[-1,1]))), coef(lm(I(1/tab[-1,i])~I(1/tab[-1,2]))), coef(lmrob(I(1/tab[-1,i])~I(1/tab[-1,2]),setting="KS2011")))
  }))
  v<-vector()
  for(i in 1:nang)
    v<-c(v,tab[-1,6+(i-1)*4])
  tab3<-rbind(tab3,c(coef(lm(I(1/v)~rep(1/tab[-1,1],nang))), coef(lm(I(1/v)~rep(1/tab[-1,2],nang))), coef(lmrob(I(1/v)~rep(1/tab[-1,2],nang),setting="KS2011"))))
  rownames(tab3)<-c(paste("LS",scat.ang,"deg"),"Combined")
  colnames(tab3)<-c("hat(2*f*A2) (Raw)","hat(f/M) (Raw)", "hat(2*f*A2) (Corr)","hat(f/M) (Corr)", "hat(2*f*A2) (RR, Corr)","hat(f/M) (RR, Corr)")
  #Organize the output and return
  if(verbose)
    cat("Returning data object\n")
  allscmemb<-rep(NA,length(stimes))
  allscmemb[keeps]<-scmemb
  allrcmemb<-rep(NA,length(rtimes))
  allrcmemb[keepr]<-rcmemb
  out<-list(results=tab, ols.coefs=tab2, param.est=tab3, include.LS=keeps, include.RI=keepr, clust.LS=allscmemb, clust.RI=allrcmemb, time.LS=stimes, time.RI=rtimes, vals.LS=sa, vals.RI=dri, trim.tol=trim.tol, instability.LS=dsdt, instability.RI=drdt, conc.corrected=conccorr, conc.correction.model=corrmod, conc.units=conc.units, time.units=time.units, scat.ang=scat.ang, intermeasure.flush=intermeasure.flush, measurement.clusters=mcl)
  class(out)<-"photoFinish"
  out
}


#A print method for photoFinish objects.  Use is self-explanatory
print.photoFinish<-function(x, ...){
  cat("\nQuality-trimmed Light Scattering Data, by Concentration\n\n")
  print(x$results,...)
  cat("\nOLS Coefficients for Median Measurements by Concentration\n\n")
  print(x$ols.coefs)
  cat("\nLS Parameter Estimates\n\n")
  print(x$param.est)
  cat("\n")
}


#A plot method for photoFinish objects.
#  Arguments:
#    x - a photoFinish object.
#    type - the type of plot to produce:
#      "results.LS" - Summary plot for light scattering results (medians by
#        concentration for each angle, w/95% CIs and an OLS line)
#      "results.RI" - Summary plot for refractive index results (medians by
#        concentration, w/95% CIs and an OLS line)
#      "params.LS" - Summary plot for conventional (?) estimation of LS parameters
#        by relating inverse scattering intensity to concentration
#      "time.LS" - Detailed display of trimmed light scatting data by time (all angles
#        overlaid).  Concentration regimes are indicated by colored clusters
#        (with lines indicating medians); the first and last clusters are taken
#        to be baselines, no medians are shown.
#      "time.RI" - Detailed display of trimmed refractive index data by time.  
#        Concentration regimes are indicated by colored clusters (with lines indicating 
#        medians); the first and last clusters are taken to be baselines, no medians 
#        are shown.
#      "raw.LS" - As "time.LS," but with untrimmed values included (black points).
#      "raw.RI" - As "time.RI," but with untrimmed values included (black points).
#    ... - additional arguments to plot()
#
plot.photoFinish<-function(x, type=c("results.LS","results.RI", "params.LS", "time.LS", "time.RI", "raw.LS", "raw.RI"), corrected.concentration=FALSE, fit.type=c("ols","robust","l1"), ...){
  require(robustbase)
  #Extract some generally useful stuff from the data object
  tab<-x$results
  conc<-tab[,1]
  concc<-tab[,2]
  nang<-length(x$scat.ang)
  ls<-tab[,6+((1:nang)-1)*4]
  #Produce the selected visualization
  if(match.arg(type)=="results.LS"){
    lslb<-tab[,7+((1:nang)-1)*4]
    lsub<-tab[,8+((1:nang)-1)*4]
    if(corrected.concentration){
      xv<-concc
      xlab<-paste("RI Corrected Concentration (",x$conc.units,")",sep="")
    }else{
      xv<-conc
      xlab<-paste("Concentration (",x$conc.units,")",sep="")
    }
    yv<-ls
    fit<-apply(yv,2,function(y,ft){
      switch(ft,
        ols=lm(y~xv),
        robust=lmrob(y~xv,setting="KS2011"),
        l1=lmrob.lar(xv,y)
      )
    },ft=match.arg(fit.type))
    pal<-hsv(((1:nang)-1)/(nang-1)*0.65)
    plot(rep(xv,nang),as.numeric(yv),xlab=xlab, ylab="Scattering Intensity", col=rep(pal,each=length(xv)), xlim=c(0,max(xv)), ylim=c(0,max(yv)), ...)
    segments(rep(xv,nang),as.numeric(lslb),rep(xv,nang), as.numeric(lsub),col=rep(pal,each=length(xv)))
    if(match.arg(fit.type)=="l1"){
      for(i in 1:nang)
        abline(a=coef(fit[[i]])[1],b=coef(fit[[i]])[2],lty=2,col=pal[i])
    }else{
      for(i in 1:nang)
        abline(fit[[i]],lty=2,col=pal[i])
    }
    abline(h=0,v=0,lty=3,col=rgb(0,0,0,0.5))
    legend("bottomright",fill=pal,legend=paste(x$scat.ang,"Degrees"),bty="n")
  }else if(match.arg(type)=="results.RI"){
    if(corrected.concentration){
      xv<-concc
      xlab<-paste("RI Corrected Concentration (",x$conc.units,")",sep="")
    }else{
      xv<-conc
      xlab<-paste("Concentration (",x$conc.units,")",sep="")
    }
    yv<-tab[,NCOL(tab)-2]
    fit<-switch(match.arg(fit.type),
      ols=lm(yv~xv),
      robust=lmrob(yv~xv,setting="KS2011"),
      l1=lmrob.lar(xv,yv)
    )
    plot(xv,yv,xlab=xlab, ylab="dRI", col=2, xlim=c(0,max(xv)), ylim=c(0,max(yv)), ...)
    segments(xv,tab[,NCOL(tab)-1],xv,tab[,NCOL(tab)],col=2)
    if(match.arg(fit.type)=="l1")
      abline(a=coef(fit)[1],b=coef(fit)[2],lty=2)
    else
      abline(fit,lty=2)
    abline(h=0,v=0,lty=3,col=rgb(0,0,0,0.5))
  }else if(match.arg(type)=="params.LS"){
    #This is based on the model 1/delta R(Q) = c f/M + A2 f/M
    pal<-hsv(((1:nang)-1)/(nang-1)*0.65)
    if(corrected.concentration){
      xv<-1/concc[-1]
      xlab<-paste("RI Corrected Inverse Concentration (1/",x$conc.units,")",sep="")
    }else{
      xv<-1/conc[-1]
      xlab<-paste("Inverse Concentration (1/",x$conc.units,")",sep="")
    }
    yv<-1/sweep(ls,2,ls[1,],"/")[-1,]  #Inverse delta R(Q)
    fit<-apply(yv,2,function(y,ft){
      switch(ft,
        ols=lm(y~xv),
        robust=lmrob(y~xv,setting="KS2011"),
        l1=lmrob.lar(cbind(rep(1,length(xv)),xv),y)
      )
    },ft=match.arg(fit.type))
    pal<-hsv(((1:nang)-1)/(nang-1)*0.65)
    plot(rep(xv,nang),as.numeric(yv),xlab=xlab, ylab="Inverse Scattering Intensity", col=rep(pal,each=length(xv)), xlim=c(0,max(xv)), ylim=range(c(as.numeric(yv),coef(fit)[1])), ...)
    if(match.arg(fit.type)!="l1"){
      for(i in 1:nang)
        abline(fit[[i]],lty=2,col=pal[i])
    }else{
      for(i in 1:nang)
        abline(a=coef(fit[[i]])[1],b=coef(fit[[i]])[2],lty=2,col=pal[i])
    }
    abline(h=0,v=0,lty=3,col=rgb(0,0,0,0.5))
    legend("bottomright",fill=pal,legend=paste(x$scat.ang,"Degrees"),bty="n")
  }else if(match.arg(type)=="time.LS"){
    sel<-x$include.LS
    cl<-unique(x$clust.LS[sel])
    pal<-hsv((cl-min(cl))/max(cl)*0.6,alpha=0.5)
    cols<-pal[match(x$clust.LS[sel],cl)]
    plot(x$time.LS[sel],x$vals.LS[sel,1],ylim=range(x$vals.LS[sel,]),type="p", xlab=paste("Time (",x$time.units,")",sep=""), ylab="Scattering Intensity", col=cols, cex=0.15, ...)
    for(i in 2:NCOL(x$vals.LS))
      points(x$time.LS[sel],x$vals.LS[sel,i],col=cols,cex=0.15)
    for(i in 1:nang)
      abline(h=ls[,i],col=pal[x$measurement.clusters],lwd=0.5,lty=3)
    if(corrected.concentration){
      xv<-concc
    }else{
      xv<-conc
    }
    legend("topleft",legend=(paste(signif(xv,2),x$conc.units)),fill=pal[x$measurement.clusters],bty="n")
  }else if(match.arg(type)=="time.RI"){
    sel<-x$include.RI
    cl<-unique(x$clust.RI[sel])
    pal<-hsv((cl-min(cl))/max(cl)*0.6,alpha=1)
    cols<-pal[match(x$clust.RI[sel],cl)]
    plot(x$time.RI[sel],x$vals.RI[sel],ylim=range(x$vals.RI[sel]),type="p", xlab=paste("Time (",x$time.units,")",sep=""), ylab="dRI", col=cols, cex=0.35, ...)
    abline(h=tab[,NCOL(tab)-2],col=pal[x$measurement.clusters],lwd=1,lty=3)
    if(corrected.concentration){
      xv<-concc
    }else{
      xv<-conc
    }
    legend("topleft",legend=(paste(signif(xv,2),x$conc.units)),fill=pal[x$measurement.clusters],bty="n")
  }else if(match.arg(type)=="raw.LS"){
    sel<-x$include.LS
    cl<-unique(x$clust.LS[sel])
    pal<-hsv((cl-min(cl))/max(cl)*0.6,alpha=0.5)
    cols<-pal[match(x$clust.LS,cl)]
    cols[is.na(cols)]<-rgb(0,0,0,0.75)
    plot(x$time.LS,x$vals.LS[,1],ylim=range(x$vals.LS,na.rm=TRUE),type="p", xlab=paste("Time (",x$time.units,")",sep=""), ylab="Scattering Intensity", col=cols, cex=0.25, ...)
    for(i in (1:NCOL(x$vals.LS))[-1])
      points(x$time.LS,x$vals.LS[,i],col=cols,cex=0.25)
    for(i in 1:nang)
      abline(h=ls[,i],col=pal[x$measurement.clusters],lwd=0.5,lty=3)
    if(corrected.concentration){
      xv<-concc
    }else{
      xv<-conc
    }
    legend("topleft",legend=(paste(signif(xv,2),x$conc.units)),fill=pal[x$measurement.clusters],bty="n")
  }else if(match.arg(type)=="raw.RI"){
    sel<-x$include.RI
    cl<-unique(x$clust.RI[sel])
    pal<-hsv((cl-min(cl))/max(cl)*0.6,alpha=0.75)
    cols<-pal[match(x$clust.RI,cl)]
    cols[is.na(cols)]<-rgb(0,0,0,0.75)
    plot(x$time.RI,x$vals.RI,ylim=range(x$vals.RI,na.rm=TRUE),type="p", xlab=paste("Time (",x$time.units,")",sep=""), ylab="dRI", col=cols, cex=0.5, ...)
    abline(h=tab[,NCOL(tab)-2],col=pal[x$measurement.clusters],lwd=1,lty=3)
    if(corrected.concentration){
      xv<-concc
    }else{
      xv<-conc
    }
    legend("topleft",legend=(paste(signif(xv,2),x$conc.units)),fill=pal[x$measurement.clusters],bty="n")
  }
}


#A slow but effective function to perform an agglomerative clustering of the points
#in vals that is constrained to a single-dimensional merge structure.  The cluster
#merge penalty is by a double median absolute difference (first, interpoint distances
#are taken to be the median of their absolute differences on each vals (column)
#dimension; and second, the penalty for merging two clusters is the median distance
#over all pairs of points in the respective clusters).  Merges are only possible
#between clusters that are adjacent in the underlying (row) dimension of vals,
#resulting in a solution that is guaranteed to be (1) contiguous, and (2) optimal
#with respect to the median merge penalty (given the number of clusters and
#contiguity constraint).
#
#  Arguments:
#    vals - matrix or vector of point values; rows must be in the same order as the
#           underlying dimension of interest.
#    nclust - number of clusters to return.
#
#  Return value:
#    A vector containing the cluster membership for each point (as an integer).
#
constrained.clust<-function(vals,nclust){
  if(length(dim(vals))<2)               #Ensure that our values are in matrix form
    vals<-matrix(vals,ncol=1)
  n<-NROW(vals)                         #Number of data vals
  cmemb<-1:n                            #Initial cluster memberships
  cl<-1:n                               #Initial cluster list
  nc<-n                                 #Initial cluster count
  mergecosts<-rep(Inf,nrow=n)           #Initial merge cost vector (i->i+1)
  for(i in 1:(nc-1)){
    diff<-median(abs(vals[i+1,]-vals[i,]),na.rm=TRUE)
    if(is.na(diff))
      diff<--Inf    #Merge NAs with any neighbor - they don't affect anything
    mergecosts[i]<-diff
  }
  while(nc>nclust){                     #Perform the clustering
    #Find the current best merge
    bm<-which(mergecosts==min(mergecosts))[1]
    #Merge cl[bm] with cl[bm+1]; first, relabel members of cl[bm+1] as cl[bm]
    cmemb[cmemb==cl[bm+1]]<-cl[bm]
    #Strike cl[bm+1] from our memory - we speak of it no more
    cl<-cl[-(bm+1)] 
    mergecosts<-mergecosts[-(bm+1)]
    nc<-nc-1
    #mergecosts[bm] used to be cl[bm] versus the unspeakable cluster; we now need to
    #update it to be equal to the merge cost for cl[bm] vs the new cl[bm+1].  We also
    #need to update the merge cost for cl[bm-1] vs cl[bm], if applicable.
    sel<-which(cmemb==(cl[bm]))         #Identify the members of cl[bm]
    if(bm<nc){                          #Find the cl[bm], cl[bm+1] merge cost
      selc<-which(cmemb==(cl[bm+1]))
      dm<-matrix(nrow=length(sel),ncol=length(selc))
      for(i in 1:length(sel))
        for(j in 1:length(selc))
          dm[i,j]<-median(abs(vals[sel[i]]-vals[selc[j]]),na.rm=TRUE)
      mergecosts[bm]<-median(dm,na.rm=TRUE)
      if(is.na(mergecosts[bm]))
        mergecosts[bm]<--Inf            #If all missing, can merge trivially
    }else
      mergecosts[bm]<-Inf               #There's no one above us!
    if(bm>1){                           #Find the cl[bm-1], cl[bm] merge cost
      selc<-which(cmemb==(cl[bm-1]))
      dm<-matrix(nrow=length(sel),ncol=length(selc))
      for(i in 1:length(sel))
        for(j in 1:length(selc))
          dm[i,j]<-median(abs(vals[sel[i]]-vals[selc[j]]),na.rm=TRUE)
      mergecosts[bm-1]<-median(dm,na.rm=TRUE)
      if(is.na(mergecosts[bm-1]))
        mergecosts[bm-1]<--Inf          #If all missing, can merge trivially
    }
  }
  #We are finished!  Return the cluster assignments
  assig<-rep(0,n)                       #Relabel for freshness!
  for(i in 1:nc)
    assig[cmemb==(cl[i])]<-i
  assig
}


#
virialEst<-function(con.meas, con.design, scat.meas, ri.meas, monomer.mass, scat.ang, 
                    n0.meas=NULL, wl=658e-7, draws=500, burn=1000, thin=250, 
                    prior.A2.sd=1, 
                    prior.agg.lambda=0.5, prior.con.lsd=log(1.05), 
                    prior.con.lsd.wgt=1, prior.con.meas.lsd=log(1.05), 
                    prior.con.meas.lsd.wgt=1, prior.scat.sd=1e-6, 
                    prior.scat.sd.wgt=1, prior.dndc.mean=0.185, 
                    prior.dndc.sd=1, prior.n0.mean=1, prior.n0.lsd=log(1.1), 
                    prior.n0.meas.sd=0.01, prior.ri.sd=1e-7, 
                    prior.ri.sd.wgt=1, max.oligomer.size=5, A2.prop.w0=1e-4, con.prop.w0=1e-4, dndc.prop.w0=1e-4, 
                    n0.prop.w0=1e-4, ri.s2.prop.w0=1e-4, 
                    param.seeds=list(con=con.meas, sigma=c(prior.con.lsd,prior.con.meas.lsd,prior.scat.sd,prior.ri.sd), 
                                     M=1, A2=0, dndc=NULL, n0=prior.n0.mean), verbose=TRUE){
  require(sna)
  require(MASS)
  #Define a semi-general function to draw a single slice sample in one dimension.
  #x0 here is the initial point, f is the log (unnormalized) target density, and
  #w0 is the initial test interval width.  The remaining arguments are passed to f,
  #which must take the candidate point as its first argument.
  drawslice<-function(x0,f,w0,...){
    #Get the initial log "height" at x0, and draw the target height
    ilh<-f(x0,...)
    th<-ilh+log(runif(1))  #log(u(0,1)*f(x0))
    #cat("Target height:",th,"\n")
    #Check to see if (x0+/-w0,exp(ilh)) are above f, and if not then expand
    w<-c(-w0,w0)
    #cat("w",w,"\n")
    while(f(x0+w[1],...)>th)  #Left-hand limit
      w[1]<-w[1]*2
    while(f(x0+w[2],...)>th)  #Right-hand limit
      w[2]<-w[2]*2
    #OK, we now have the x range of f<=th bracketed by w.  Draw x by rejection.
    x<-runif(1,x0+w[1],x0+w[2])
    while(f(x,...)<th){               #If we rejected, tighten the interval
      if(x<x0)
        w[1]<-x-x0
      if(x>x0)
        w[2]<-x-x0
      x<-runif(1,x0+w[1],x0+w[2])
    }
    #Done!  Now return x.
    x
  }
  #Define a truncated normal pdf
  dtnorm<-function(x,mean=0,sd=1,range=c(-Inf,Inf),log=FALSE){
    if(all(range==c(-Inf,Inf)))
      return(dnorm(x,mean=mean,sd=sd,log=log))
    else{
      lp<-dnorm(x,mean=mean,sd=sd,log=TRUE)
      lp<-lp-log(pnorm((range[2]-mean)/sd)-pnorm((range[1]-mean)/sd))
      lp[(x<range[1])|(x>range[2])]<--Inf
    }
    if(log)
      lp
    else
      exp(lp)    
  }
  #Define the conditional unnormalized log-posterior for A2, suitable for use with
  #drawslice.
  culp.A2<-function(x,M,con,dndc,n0,s2.ls,...){
    rQ.x<-rQ(M,con[has.scat],x,dndc,n0)
    dnorm(x,sd=prior.A2.sd,log=TRUE)+sum(dnorm(scat.meas[has.scat,], mean=rQ.x, sd=s2.ls^0.5, log=TRUE))
  }
  #Define the conditional unnormalized log-posterior for concentration elements,
  #suitable for use with drawslice.
  culp.con<-function(x,M,A2,dndc,n0,ind,s2.ls,s2.ri,s2.con,s2.con.meas,...){
    if(x<=0)
      return(-Inf)
    lp<-dlnorm(x,meanlog=log(con.design[ind]),sdlog=s2.con^0.5,log=TRUE) + dlnorm(con.meas[ind],meanlog=log(x),sd=s2.con.meas^0.5, log=TRUE)
    if(has.scat[ind])
      lp<-lp+sum(dnorm(scat.meas[ind,], mean=rQ(M,x,A2,dndc,n0), sd=s2.ls^0.5, log=TRUE))
    if(has.ri[ind])
      lp<-lp+dtnorm(ri.meas[ind],mean=x*dndc,sd=s2.ri^0.5,log=TRUE,range=c(0,Inf))
    lp
  }
  #Define the conditional unnormalized log-posterior for the refractive index
  #increment, suitable for use with drawslice.
  culp.dndc<-function(x,M,con,A2,n0,s2.ls,s2.ri,...){
    if(x<=0)
      return(-Inf)
    lp<-dtnorm(x,mean=prior.dndc.mean,sd=prior.dndc.sd,log=TRUE,range=c(0,Inf))
    rQ.x<-rQ(M,con[has.scat],A2,x,n0)
    lp<-lp+sum(dnorm(scat.meas[has.scat,],mean=rQ.x,sd=s2.ls^0.5,log=TRUE))
    lp<-lp+sum(dtnorm(ri.meas[has.ri],mean=con[has.ri]*x,sd=s2.ri^0.5,log=TRUE, range=c(0,Inf)))
    lp
  }
  #Define the conditional unnormalized log-posterior for the solvent refractive index,
  #suitable for use with drawslice.
  culp.n0<-function(x,M,con,A2,dndc,s2.ls,...){
    if(x<=0)
      return(-Inf)
    lp<-dlnorm(x,mean=log(prior.n0.mean),sd=prior.n0.lsd,log=TRUE)
    if(!is.null(n0.meas))
      lp<-lp+sum(dtnorm(n0.meas,mean=x,sd=prior.n0.meas.sd,log=TRUE,range=c(1,Inf)))
    rQ.x<-rQ(M,con[has.scat],A2,dndc,x)
    lp<-lp+sum(dnorm(scat.meas[has.scat,],mean=rQ.x,sd=s2.ls^0.5,log=TRUE))
    lp
  }
  #Define the conditional unnormalized log-posterior for refractive index measurement
  #variance, suitable for use with drawslice.
  culp.s2ri<-function(x,con,dndc,...){
    if((x<=0)|(x>1e300))
      return(-Inf)
    lp<-dgamma(1/x,prior.ri.sd.wgt/2,prior.ri.sd.wgt*prior.ri.sd/2,log=TRUE)
    lp<-lp+sum(dtnorm(ri.meas[has.ri],mean=con[has.ri]*dndc,sd=x^0.5,log=TRUE, range=c(0,Inf)))
    lp
  }
  #Create variables and whatnot for the simulation
  n<-length(con.meas)
  if(!is.matrix(scat.meas))
    scat.meas<-matrix(scat.meas,ncol=1)
  has.scat<-!is.na(scat.meas[,1])
  nls<-sum(has.scat)
  if(nls==0)
    stop("No light scattering measurements provided.\n")
  if((length(con.design)!=n)||(NROW(scat.meas)!=n))
    stop("Number of concentration measurements must match design values and scattering measurement points.\n") # n: # of concentration measurements
    # con.design = design values; scat.meas = scattering measurement points, now we have more concentration deisgn values
  has.ri<-!is.na(ri.meas)
  nr<-sum(has.ri)
  if(nr==0)  #Make sure we have some RI measurements
    stop("No RI measurements provided.\n")
  if(length(ri.meas)!=n)
    stop("Number of RI measurement entries doesn't match number of concentrations.\n")
  nang<-length(scat.ang)
  if(NCOL(scat.meas)!=nang)
    stop("Number of scattering angles in scat.meas does not match number of angles provided.\n")
  A2<-rep(0,draws)               #Second virial coefficient
  M<-rep(0,draws)                #Number of monomers per particle
  con<-matrix(0,draws,n)         #True concentrations
  s2<-matrix(0,draws,4)          #Variance parameters (con prior, con obs, scat obs, RI)
  dndc<-rep(0,draws)             #Refractive index increment
  n0<-rep(0,draws)               #Baseline solvent RI
  #Initialize the variables
  f<-(4*pi^2)/(wl^4*6.022e23)                     #Const part of Rayleigh "constant"
  rQ<-function(M,con,A2,dndc,n0){                 #Predicted adj scattering meas
    k<-f*(n0*dndc)^2                              #k=f* Solvent RI^2 * dndc^2
    dr<-k*con*M*monomer.mass*(1-2*A2*con*M*monomer.mass)  #Eq 8 from Wyatt, 1993
    t(sapply(dr,rep,nang))                        #Repeat for each angle
  }
  if(!is.null(param.seeds$con))
    con[1,]<-param.seeds$con
  else
    con[1,]<-con.meas
  if(!is.null(param.seeds$sigma))
    s2[1,]<-param.seeds$sigma^2
  else
    s2[1,]<-c(prior.con.lsd^2,prior.con.meas.lsd^2,prior.scat.sd^2,prior.ri.sd^2)
  if(!is.null(param.seeds$M))
    M[1]<-param.seeds$M
  else
    M[1]<-1
  if(!is.null(param.seeds$A2))
    A2[1]<-param.seeds$A2
  else
    A2[1]<-0
  if(!is.null(param.seeds$dndc))
    dndc[1]<-param.seeds$dndc
  else
    dndc[1]<-coef(lm(I(ri.meas[has.ri])~I(con.meas[has.ri])-1))
  if(!is.null(param.seeds$n0))
    n0[1]<-param.seeds$n0
  else{
    if(!is.na(n0.meas))
      n0[1]<-n0.meas
    else
      n0[1]<-prior.n0.mean
  }  
  #Run the MCMC
  i<-1
  bc<-1
  tc<-1
  while(i<=draws){
    if(verbose){
      if((bc<burn)&&(bc%%100==0))
        cat("Burn-in iteration",bc,"of",burn,"\n")
      else if((bc>=burn)&&(tc%%thin==1)&&(i%%50==0))
        cat("Taking draw",i,"of",draws,"\n")
    }
    #Update each parameter 
    #A|M,c,n0,dndc,s2r via slice sampling
    A2[i]<-drawslice(A2[i], f=culp.A2, w0=A2.prop.w0, M=M[i], con=con[i,], dndc=dndc[i], n0=n0[i], s2.ls=s2[i,3])
    #M-1|A,c,n0,dndc,s2r via inverse pmf
    Mlup<-dpois((1:max.oligomer.size)-1,prior.agg.lambda,log=TRUE)
    for(j in 1:max.oligomer.size)
      Mlup[j]<-Mlup[j]+sum(dnorm(scat.meas[has.scat,], mean=rQ(j,con[i,has.scat],A2[i],dndc[i],n0[i]), sd=s2[i,3]^0.5, log=TRUE))
    Mlup<-Mlup-logSum(Mlup)
    M[i]<-min(which(cumsum(exp(Mlup))>=runif(1)))
    #Draw the concentration values, some of which have scattering and some RI vals,
    #using slice sampling
    for(j in 1:n)
      con[i,j]<-drawslice(con[i,j], f=culp.con, w0=con.prop.w0, M=M[i], A2=A2[i], dndc=dndc[i], n0=n0[i], ind=j, s2.ls=s2[i,3], s2.ri=s2[i,4], s2.con=s2[i,1], s2.con.meas=s2[i,2])
    #s2cd|c ~ IG((ncd0+n)/2,s2cd0+sum_i (ci-cdi)^2/2)
    ssq<-sum((log(con[i,])-log(con.design))^2)
    s2[i,1]<-1/rgamma(1,(prior.con.lsd.wgt+n)/2, prior.con.lsd.wgt*prior.con.lsd^2/2 + ssq/2)
    #s2c|c ~ IG((nc0+n)/2,s2c0+sum_i (coi-ci)^2/2)
    ssq<-sum((log(con[i,])-log(con.meas))^2)
    s2[i,2]<-1/rgamma(1,(prior.con.meas.lsd.wgt+n)/2, prior.con.meas.lsd.wgt*prior.con.meas.lsd^2/2 + ssq/2)
    #s2r|c,A,M ~ IG((nr0+n)/2,s2r0+sum_i (roi-1/(f/(mMci) + 2fA))^2/2)
    rQ.cur<-rQ(M[i],con[i,has.scat],A2[i],dndc[i],n0[i])
    s2[i,3]<-1/rgamma(1,(prior.scat.sd.wgt+nang*nls)/2, prior.scat.sd.wgt*prior.scat.sd^2/2 + sum((rQ.cur-scat.meas[has.scat,])^2)/2)
    #Draw dndc using slice sampling
    dndc[i]<-drawslice(dndc[i], f=culp.dndc, w0=dndc.prop.w0, M=M[i], con=con[i,], A2=A2[i], n0=n0[i], s2.ls=s2[i,3], s2.ri=s2[i,4])
    #Draw n0 using slice sampling
    n0[i]<-drawslice(n0[i], f=culp.n0, w0=n0.prop.w0, M=M[i], con=con[i,], A2=A2[i], dndc=dndc[i], s2.ls=s2[i,3])
    #Draw s2n using slice sampling
    #The code below is now obsolete; keeping it in case we go back to unconstrained
    #Gaussian errors:
    ##s2b|c,bo,d ~ IG(nb0/2+n/2,nb0*s2b0/2+sum((bo-x*d)^2))
    #s2[i,4]<-1/rgamma(1,prior.ri.sd.wgt/2+nr/2, prior.ri.sd.wgt*prior.ri.sd^2/2 + sum((ri.meas[has.ri] - con[i,has.ri]*dndc[i])^2)/2)
    s2[i,4]<-drawslice(s2[i,4], f=culp.s2ri, w0=ri.s2.prop.w0, con=con[i,], dndc=dndc[i])
    #Update counters
    if(bc<burn)
      bc<-bc+1
    else{
      if(tc%%thin==0){
        if(i<draws){
          A2[i+1]<-A2[i]
          M[i+1]<-M[i]
          con[i+1,]<-con[i,]
          dndc[i+1]<-dndc[i]
          n0[i+1]<-n0[i]
          s2[i+1,]<-s2[i,]
        }
        i<-i+1
      }
      tc<-tc+1
    }
  }
  #Organize and return the results
  out<-list(A2=A2, M=M*monomer.mass, aggregate.size=M, concentration=con, dndc=dndc, n0=n0, sigma2.con=s2[,1], 
            sigma2.con.meas=s2[,2], sigma2.scat=s2[,3], sigma2.RI=s2[,4], concentration.design=con.design, 
            concentration.measured=con.meas, RI.measured=ri.meas, scat.ang=scat.ang, wavelength=wl, 
            monomer.mass=monomer.mass, scattering.meas=scat.meas, param.seeds=param.seeds, draws=draws, thin=thin, burn=burn)
  class(out)<-"virialEst"
  out
}

# inside function: M how many particles we have
# 

#K(theta)c/(RA(theta)-R0(theta)) = 
#             1/W (1 + q(theta)^2 R^2/3 + O(q(theta)^4))(1 + 2WAc + O(c^2))
#         apx 1/W (1 + q(theta)^2 R^2/3)(1 + 2WAc)
#         apx 1/W + 2 A c
#  1/(RA(theta)-R0(theta)) apx 1/K(theta) 1/(Wc) + 2 A/K(theta)
#                          = f/(Wc) + 2 f A
#  RA(thea)-R0(theta) apx K(theta)/(1/(Wc) + 2A)
#K(theta) = 4 pi^2 n0^2 (dn/dc)^2/(6.022e23 lambda^4) = 1/f
#q(theta) = 4*pi*n0*sin(theta/2)/lambda

#ci ~ N(cdi,s2cd)
#s2cd ~ IG(ncd0/2,ncd0*s2cd0/2)
#coi | ci ~ N(ci,s2c)
#s2c ~ IG(nc0/2,nc0*s2c0/2)
#roi | ci,M,A,s2r ~ N(1/(/f(mMci) + 2fA),s2r)
#s2r ~ IG(nr0/2,nr0*s2r0/2)
#A ~ N(0,s2a)
#M-1 ~ Pois(l)
#A|M,c,s2r propto N(A|0,s2a) prod_i N(roi|1/(f/(mMci) + 2fA),s2r)
#   ~ N(A|)
#M-1|A,c,s2r propto Pois(M|l) prod_i N(roi|1/(f/(mMci) + 2fA),s2r)
#ci|d2cd,s2c,roi,M,A,s2r propto N(ci|cdi,s2cd) N(coi|ci,s2c) N(roi|1/(f/(mMci) + 2fA),s2r)
#s2cd|c ~ IG((ncd0+n)/2,ncd0*s2cd0/2+sum_i (ci-cdi)^2/2)
#s2c|c ~ IG((nc0+n)/2,nc0*s2c0/2+sum_i (coi-ci)^2/2)
#s2r|c,A,M ~ IG((nr0+n)/2,nr0*s2r0/2+sum_i (roi-1/(f/(mMci) + 2fA))^2/2)
#boi ~ N(d*ci+do, s2b)
#(d,do) ~ N((d0,do0),s2d0)
#s2b ~ IG(nd0/2,nd0*s2b0/2)
#x=(1,c), s0inv=diag(1/s2d0), sm=(s0inv+x^Tx)^-1
#d|bo,c,s2b ~ mvn(sm*(s0inv*(do0,d0) + x^T*ri),s2b*sm)
#s2b|c,bo,d ~ IG(nd0/2+n/2,nd0*s2b0/2+1/2(ri^T*ri + d0^T*s0inv*d0 + d^T*sm^-1*d))
#R ~ lnorm(mR0,s2r0)
virialEst.old<-function(con.meas, con.design, scat.meas, ri.meas, monomer.mass, scat.ang, wl=658e-7, draws=500, burn=1000, thin=250, prior.A2.sd=1, prior.A2.bias.sd=1, prior.agg.lambda=0.5, prior.con.lsd=log(1.05), prior.con.lsd.wgt=1, prior.con.meas.lsd=log(1.05), prior.con.meas.lsd.wgt=1, prior.scat.sd=1e-6, prior.scat.sd.wgt=1, prior.Rg.lmean=log(1.5e-7), prior.Rg.lsd=1, prior.dndc.mean=0.185, prior.dndc.sd=1, prior.n0.mean=1, prior.n0.lsd=log(1.05), prior.ri.sd=1e-7, prior.ri.sd.wgt=1, max.oligomer.size=5, A2.prop.w0=1e-4, lnRg.prop.w0=1e-4, con.prop.w0=1e-4, dndc.prop.w0=1e-4, n0.prop.w0=1e-4, verbose=TRUE){
  require(sna)
  require(MASS)
  #Define a semi-general function to draw a single slice sample in one dimension.
  #x0 here is the initial point, f is the log (unnormalized) target density, and
  #w0 is the initial test interval width.  The remaining arguments are passed to f,
  #which must take the candidate point as its first argument.
  drawslice<-function(x0,f,w0,...){
    #Get the initial log "height" at x0, and draw the target height
    ilh<-f(x0,...)
    th<-ilh+log(runif(1))  #log(u(0,1)*f(x0))
    #Check to see if (x0+/-w0,exp(ilh)) are above f, and if not then expand
    w<-c(-w0,w0)
    while(f(x0+w[1],...)>th)  #Left-hand limit
      w[1]<-w[1]*2
    while(f(x0+w[2],...)>th)  #Right-hand limit
      w[2]<-w[2]*2
    #OK, we now have the x range of f<=th bracketed by w.  Draw x by rejection.
    x<-runif(1,x0+w[1],x0+w[2])
    while(f(x,...)<th){               #If we rejected, tighten the interval
      if(x<x0)
        w[1]<-x-x0
      if(x>x0)
        w[2]<-x-x0
      x<-runif(1,x0+w[1],x0+w[2])
    }
    #Done!  Now return x.
    x
  }
  #Define the conditional unnormalized log-posterior for A2, suitable for use with
  #drawslice.
  culp.A2<-function(x,lnRg,M,con,dndc,n0,A2b,s2.ls,...){
    rQ.x<-rQ(lnRg,M,con[has.scat],x,dndc,n0,A2b)
    dnorm(x,sd=prior.A2.sd,log=TRUE)+sum(dnorm(scat.meas[has.scat,], mean=rQ.x, sd=s2.ls^0.5, log=TRUE))
  }
  #Define the conditional unnormalized log-posterior for lnRg, suitable for use with
  #drawslice.
  culp.lnRg<-function(x,M,con,A2,dndc,n0,A2b,s2.ls,...){
    rQ.x<-rQ(x,M,con[has.scat],A2,dndc,n0,A2b)
    dlnorm(exp(x),meanlog=prior.Rg.lmean,sdlog=prior.Rg.lsd,log=TRUE)+sum(dnorm(scat.meas[has.scat,], mean=rQ.x, sd=s2.ls^0.5, log=TRUE))
  }
  #Define the conditional unnormalized log-posterior for concentration elements,
  #suitable for use with drawslice.
  culp.con<-function(x,lnRg,M,A2,dndc,n0,A2b,ind,s2.ls,s2.ri,s2.con,s2.con.meas,...){
    if(x<=0)
      return(-Inf)
    lp<-dlnorm(x,meanlog=log(con.design[ind]),sdlog=s2.con^0.5,log=TRUE) + dlnorm(con.meas[ind],meanlog=log(x),sd=s2.con.meas^0.5, log=TRUE)
    if(has.scat[ind])
      lp<-lp+sum(dnorm(scat.meas[ind,], mean=rQ(lnRg,M,x,A2,dndc,n0,A2b), sd=s2.ls^0.5, log=TRUE))
    if(has.ri[ind])
      lp<-lp+dnorm(ri.meas[ind],mean=x*dndc,sd=s2.ri^0.5,log=TRUE)
    lp
  }
  #Define the conditional unnormalized log-posterior for the refractive index
  #increment, suitable for use with drawslice.
  culp.dndc<-function(x,lnRg,M,con,A2,n0,A2b,s2.ls,s2.ri,...){
    if(x<=0)
      return(-Inf)
    lp<-dnorm(x,mean=prior.dndc.mean,sd=prior.dndc.sd,log=TRUE)
    rQ.x<-rQ(lnRg,M,con[has.scat],A2,x,n0,A2b)
    lp<-lp+sum(dnorm(scat.meas[has.scat,],mean=rQ.x,sd=s2.ls^0.5,log=TRUE))
    lp<-lp+sum(dnorm(ri.meas[has.ri],mean=con[has.ri]*x,sd=s2.ri^0.5,log=TRUE))
    lp
  }
  #Define the conditional unnormalized log-posterior for the solvent refractive index,
  #suitable for use with drawslice.
  culp.n0<-function(x,lnRg,M,con,A2,dndc,A2b,s2.ls,...){
    if(x<=0)
      return(-Inf)
    lp<-dlnorm(x,mean=log(prior.n0.mean),sd=prior.n0.lsd,log=TRUE)
    rQ.x<-rQ(lnRg,M,con[has.scat],A2,dndc,x,A2b)
    lp<-lp+sum(dnorm(scat.meas[has.scat,],mean=rQ.x,sd=s2.ls^0.5,log=TRUE))
    lp
  }
  #Create variables and whatnot for the simulation
  n<-length(con.meas)
  if(!is.matrix(scat.meas))
    scat.meas<-matrix(scat.meas,ncol=1)
  has.scat<-!is.na(scat.meas[,1])
  nls<-sum(has.scat)
  if(nls==0)
    stop("No light scattering measurements provided.\n")
  if((length(con.design)!=n)||(NROW(scat.meas)!=n))
    stop("Number of concentration measurements must match design values and scattering measurement points.\n")
  has.ri<-!is.na(ri.meas)
  nr<-sum(has.ri)
  if(nr==0)  #Make sure we have some RI measurements
    stop("No RI measurements provided.\n")
  if(length(ri.meas)!=n)
    stop("Number of RI measurement entries doesn't match number of concentrations.\n")
  nang<-length(scat.ang)
  if(NCOL(scat.meas)!=nang)
    stop("Number of scattering angles in scat.meas does not match number of angles provided.\n")
  A2<-rep(0,draws)               #Second virial coefficient
  M<-rep(0,draws)                #Number of monomers per particle
  lnRg<-rep(0,draws)             #Log radius of gyration
  con<-matrix(0,draws,n)         #True concentrations
  s2<-matrix(0,draws,4)          #Variance parameters (con prior, con obs, scat obs, RI)
  dndc<-rep(0,draws)             #Refractive index increment
  n0<-rep(0,draws)               #Baseline solvent RI
  #Initialize the variables
  qon0<-4*pi*sin(scat.ang/360*pi)/wl              #q/n0=4*pi*sin(theta/2)/lambda
  f<-(4*pi^2)/(wl^4*6.022e23)                     #My version of const for Zimm model
  rQ<-function(lnRg,M,con,A2,dndc,n0){            #Predicted adj scattering meas
    fd<-f*(n0*dndc)^2                             #Solvent RI^2 * dndc^2
    a2lim<-outer(con,qon0*n0,function(x,y){fd/((1+y^2*exp(2*lnRg)/3)*(1/(monomer.mass*M*x)+2*A2))})
    a2lim+t(sapply(A2b/con,rep,nang))             #Add the low-conc bias effect
  }
  con[1,]<-con.meas
  s2[1,]<-c(prior.con.lsd^2,prior.con.meas.lsd^2,prior.scat.sd^2,prior.ri.sd^2)
  M[1]<-1
  A2[1]<-0
  A2b[1]<-0
  lnRg[1]<-prior.Rg.lmean
  dndc[1]<-coef(lm(I(ri.meas[has.ri])~I(con.meas[has.ri])-1))
  n0[1]<-prior.n0.mean
  #Run the MCMC
  i<-1
  bc<-1
  tc<-1
  while(i<=draws){
    if(verbose){
      if((bc<burn)&&(bc%%100==0))
        cat("Burn-in iteration",bc,"of",burn,"\n")
      else if((bc>=burn)&&(tc%%thin==1)&&(i%%50==0))
        cat("Taking draw",i,"of",draws,"\n")
    }
    #Update each parameter 
    #A|M,c,s2r propto N(A|0,s2a) prod_i N(roi|1/(f/(mMci) + 2fA),s2r)
    A2[i]<-drawslice(A2[i], f=culp.A2, w0=A2.prop.w0, lnRg=lnRg[i], M=M[i], con=con[i,], dndc=dndc[i], n0=n0[i], s2.ls=s2[i,3])
    #cand<-A2[i]+rnorm(1,sd=A2.prop.sd)
    #rQ.cur<-rQ(lnRg[i],M[i],con[i,],A2[i],d[i,])
    #rQ.cand<-rQ(lnRg[i],M[i],con[i,],cand,d[i,])
    #sel<-con[i,]>0
    #a2lup<-dnorm(A2[i],sd=prior.A2.sd,log=TRUE)+sum(dnorm(scat.meas[sel,],mean=rQ.cur[sel,], sd=sqrt(s2[i,3]), log=TRUE))
    #cand.a2lup<-dnorm(cand,sd=prior.A2.sd,log=TRUE)+sum(dnorm(scat.meas[sel,],mean=rQ.cand[sel,], sd=sqrt(s2[i,3]), log=TRUE))
    ##cat("A2=",A2[i],"cand=",cand,"a2lup=",a2lup,"cand.a2lup=",cand.a2lup,"diff=", cand.a2lup-a2lup,"errdiff",sum((scat.meas-rQ.cand)[sel,]^2)-sum((scat.meas-rQ.cur)[sel,]^2),"\n")
    #if(cand.a2lup-a2lup>log(runif(1))){
    #  #cat("\tAccepted - err=",sum((scat.meas-rQ.cand)[sel,]^2),"\n")
     # A2[i]<-cand
    #}
    #M-1|A,c,s2r propto Pois(M|l) prod_i N(roi|1/(f/(mMci) + 2fA),s2r)
    Mlup<-dpois((1:max.oligomer.size)-1,prior.agg.lambda,log=TRUE)
    for(j in 1:max.oligomer.size)
      Mlup[j]<-Mlup[j]+sum(dnorm(scat.meas[has.scat,], mean=rQ(lnRg[i],j,con[i,has.scat],A2[i],dndc[i],n0[i]), sd=s2[i,3]^0.5, log=TRUE))
    Mlup<-Mlup-logSum(Mlup)
    #cat("Mlup\n")
    #print(Mlup)
    #cat("Normalized M post\n")
    #print(cumsum(exp(Mlup)))
    M[i]<-min(which(cumsum(exp(Mlup))>=runif(1)))
    #Rg|A,c,M,s2r propto lnorm(Rg|mR0,s2R0) prod_i N(roi|1/(f/(mMci) + 2fA),s2r)
    lnRg[i]<-drawslice(lnRg[i], f=culp.lnRg, w0=lnRg.prop.w0, M=M[i], con=con[i,], A2=A2[i], dndc=dndc[i], n0=n0[i], s2.ls=s2[i,3])
    #cand<-lnRg[i]+rnorm(1,sd=Rg.prop.sd)
    #rQ.cur<-rQ(lnRg[i],M[i],con[i,],A2[i],d[i,])
    #rQ.cand<-rQ(cand,M[i],con[i,],A2[i],d[i,])
    #sel<-con[i,]>0
    ##cat("cur",lnRg[i],"\n")
    ##print(rQ.cur)
    ###cat("cand",cand,"\n")
    ##print(rQ.cand)
    #rglup<-dlnorm(exp(lnRg[i]),meanlog=prior.Rg.lmean,sdlog=prior.Rg.lsd,log=TRUE)+sum(dnorm(scat.meas[sel,],mean=rQ.cur[sel,], sd=sqrt(s2[i,3]), log=TRUE))
    #cand.rglup<-dlnorm(exp(cand),meanlog=prior.Rg.lmean,sdlog=prior.Rg.lsd,log=TRUE)+sum(dnorm(scat.meas[sel,],mean=rQ.cand[sel,], sd=sqrt(s2[i,3]), log=TRUE))
    #if(cand.rglup-rglup>log(runif(1))){
    #  lnRg[i]<-cand
    #}
    #Draw the concentration values, some of which have scattering and some RI vals
    for(j in 1:n)
      con[i,j]<-drawslice(con[i,j], f=culp.con, w0=con.prop.w0, lnRg=lnRg[i], M=M[i], A2=A2[i], dndc=dndc[i], n0=n0[i], ind=j, s2.ls=s2[i,3], s2.ri=s2[i,4], s2.con=s2[i,1], s2.con.meas=s2[i,2])
    #if(common.con.meas){  #Concentration values are linked across scattering, RI
    #  #ci|d2cd,s2c,roi,M,A,s2r propto N(ci|cdi,s2cd) N(coi|ci,s2c) 
    #  #  * N(roi|1/(f/(mMci) + 2fA),s2r) N(bi|d[1]+d[2]*ci,s2b)
    #  cand<-con[i,]*exp(rnorm(n,sd=con.prop.lsd))
    #  rQ.cur<-rQ(lnRg[i],M[i],con[i,],A2[i],d[i,])
    #  rQ.cand<-rQ(lnRg[i],M[i],cand,A2[i],d[i,])
    #  RI.cur<-d[i,1]+d[i,2]*con[i,]
    #  RI.cand<-d[i,1]+d[i,2]*cand
    #  sel<-con.design>0
    #  conlup<-c(0,dlnorm(con[i,sel],meanlog=log(con.design[sel]),sdlog=s2[i,1]^0.5,log=TRUE)) + c(0,dlnorm(con.meas[sel],meanlog=log(con[i,sel]),sd=s2[i,2]^0.5, log=TRUE)) + c(0,rowSums(dnorm(scat.meas[sel,],mean=rQ.cur[sel,],sd=s2[i,3]^0.5,log=TRUE))) + dnorm(ri.meas,mean=RI.cur,sd=s2[i,4]^0.5,log=TRUE)
    #  conlup.cand<-c(0,dlnorm(cand[sel],mean=log(con.design[sel]),sd=s2[i,1]^0.5,log=TRUE)) + c(0,dlnorm(con.meas[sel],mean=log(cand[sel]),sd=s2[i,2]^0.5,log=TRUE)) + c(0,rowSums(dnorm(scat.meas[sel,],mean=rQ.cand[sel,], sd=s2[i,3]^0.5,log=TRUE))) + dnorm(ri.meas,mean=RI.cand,sd=s2[i,4]^0.5,log=TRUE)
    #  #cat("con stuff:\n")
    #  #print(conlup)
    #  #print(conlup.cand)
    #  sel<-(conlup.cand-conlup)>log(runif(n))
    #  if(any(sel)){
    #    #cat("Accepted",sum(sel),"\n")
    #    con[i,sel]<-cand[sel]
    #  }
    #}else{                #Different concentration values for scattering, RI
    #  #First, the scattering concentrations
    #  #ci|d2cd,s2c,roi,M,A,s2r propto N(ci|cdi,s2cd) N(coi|ci,s2c) N(roi|1/(f/(mMci) + 2fA),s2r)
    #  cand<-con[i,]*exp(rnorm(n,sd=con.prop.lsd))
    #  rQ.cur<-rQ(lnRg[i],M[i],con[i,],A2[i],d[i,])
    #  rQ.cand<-rQ(lnRg[i],M[i],cand,A2[i],d[i,])
    #  sel<-con[i,]>0
    #  conlup<-ldcon(con[i,],mean=con.design,sd=s2[i,1]^0.5) + ldcon(con.meas,mean=con[i,],sd=s2[i,2]^0.5) + c(0,rowSums(dnorm(scat.meas[sel,],mean=rQ.cur[sel,],sd=s2[i,3]^0.5,log=TRUE)))
    #  conlup.cand<-ldcon(cand,mean=con.design,sd=s2[i,1]^0.5) + ldcon(con.meas,mean=cand,sd=s2[i,2]^0.5) + c(0,rowSums(dnorm(scat.meas[sel,],mean=rQ.cand[sel,], sd=s2[i,3]^0.5,log=TRUE)))
    #  #print(con[i,])
    #  #print(cand)
    #  #print(conlup)
    #  #print(conlup.cand)
    #  sel<-(conlup.cand-conlup)>log(runif(n))
    #  if(any(sel))
    #    con[i,sel]<-cand[sel]
    #  #Now, the RI concentrations
    #  #ci|d2cd,s2c,roi,M,A,s2r propto N(ci|cdi,s2cd) N(coi|ci,s2c) N(bi|d[1]+d[2]*ci),s2b)
    #  cand<-conri[i,]*exp(rnorm(nr,sd=con.prop.lsd))
    #  RI.cur<-d[i,1]+d[i,2]*conri[i,]
    #  RI.cand<-d[i,1]+d[i,2]*cand
    #  conlup<-ldcon(conri[i,],mean=con.design.ri,sd=s2[i,1]^0.5) + ldcon(con.meas.ri,mean=conri[i,],sd=s2[i,2]^0.5) + dnorm(ri.meas,mean=RI.cur,sd=s2[i,4]^0.5,log=TRUE)
    #  conlup.cand<-ldcon(cand,mean=con.design.ri,sd=s2[i,1]^0.5) + ldcon(con.meas.ri,mean=cand,sd=s2[i,2]^0.5) + dnorm(ri.meas,mean=RI.cand,sd=s2[i,4]^0.5,log=TRUE)
    #  sel<-(conlup.cand-conlup)>log(runif(nr))
    #  if(any(sel))
    #    conri[i,sel]<-cand[sel]
    #}
    #s2cd|c ~ IG((ncd0+n)/2,s2cd0+sum_i (ci-cdi)^2/2)
    ssq<-sum((log(con[i,])-log(con.design))^2)
    s2[i,1]<-1/rgamma(1,(prior.con.lsd.wgt+n)/2, prior.con.lsd.wgt*prior.con.lsd^2/2 + ssq/2)
    #if(common.con.meas){
    #  ssq<-(log(con[i,])-log(con.design))^2
    #  sel<-(con[i,]*con.design)>0
    #  ssq<-sum(ssq[sel])
    #  s2[i,1]<-1/rgamma(1,(prior.con.lsd.wgt+sum(sel))/2, prior.con.lsd.wgt*prior.con.lsd^2/2 + ssq/2)
    #}else{
    #  ssq<-c((log(con[i,])-log(con.design))^2,(log(conri[i,])-log(con.design.ri))^2)
    #  sel<-c((con[i,]*con.design)>0,(conri[i,]*con.design.ri)>0)
    #  ssq<-sum(ssq[sel])
    #  s2[i,1]<-1/rgamma(1,(prior.con.lsd.wgt+sum(sel))/2, prior.con.lsd.wgt*prior.con.lsd^2/2 + ssq/2)
    #}
    #s2c|c ~ IG((nc0+n)/2,s2c0+sum_i (coi-ci)^2/2)
    ssq<-sum((log(con[i,])-log(con.meas))^2)
    s2[i,2]<-1/rgamma(1,(prior.con.meas.lsd.wgt+n)/2, prior.con.meas.lsd.wgt*prior.con.meas.lsd^2/2 + ssq/2)
    #if(common.con.meas){
    #  ssq<-(log(con[i,])-log(con.meas))^2
    #  sel<-(con[i,]*con.meas)>0
    #  ssq<-sum(ssq[sel])
    #  s2[i,2]<-1/rgamma(1,(prior.con.meas.lsd.wgt+sum(sel))/2, prior.con.meas.lsd.wgt*prior.con.meas.lsd^2/2 + ssq/2)
    #}else{
    #  ssq<-c((log(con[i,])-log(con.meas))^2,(log(conri[i,])-log(con.meas.ri))^2)
    #  sel<-c((con[i,]*con.meas)>0,(conri[i,]*con.meas.ri)>0)
    #  ssq<-sum(ssq[sel])
    #  s2[i,2]<-1/rgamma(1,(prior.con.meas.lsd.wgt+sum(sel))/2, prior.con.meas.lsd.wgt*prior.con.meas.lsd^2/2 + ssq/2)
    #}
    #s2r|c,A,M ~ IG((nr0+n)/2,s2r0+sum_i (roi-1/(f/(mMci) + 2fA))^2/2)
    rQ.cur<-rQ(lnRg[i],M[i],con[i,has.scat],A2[i],dndc[i],n0[i])
    s2[i,3]<-1/rgamma(1,(prior.scat.sd.wgt+nang*nls)/2, prior.scat.sd.wgt*prior.scat.sd^2/2 + sum((rQ.cur-scat.meas[has.scat,])^2)/2)
    #cat("Scat: ssq/df=",sum((rQ.cur-scat.meas)[sel,]^2)/(nang*sum(sel)),"sigma2=",s2[i,3],"\n")
    #Draw dndc using slice sampling
    dndc[i]<-drawslice(dndc[i], f=culp.dndc, w0=dndc.prop.w0, lnRg=lnRg[i], M=M[i], con=con[i,], A2=A2[i], n0=n0[i], s2.ls=s2[i,3], s2.ri=s2[i,4])
    #if(common.con.meas)
    #  x<-cbind(rep(1,nr),con[i,])
    #else
    #  x<-cbind(rep(1,nr),conri[i,])
    #cand<-d[i,]+c(rnorm(1,sd=dndc.prop.sd),0)
    #rQ.cur<-rQ(lnRg[i],M[i],con[i,],A2[i],d[i,])
    #rQ.cand<-rQ(lnRg[i],M[i],con[i,],A2[i],cand)
    #RI.cur<-x%*%d[i,]
    #RI.cand<-x%*%cand
    #sel<-con.design>0
    #dlup<-sum(dnorm(d[i,],mean=prior.dndc.mean,sd=prior.dndc.sd,log=TRUE)) + sum(dnorm(scat.meas[sel,],mean=rQ.cur[sel,],sd=s2[i,3]^0.5,log=TRUE)) + sum(dnorm(ri.meas,mean=RI.cur,sd=s2[i,4]^0.5,log=TRUE))
    #cand.dlup<-sum(dnorm(cand,mean=prior.dndc.mean,sd=prior.dndc.sd,log=TRUE)) + sum(dnorm(scat.meas[sel,],mean=rQ.cand[sel,],sd=s2[i,3]^0.5,log=TRUE)) + sum(dnorm(ri.meas,mean=RI.cand,sd=s2[i,4]^0.5,log=TRUE))
    #if(cand.dlup-dlup>log(runif(1))){
    #  d[i,]<-cand
    #}
    #cand<-d[i,]+c(0,rnorm(1,sd=dndc.prop.sd))
    #rQ.cur<-rQ(lnRg[i],M[i],con[i,],A2[i],d[i,])
    #rQ.cand<-rQ(lnRg[i],M[i],con[i,],A2[i],cand)
    #RI.cur<-x%*%d[i,]
    #RI.cand<-x%*%cand
    #sel<-con.design>0
    #dlup<-sum(dnorm(d[i,],mean=prior.dndc.mean,sd=prior.dndc.sd,log=TRUE)) + sum(dnorm(scat.meas[sel,],mean=rQ.cur[sel,],sd=s2[i,3]^0.5,log=TRUE)) + sum(dnorm(ri.meas,mean=RI.cur,sd=s2[i,4]^0.5,log=TRUE))
    #cand.dlup<-sum(dnorm(cand,mean=prior.dndc.mean,sd=prior.dndc.sd,log=TRUE)) + sum(dnorm(scat.meas[sel,],mean=rQ.cand[sel,],sd=s2[i,3]^0.5,log=TRUE)) + sum(dnorm(ri.meas,mean=RI.cand,sd=s2[i,4]^0.5,log=TRUE))
    #if(cand.dlup-dlup>log(runif(1))){
    #  d[i,]<-cand
    #}
    #Draw dndc using slice sampling
    n0[i]<-drawslice(n0[i], f=culp.n0, w0=n0.prop.w0, lnRg=lnRg[i], M=M[i], con=con[i,], A2=A2[i], dndc=dndc[i], s2.ls=s2[i,3])
    #s2b|c,bo,d ~ IG(nb0/2+n/2,nb0*s2b0/2+sum((bo-x*d)^2))
    s2[i,4]<-1/rgamma(1,prior.ri.sd.wgt/2+nr/2, prior.ri.sd.wgt*prior.ri.sd^2/2 + sum((ri.meas[has.ri] - con[i,has.ri]*dndc[i])^2)/2)
    #cat("RI: ssq/df=",sum((ri.meas-x%*%d[i,])^2)/nr,"sigma2=",s2[i,3],"\n")
    #Update counters
    if(bc<burn)
      bc<-bc+1
    else{
      if(tc%%thin==0){
        if(i<draws){
          A2[i+1]<-A2[i]
          M[i+1]<-M[i]
          con[i+1,]<-con[i,]
          dndc[i+1]<-dndc[i]
          n0[i+1]<-n0[i]
          s2[i+1,]<-s2[i,]
          lnRg[i+1]<-lnRg[i]
        }
        i<-i+1
      }
      tc<-tc+1
    }
  }
  #Organize and return the results
  out<-list(A2=A2,lnRg=lnRg,M=M*monomer.mass,aggregate.size=M, concentration=con, dndc=dndc, n0=n0, sigma2.con=s2[,1], sigma2.con.meas=s2[,2], sigma2.scat=s2[,3], sigma2.RI=s2[,4], concentration.design=con.design, concentration.measured=con.meas, RI.measured=ri.meas, scat.ang=scat.ang, wavelength=wl, monomer.mass=monomer.mass, scattering.meas=scat.meas)
  class(out)<-"virialEst"
  out
}

#Simple plotting method for virialEst objects.
plot.virialEst<-function(x, type=c("ppc.con","ppc.RI","ppc.scat"), ...){
  #Some basic infrastructure for calculating useful things...
  qon0<-4*pi*sin(x$scat.ang/360*pi)/x$wavelength  #q/n0=4*pi*sin(theta/2)/lambda
  f<-(x$wavelength^4*6.022e23)/(4*pi^2)           #My version of const for Zimm model
  rQ<-function(lnRg,M,con,A2,d){                  #Predicted adj scattering meas
    fd<-f/((x$solvent.base.RI+d[1])*d[2])^2  #solvent RI^2 * dndc^2
    outer(con,qon0*(x$solvent.base.RI+d[1]),function(z,y){1/((1+y^2*exp(2*lnRg)/3)*(fd/(M*z)+2*fd*A2))})
  }
  #Create the requested plot
  if(match.arg(type)=="ppc.con"){ #Posterior predictive check for concentrations
    nc<-length(x$concentration.measured)
    nr<-length(x$concentration.measured.RI)
    #Set up the plot
    plot(1,0,type="n",xlim=c(1,nc+nr),ylim=range(c(as.numeric(x$concentration),as.numeric(x$concentration.RI))), xlab="Concentration Measurements",ylab="Concentration", axes=FALSE,...)
    #Plot posterior draws, measurements, and design concentrations
    for(i in 1:nc){
      points(jitter(rep(i,NROW(x$concentration))),x$concentration[,i],col=rgb(1,0,0,0.5))
      points(i,x$concentration.measured[i],pch=19)
      points(i,x$concentration.design[i],pch=19,col=rgb(0,0,1))
    }
    if(nr>0){
      for(i in 1:nr){
        points(jitter(rep(i+nc,NROW(x$concentration.RI))),x$concentration.RI[,i],col=rgb(1,0,0,0.5))
        points(i+nc,x$concentration.measured.RI[i],pch=19)
        points(i+nc,x$concentration.design.RI[i],pch=19,col=rgb(0,0,1))
      }
    }
    #Create axes
    axis(2)
    axlab<-paste("S",1:nc)
    if(nr>0)
      axlab<-c(axlab,paste("R",1:nr))
    axis(1,at=1:(nc+nr),labels=axlab,las=3)
    #Create legend
    legend("topleft",pch=c(1,19,19),col=c(rgb(1,0,0,0.5),rgb(0,0,0),rgb(0,0,1)), legend=c("Posterior Draws","Measurements","Design Values"))
  }else if(match.arg(type)=="ppc.RI"){ #Posterior predictive check for RI measurements
    nr<-length(x$RI.measured)
    #Compute the posterior predictive RI values
    if(is.null(x$concentration.RI))
      con<-x$concentration
    else
      con<-x$concentration.RI
    ppri<-sweep(sweep(con,1,x$RI.reg.param[,2],"*"),1,x$RI.reg.param[,1],"+")
    #Set up the plot
    plot(1,0,type="n",xlim=c(1,nr),ylim=range(ppri),xlab="Measurements", ylab="Refractive Index Increment",axes=F,...)
    #Plot the observations
    for(i in 1:nr){
      points(jitter(rep(i,NROW(ppri[,i]))),ppri[,i],col=rgb(1,0,0,0.5))
      points(i,x$RI.measured[i],pch=19)
    }
    #Create axes
    axis(2)
    axlab<-paste("R",1:nr)
    axis(1,at=1:nr,labels=axlab,las=3)
    #Create legend
    legend("topleft",pch=c(1,19),col=c(rgb(1,0,0,0.5),rgb(0,0,0)), legend=c("Posterior Predictive Draws","Measurements"))
  }else if(match.arg(type)=="ppc.scat"){ #Posterior predictive check for scattering
    #Compute posterior predictive scattering values
    scat<-x$scattering.meas
    nang<-length(x$scat.ang)
    nc<-NROW(scat)
    n<-length(x$lnRg)
    pprq<-array(dim=c(n,nc,nang))
    for(i in 1:n)
      pprq[i,,]<-rQ(x$lnRg[i],x$M[i],x$concentration[i,],x$A2[i],x$RI.reg.param[i,])
    #Form the plot
    pal<-hsv((0:(nang-1))/(nang-1)*0.65,alpha=0.25)
    lpal<-hsv((0:(nang-1))/(nang-1)*0.65)
    plot(1,0,type="n",xlim=c(1,nc*(nang+1)),ylim=range(pprq),xlab="", ylab=expression(Delta*R(theta,c)), axes=FALSE, ...)
    #Plot the draws  (two passes to prevent overlap issues w/many angles)
    for(i in 1:nc)
      for(j in 1:nang){
        points(jitter(rep((i-1)*(nang+1)+j,n)),pprq[,i,j],col=pal[j],cex=0.25)
      }
    for(i in 1:nc)
      for(j in 1:nang){
        points((i-1)*(nang+1)+j,scat[i,j],pch=19)
      }
    #Create axes
    axis(2)
    axlab<-paste("Conc",x$concentration.design)
    axis(1,at=(1:nc)*(nang+1)-(nang+1)/2,labels=axlab,las=3)
    #Create legend
    legend("topleft",pch=c(rep(1,nang),19),col=c(lpal,rgb(0,0,0)), legend=c(paste("PPD",x$scat.ang),"Measurements"))
  }
}


#Simulation tool for virialEst data - useful for testing virialEst
simulate.virialEst<-function(con.design, solvent.base.ri, con.design.ri=NULL, monomer.mass, scat.ang, A2, aggregate.size=1, wl=658e-7, dndc=0.188, con.lsd=log(1.05), con.meas.lsd=log(1.05), scat.sd=1e-6, ri.sd=1e-7){
  #First, define a function to produce scattering values (assume no angular effects)
  nang<-length(scat.ang)
  f<-(4*pi^2)/(wl^4*6.022e23)                     #Const part of Rayleigh "constant"
  rQ<-function(M,con,A2,dndc,n0){                 #Predicted adj scattering meas
    k<-f*(n0*dndc)^2                              #k=f* Solvent RI^2 * dndc^2
    dr<-k*con*M*monomer.mass*(1-2*A2*con*M*monomer.mass)  #Eq 8 from Wyatt, 1993
    t(sapply(dr,rep,nang))                        #Repeat for each angle
  }
  #Function for drawing from a truncated normal
  rtnorm<-function(n,mean,sd,lim=c(-Inf,Inf)){
    plim<-pnorm(lim,mean=mean,sd=sd)
    qnorm(runif(n,plim[1],plim[2]),mean=mean,sd=sd)
  }
  #Now, produce real and measured concentrations from the design concentrations
  con<-con.design*exp(rnorm(length(con.design),sd=con.lsd))
  con.meas<-con*exp(rnorm(length(con),sd=con.meas.lsd))
  if(!is.null(con.design.ri)){
    conri<-con.design.ri*exp(rnorm(length(con.design.ri),sd=con.lsd))
    con.meas.ri<-conri*exp(rnorm(length(conri),sd=sd.con.meas.lsd))
  }
  #Produce refractive index values from the real concentrations
  if(is.null(con.design.ri))
    ri.meas<-rtnorm(length(con),mean=dndc*con,sd=ri.sd,lim=c(0,Inf))
  else
    ri.meas<-rtnorm(length(conri),mean=dndc*conri,sd=ri.sd,lim=c(0,Inf))
  base.ri.meas<-rtnorm(1,mean=solvent.base.ri,sd=ri.sd,lim=c(1,Inf))
  #Produce the scattering values
  scat.meas<-apply(rQ(aggregate.size,con,A2,dndc,solvent.base.ri),c(1,2),function(z){rnorm(1,mean=z,sd=scat.sd)})
  #Finally, return everything
  out<-list(concentration=con, concentration.measured=con.meas, concentration.design=con.design, solvent.base.RI=solvent.base.ri, dndc=dndc, RI.measured=ri.meas, base.RI.measured=base.ri.meas, scattering.measured=scat.meas, scat.ang=scat.ang, aggregate.size=aggregate.size, monomer.mass=monomer.mass, M=aggregate.size*monomer.mass, A2=A2, wavelength=wl, concentration.lsd=con.lsd, concentration.measured.lsd=con.meas.lsd, scattering.sd=scat.sd, RI.sd=ri.sd)
  if(!is.null(con.design.ri)){
    out$concentration.design.RI=con.design.ri
    out$concentration.RI<-conri
  }
  out
}


#Function to perform various heuristic or other point estimates of A2 and friends,
#with uncertainties estimated by nonparametric bootstrap.  Note that we here
#assume the radius of gyration to be extremely small, and hence do not factor in
#angle-dependent effects.
virialEstBoot<-function(dat,solvent.base.RI,monomer.mass,aggregate.size=1,RI.dat=NULL,method=c("Zimm","inversion","MLE"),inv.est=c("extrapolation","median","mean"), mle.est=c("a2","a3"), wl=658e-7, use.RI.correction=TRUE, correct.RI=TRUE, bootstrap.replicates=1e3,  max.con.for.RI=Inf, min.con.for.LS=0, LS.lognorm=FALSE, monodispersed=TRUE){
  require(boot)
  if(match.arg(method)=="MLE"){
    require(trust)
    require(numDeriv)
  }
  inv.est<-match.arg(inv.est)
  #Extract data from the results object
  tab<-dat$results
  conc<-tab[,1]
  nang<-length(dat$scat.ang)
  ls<-tab[,6+((1:nang)-1)*4]
  ls<-sweep(ls,2,ls[1,],"-")
  ri<-tab[,NCOL(tab)-2]
  if(correct.RI)                 #If desired, correct RI for baseline measurement
    ri<-ri-ri[1]
  ri[conc>max.con.for.RI]<-NA    #Remove unreliable high-conc RI measurements
  ls[conc<min.con.for.LS,]<-NA   #Remove unreliable low-conc LS measurements
  if(!is.null(RI.dat)){          #Were we given additional data for RI purposes?
    ri.add<-vector()
    conc.add<-vector()
    for(i in 1:length(RI.dat)){
      conc.add<-c(conc.add,RI.dat[[i]]$results[-1,1])
      ritemp<-RI.dat[[i]]$results[,NCOL(RI.dat[[i]]$results)-2]
      if(correct.RI)
        ritemp<-ritemp-ritemp[1]
      ritemp[RI.dat[[i]]$results[,1]>max.con.for.RI]<-NA
      ri.add<-c(ri.add,ritemp[-1])
    }
    sel<-!is.na(ri.add)
    conc.add<-conc.add[sel]
    ri.add<-ri.add[sel]
  }
  scat.ang<-dat$scat.ang
  #Function to infer A2 directly in the ideal case...
  a2fun<-function(dndc,dR,m,c){
    ((2*pi*solvent.base.RI*dndc)^2/(dR*wl^4*6.022e23)-1/(monomer.mass*m*c))*0.5
  }
  #Function to calculate dR
  dRfun<-function(dndc,a2,m,c,quad.expansion=TRUE){
    #Commented out is a frequently used version - I am using a more fundamental one.
    #(2*pi*solvent.base.RI*dndc)^2/((wl^4*6.022e23)*(1/(monomer.mass*m*c)+2*a2))
    if(quad.expansion){ 
      #Eq 8 of Wyatt, 1993, following Zimm (2nd order expansion, in P->1 limit)
      (2*pi*solvent.base.RI*dndc)^2/(wl^4*6.022e23)*m*c*monomer.mass*(1-2*a2*m*monomer.mass*c)
    }else{
      #Another Zimm formula, from 2nd order expansion of P^-1(theta)
      (2*pi*solvent.base.RI*dndc)^2/((wl^4*6.022e23)*(1/(monomer.mass*m*c)+2*a2))
    }
  }
  dRfunA3<-function(dndc,a2,a3,m,c){ #From Eq9 of Wyatt, 1993, assuming P->1, Q->1
    c*(2*pi*solvent.base.RI*dndc)^2/(wl^4*6.022e23)*(monomer.mass*m)*((1-2*a2*monomer.mass*m*c)+monomer.mass*m*c^2*(4*a2^2*monomer.mass*m+3*a3))
  }
  #Function to calculate deviance related quantities for the dRI data.
  #By default, Gaussian errors are assumed (although lognormal errors may be selected).
  devdRI<-function(con,dri,dndc,s2,lognormal=FALSE){
    #Prep the data
    sel<-!(is.na(con)|is.na(dri))
    con<-con[sel]
    dri<-dri[sel]
    if((length(con)!=length(dri))||(length(con)==0))
      return(list(deviance=Inf))
    #Compute the deviance
    if(lognormal){
      deviance<--2*sum(dlnorm(dri,meanlog=log(con*dndc),sdlog=s2^0.5,log=TRUE))
    }else{
      deviance<--2*sum(dnorm(dri,mean=con*dndc,sd=s2^0.5,log=TRUE))
    }
    if(!is.finite(deviance))
      return(list(deviance=Inf))
    #Compute the gradient
    if(lognormal){
      gradient<-c(sum(2*(log(con*dndc)-log(dri))/(dndc*s2)), sum((s2-(log(con*dndc)-log(dri))^2)/s2^2))
    }else{
      gradient<-c(sum(2*con*(con*dndc-dri)/s2), sum((s2-(dri-con*dndc)^2)/s2^2))
    }
    #Compute the hessian
    if(lognormal){
      hessian<-rbind(
        c(sum(2*(1-log(con*dndc)+log(dri))/(dndc^2*s2)), sum((2*(log(dri)-log(con*dndc)))/(dndc*s2^2))),
        c(sum((2*(log(dri)-log(con*dndc)))/(dndc*s2^2)), sum((2*(log(con*dndc)-log(dri))^2-s2)/s2^3))
      )
    }else{
      hessian<-rbind(
        c(sum(2*con^2/s2), sum(2*con*(dri-con*dndc)/s2^2)),
        c(sum(2*con*(dri-con*dndc)/s2^2), sum(2*((dri-con*dndc)^2-s2)/s2^3))
      )
    }
    #Return the result
    list(deviance=deviance,gradient=gradient,hessian=hessian)
  }
  #Function to calculate deviance related quantities for the dR data under the second
  #order Zimm expansion.  M here can be treated as constant or not.  Note that we
  #assume that Rg<<wl/20, and that P(theta)->1.
  #By default, Gaussian errors are assumed (although lognormal errors may be selected).
  devdR2o<-function(con,dr,m,a2,dndc,n0,s2,const.m=TRUE,lognormal=FALSE){
    #Prep the data
    sel<-!(is.na(con)|is.na(dr))
    con<-con[sel]
    dr<-dr[sel]
    if((length(con)!=length(dr))||(length(con)==0)||(dndc<=0)||(s2<=0)||(m<1))
      return(list(deviance=Inf))
    #Compute the deviance
    drpred<-dRfun(dndc=dndc,a2=a2,m=m,c=con)
    if(lognormal&&any(drpred<=0))
      return(list(deviance=Inf))
    #plot(drpred,dr); abline(0,1)
    if(lognormal){
      deviance<- -2*sum(dlnorm(dr,meanlog=log(drpred),sdlog=s2^0.5,log=TRUE))
    }else{
      deviance<- -2*sum(dnorm(dr,mean=drpred,sd=s2^0.5,log=TRUE))
    }
    #cat("Sanity check: sd=",sqrt(s2),"ssq=",sum((dr-drpred)^2),"ss2/s2",sum((dr-drpred)^2)/s2,"mydev=",sum((dr-drpred)^2)/s2+length(dr)*log(2*pi*s2),"dev=",deviance,"\n")
    if(!is.finite(deviance))
      return(list(deviance=Inf))
    #Compute the gradient
    if(lognormal){
      lnawl<-log(6.022e23)+4*log(wl)
      lmmc<-log(con)+log(m)+log(monomer.mass)
      mmc<-exp(lmmc)
      lr1<-log(4)+lmmc+2*log(dndc*n0*pi)+log(1-2*a2*con*m*monomer.mass)-lnawl
      gradient<-c(
        sum(-4*(log(dr)-lr1)/(dndc*s2)),  #dndc
        sum(4*mmc*(log(dr)-lr1)/(s2-2*a2*s2*mmc)),  #a2
        sum(-2*(4*a2*mmc)*(log(dr)-lr1)/(m*s2*(2*a2*mmc-1))),  #m
        sum((s2-(log(dr)-lr1)^2)/s2^2)   #s2
      )
    }else{
      lnawl<-log(6.022e23)+4*log(wl)
      nawl<-exp(lnawl)
      mmc<-exp(log(con)+log(m)+log(monomer.mass))
      nawl2s2<-exp(2*lnawl+log(s2))
      gradient<-c(
        sum(16*dndc*mmc*(2*a2*mmc-1)*(n0*pi)^2*(4*(dndc*n0*pi)^2*mmc*(2*a2*mmc-1)+dr*nawl)/nawl2s2),  #dndc
        sum(16*(dndc*mmc*n0*pi)^2*(4*(dndc*n0*pi)^2*mmc*(2*a2*mmc-1)+dr*nawl)/nawl2s2),  #a2
        sum(8*exp(log(con)+2*log(dndc*n0*pi)+log(monomer.mass))*(4*a2*mmc-1)*(4*(dndc*n0*pi)^2*mmc*(2*a2*mmc-1)+dr*nawl)/nawl2s2),  #m
        sum((s2-(dr+4*(dndc*n0*pi)^2*mmc*(2*a2*mmc-1)/nawl)^2)/s2^2)   #s2
      )
    }
    if(const.m)
      gradient<-gradient[-3]
    #Compute the hessian
    if(lognormal){
      lnawl<-log(6.022e23)+4*log(wl)
      lmmc<-log(con)+log(m)+log(monomer.mass)
      mmc<-exp(lmmc)
      lr1<-log(4)+lmmc+2*log(dndc*n0*pi)+log(1-2*a2*mmc)-lnawl
      hessian<-rbind(
        c(  #dndc x 
          sum((4*(2+log(dr)-lr1))/(dndc^2*s2)), #dndc
          sum(-(8*mmc)/(dndc*s2*(1-2*a2*mmc))), #a2
          sum((4-16*a2*mmc)/(dndc*m*s2*(1-2*a2*mmc))),  #m
          sum((4*(log(dr)-lr1))/(dndc*s2^2)) #s2
        ),
        c( #a2 x
          NA, #dndc 
          sum((8*mmc^2*(1+log(dr)-lr1))/(s2*(1-2*a2*mmc)^2)), #a2
          sum((4*mmc*(4*a2*mmc+log(dr)-lr1-1))/(s2*(1-2*a2*mmc)^2)), #m
          sum((4*mmc*(log(dr)-lr1))/(s2^2*(2*a2*mmc-1)))  #s2
        ),
        c( #m x
          NA, #dndc 
          NA, #a2
          sum((2*((1-4*a2*mmc)^2+(1+4*a2*mmc*(-1+2*a2*mmc))*log(dr)+(-1+4*a2*mmc*(1-2*a2*mmc))*lr1))/(m^2*s2*(1-2*a2*mmc)^2)), #m
          sum((2*(-1+4*a2*mmc)*(log(dr)-lr1))/(m*s2^2*(-1+2*a2*mmc)))  #s2
        ),
        c( #s2 x
          NA, #dndc 
          NA, #a2
          NA, #m
          sum((-s2+2*(log(dr)-lr1)^2)/s2^3)  #s2
        )
      )
    }else{
      lmmc<-log(con)+log(m)+log(monomer.mass)
      mmc<-exp(lmmc)
      lnawl<-log(6.022e23)+4*log(wl)
      mmconawl2s2<-exp(lmmc-2*lnawl-log(s2))
      mmc2onawl2s2<-exp(2*lmmc-2*lnawl-log(s2))
      mmc2onawl2s4<-exp(2*lmmc-2*lnawl-2*log(s2))
      hessian<-rbind(
        c(  #dndc x 
          sum(16*n0^2*pi^2*mmconawl2s2*(2*a2*mmc-1)*(12*(dndc*n0*pi)^2*mmc*(2*a2*mmc-1)+dr*nawl)), #dndc
          sum(32*dndc*n0^2*pi^2*mmc2onawl2s2*(8*mmc*(dndc*n0*pi)^2*(2*a2*mmc-1)+dr*nawl)), #a2
          sum((16*dndc*n0^2*pi^2*mmc2onawl2s2*(4*a2*mmc-1)*(8*mmc*(dndc*n0*pi)^2*(2*a2*mmc-1)+dr*nawl))/m),  #m
          sum((16*dndc*n0^2*pi^2*mmconawl2s2*(2*a2*mmc-1)*(4*(dndc*n0*pi)^2*mmc*(1-2*a2*mmc)-dr*nawl))/s2) #s2
        ),
        c( #a2 x
          NA, #dndc 
          sum(exp(log(128)+4*(log(dndc)+log(n0)+log(pi)+lmmc)-2*lnawl-log(s2))), #a2
          sum((32*mmc2onawl2s2*(dndc*n0*pi)^2*(2*(dndc*n0*pi)^2*mmc*(8*a2*mmc-3)+dr*nawl))/m), #m
          sum((4*dndc*n0*pi)^2*mmc2onawl2s4*(4*(dndc*n0*pi)^2*mmc*(1-2*a2*mmc)-dr*nawl))  #s2
        ),
        c( #m x
          NA, #dndc 
          NA, #a2
          sum((32*mmc2onawl2s2*(dndc*n0*pi)^2*((dndc*n0*pi)^2*(1+12*a2*mmc*(2*a2*mmc-1))+a2*dr*nawl))/(m^2)), #m
          sum(-exp(lmmc-log(m)-2*lnawl-2*log(s2))*8*(dndc*n0*pi)^2*(4*a2*mmc-1)*((2*dndc*n0*pi)^2*mmc*(2*a2*mmc-1)+dr*nawl))  #s2
        ),
        c( #s2 x
          NA, #dndc 
          NA, #a2
          NA, #m
          sum((-s2+2*(dr+((2*dndc*n0*pi)^2*mmc*(2*a2*mmc-1))/nawl)^2)/s2^3)  #s2
        )
      )
    }
    hessian[lower.tri(hessian)]<-t(hessian)[lower.tri(hessian)]
    if(const.m)
      hessian<-hessian[-3,-3]
    #Return the result
    list(deviance=deviance,gradient=gradient,hessian=hessian)
  }
  #Function to obtain estimates by a semi-traditional application of Zimm's method
  #(regression of Kc/dR on the RHS of Zimm's formula).  The main exception is that
  #we optionally constrain aggregate size to be an integer, and optimize choice of M by 
  #direct search.  (In the polydispersed case, we instead perform a search over the
  #weighted mean oligomer size, in both cases minimizing the squared error of the 
  #scattering measurements.)  The data set here should be concentration, dRI, and LS
  #vectors (in that order), with the zero concentration removed.
  zimmest<-function(dat,ind,conc.corr=FALSE){
    #First, estimate dndc
    x<-dat[ind,1]             #Concentration
    y<-dat[ind,2]             #dRI
    dndc<-coef(lm(y~x-1))[1]  #Regress dRI on concentration
    #If desired, correct concentrations using the estimated RI increment
    if(conc.corr)
      x[!is.na(y)]<-predict(lm(x~y-1))
    #Now, estimate A2 at each non-zero concentration on each angle
    k<-(2*pi*solvent.base.RI*dndc)^2/(wl^4*6.022e23)
    sel<-!is.na(dat[ind,3])                       #Select only cases with LS measurements
    y<-k*rep(x[sel],nang)/unlist(dat[ind[sel],-(1:2)]) #Kc/dR
    x2<-2*rep(x[sel],nang)                             #2c
    if(monodispersed){   #Sample is monodisperse - search given sizes for best fit
      a2est<-vector()
      lsrmse<-vector()
      for(i in 1:length(aggregate.size)){
        ym<-y-1/(monomer.mass*aggregate.size[i])  #Adjust for molecular weight
        a2est[i]<-coef(lm(ym~x2-1))[1]            #Slope estimates A2
        dRpred<-sapply(1:nang,function(z,a2,m){   #Get the RMSE 
          dRfun(dndc,a2,m,x[sel])
        },m=aggregate.size[i],a2=a2est[i])
        lsrmse[i]<-mean((as.matrix(dRpred-dat[ind[sel],-(1:2)])^2))^0.5
      }
      lsrmse[is.na(lsrmse)]<-Inf
      sel<-which(lsrmse==min(lsrmse,na.rm=TRUE))[1]
      m<-aggregate.size[sel]
      a2<-a2est[sel]
    }else{  #Sample is polydispersed - search over input range
      #Get the best fitting size-weighted mean monomer size
      m<-optimize(function(m){
          ym<-y-1/(monomer.mass*m)                  #Adjust for molecular weight
          a2est<-coef(lm(ym~x2-1))[1]            #Slope estimates A2
          dRpred<-sapply(1:nang,function(z,a2,m){   #Get the RMSE 
            dRfun(dndc,a2,m,x[sel])
          },m=m,a2=a2est)
          lsrmse<-mean((as.matrix(dRpred-dat[ind[sel],-(1:2)])^2))^0.5
          if(is.na(lsrmse))
            lsrmse<-1e300
          lsrmse                    #Return the RMSE
        },lower=max(1,min(aggregate.size)),upper=max(aggregate.size))$minimum
      #Now, find the a2 estimate at m
      ym<-y-1/(monomer.mass*m)                  #Adjust for molecular weight
      a2<-coef(lm(ym~x2-1))[1]                  #Slope estimates A2
    }
    out<-c(a2,m,dndc)
    names(out)<-c("A2","M","dndc")
    out
  }
  #Function to obtain estimates by inversion (solving Zimm's formula for A2 at each
  #concentration and finding the aggregate size whose mean A2 prediction minimizes the
  #squared error).  The data set here should be concentration, dRI, and LS vectors (in
  #that order), with the zero concentration removed.
  invest<-function(dat,ind,conc.corr=FALSE){
    #First, estimate dndc
    x<-dat[ind,1]          #Concentration
    y<-dat[ind,2]          #dRI
    dndc<-coef(lm(y~x-1))  #Regress dRI on concentration
    #If desired, correct concentrations using the estimated RI increment
    if(conc.corr)
      x[!is.na(y)]<-predict(lm(x~y-1))
    #Now, estimate A2 at each non-zero concentration on each angle
    a2est<-vector()
    lsmae<-vector()
    sel<-!is.na(dat[ind,3])                      #Select cases w/LS data
    for(i in 1:length(aggregate.size)){
      a2inv<-sapply(1:nang,function(z,m){
        a2fun(dndc,dat[ind[sel],2+z],m,x[sel])
      },m=aggregate.size[i])
      if(inv.est=="median")
        a2est[i]<-median(a2inv)
      else if(inv.est=="mean")
        a2est[i]<-mean(a2inv)
      else if(inv.est=="extrapolation"){
        lx<-1/rep(x[sel],nang)
        ly<-as.numeric(a2inv)
        a2est[i]<-coef(lm(ly~lx+I(lx^2)))[1]
      }
      dRpred<-sapply(1:nang,function(z,a2,m){
        dRfun(dndc,a2,m,x[sel])
      },m=aggregate.size[i],a2=a2est[i])
      lsmae[i]<-median(abs(as.matrix(dRpred-dat[ind[sel],-(1:2)])))
    }
    lsmae[is.na(lsmae)]<-Inf
    sel<-which(lsmae==min(lsmae,na.rm=TRUE))[1]
    m<-aggregate.size[sel]
    a2<-a2est[sel]
    out<-c(a2,m,dndc)
    names(out)<-c("A2","M","dndc")
    out
  }
  #Function to obtain estimates by maximum likelihood (jointly fitting dndc, M,
  #and A2 assuming Gaussian errors in the RI and scattering data).  The data set here
  #should be concentration, dRI, and LS vectors (in that order), with the zero 
  #concentration removed.
  mlest<-function(dat,ind,useA3=TRUE,conc.corr=FALSE,lognorm.ls=TRUE,monodispersed=TRUE){
    #First, create the basic prediction variables
    xri<-dat[ind,1]                #Concentration for RI
    yri<-dat[ind,2]                #dRI
    xls<-xri
    if(conc.corr)                  #If using correction, correct w/RI vals
      xls[!is.na(yri)]<-predict(lm(xri~yri))
    sel<-!is.na(dat[ind,3])             #Get concentrations with LS data
    xls<-rep(xls[sel],nang)             #Concentration for LS
    yls<-unlist(dat[ind[sel],-(1:2)])   #dR
    dndc0<-coef(lm(yri~xri-1))[1]       #Seed value for dndc
    #Define a function to calculate the deviance; par=(a2,dndc,s2ri,s2ls,A3,m)
    dev<-function(par,mm,xri,yri,xls,yls,lognorm.ls=TRUE,useA3=TRUE,monodispersed=TRUE){
      if(useA3)
        A3<-par[5]
      if(!monodispersed)
        m<-par[5+useA3]
      else
        m<-mm
      devRI<-devdRI(con=xri,dri=yri,dndc=par[2],s2=par[3],lognormal=FALSE) #dndc, s2ri
      #print(devRI)
      if(is.finite(devRI$deviance)){
        if(useA3)
          devLS<-devdR3o(con=xls,dr=yls,m=m,a2=par[1],dndc=par[2],n0=solvent.base.RI,s2=par[4],a3=A3,const.m=monodispersed,lognormal=lognorm.ls) #Returns as dndc, a2, m, s2ls, a3
        else
          devLS<-devdR2o(con=xls,dr=yls,m=m,a2=par[1],dndc=par[2],n0=solvent.base.RI,s2=par[4],const.m=monodispersed,lognormal=lognorm.ls) #Returns as dndc, a2, m, s2ls
      }else
        return(list(value=Inf))
      #print(devLS)
      value<-devRI$deviance+devLS$deviance
      if(!is.finite(value))
        return(list(value=Inf))
      gradient<-devLS$gradient
      gradient[1]<-gradient[1]+devRI$gradient[1]
      gradient<-c(gradient[1:2],devRI$gradient[2],gradient[-(1:2)])
      hessian<-rbind(devLS$hessian[1:2,],rep(0,NCOL(devLS$hessian)),devLS$hessian[-(1:2),])
      hessian<-cbind(hessian[,1:2],rep(0,NCOL(devLS$hessian)+1),hessian[,-(1:2)])
      hessian[c(1,3),c(1,3)]<-hessian[c(1,3),c(1,3)]+devRI$hessian
      #Order now dndc, a2, s2ri, m, s2ls, a3; map to a2, dndc, s2ri, s2ls, m, a3
      if(useA3){
        if(monodispersed){
          gradient<-gradient[c(2,1,3,4,5)]
          hessian<-hessian[c(2,1,3,4,5),c(2,1,3,4,5)]
        }else{
          gradient<-gradient[c(2,1,3,5,6,4)]
          hessian<-hessian[c(2,1,3,5,6,4),c(2,1,3,5,6,4)]
        }
      }else{
        if(monodispersed){  #dndc,a2,s2ri,s2ls -> a2,dndc,s2ri,s2ls
          gradient<-gradient[c(2,1,3,4)]
          hessian<-hessian[c(2,1,3,4),c(2,1,3,4)]
        }else{               #dndc,a2,s2ri,m,s2ls -> a2,dndc,s2ri,s2ls,m
          gradient<-gradient[c(2,1,3,5,4)]
          hessian<-hessian[c(2,1,3,5,4),c(2,1,3,5,4)]
        }
      }
      #cat("COMBINED DEV AT PAR",par,"\n")
      #print(value)
      #print(gradient)
      #print(hessian)
      list(value=value,gradient=gradient,hessian=hessian)
    }
    devold<-function(par,m,lognorm.ls=TRUE,useA3=TRUE,monodispersed=TRUE){
      #cat("entering dev:",par,"\n")
      #Deviances within deviances....
      intdev<-function(par,m,lognorm.ls=TRUE,useA3=TRUE,monodispersed=TRUE){
        #cat("\tin intdev:",par,"\n")
        if(useA3)
          A3<-par[5]
        else
          A3<-0
        if(!monodispersed)
          m<-par[5+useA3]
        dRIpred<-xri*par[2]
        if(useA3)
          dRpred<-dRfunA3(par[2],par[1],A3,m,xls)
        else
          dRpred<-dRfun(par[2],par[1],m,xls)
        #par(mfrow=c(1,2)); plot(dRIpred,yri); abline(0,1); plot(dRpred,yls); abline(0,1)
        if(any(is.na(dRpred))||any(dRpred<=0)||any(par[2:4]<=1e-2)||(m<1)|| ((!monodispersed)&&(par[5+useA3]<1))){
          Inf
        }else{
          if(lognorm.ls)
            -2*(sum(dnorm(yri,mean=dRIpred,sd=par[4]*1e-7,log=TRUE),na.rm=TRUE) + sum(dlnorm(yls,meanlog=log(dRpred),sdlog=par[3]*1e-7,log=TRUE),na.rm=TRUE))
          else
            -2*(sum(dnorm(yri,mean=dRIpred,sd=par[4]*1e-7,log=TRUE),na.rm=TRUE) + sum(dnorm(yls,mean=dRpred,sd=par[3]*1e-7,log=TRUE),na.rm=TRUE))
        }
      }
      deviance<-intdev(par,m=m,lognorm.ls=lognorm.ls,useA3=useA3, monodispersed=monodispersed)
      if(is.finite(deviance)){
        gradient<-grad(function(p,as,lognorm.ls,useA3,monodispersed){
          min(intdev(p,m=as,lognorm.ls=lognorm.ls,useA3=useA3,monodispersed=monodispersed), 1e200)
        },par,as=m,lognorm.ls=lognorm.ls,useA3=useA3,monodispersed=monodispersed)
      }else
        gradient<-NULL
      if(is.finite(deviance)){
        hess<-hessian(function(p,as,lognorm.ls,useA3,monodispersed){
          min(intdev(p,m=as,lognorm.ls=lognorm.ls,useA3=useA3,monodispersed=monodispersed), 1e200)
        },par,as=m,lognorm.ls=lognorm.ls,useA3=useA3,monodispersed=monodispersed)
      }else
        hess<-NULL 
      list(value=deviance,gradient=gradient,hessian=hess)
    }
    #Now, estimate for each aggregate size
    estvals<-vector()
    devvals<-vector()
    if(monodispersed){
      for(i in 1:length(aggregate.size)){
        #locfit<-optim(c(0,dndc0,1,1), dev, m=aggregate.size[i], method="L-BFGS-B", lower=c(-Inf,1e-9,1e-12,1e-12))
        #par<-locfit$par
        #par[3:4]<-par[3:4]*1e-3  #Correct values
        ipar<-c(0,dndc0,var(yri,na.rm=TRUE),var(yls,na.rm=TRUE))
        if(useA3)
          ipar<-c(ipar,0)
        if(lognorm.ls)
          ipar[4]<-var(log(yls),na.rm=TRUE)
        parscale<-rep(1,length(ipar))
        parscale[3:4]<-1e3
        locfit<-trust(dev,ipar,10,100,parscale,mterm=sqrt(.Machine$double.eps), minimize=TRUE, xri=xri,yri=yri,xls=xls,yls=yls,mm=aggregate.size[i],useA3=useA3, lognorm.ls=lognorm.ls, monodispersed=monodispersed)
        #locfit<-trust(olddev,ipar,10,100,rep(1,length(ipar)),mterm=sqrt(.Machine$double.eps), minimize=TRUE, m=aggregate.size[i],useA3=useA3,lognorm.ls=lognorm.ls, monodispersed=monodispersed)
        par<-locfit$argument
        #par[3:4]<-par[3:4]*1e-7  #Correct values
        estvals<-cbind(estvals,par)
        devvals[i]<-locfit$value
      }
      sel<-which(devvals==min(devvals,na.rm=TRUE))[1]
      m<-aggregate.size[sel]
      est<-estvals[,sel]
      out<-c(est,m)
    }else{
      ipar<-c(0,dndc0,var(yri,na.rm=TRUE),var(yls,na.rm=TRUE))
      if(useA3)
        ipar<-c(ipar,0)
      ipar<-c(ipar,1.05)
      if(lognorm.ls)
        ipar[4]<-var(log(yls),na.rm=TRUE)
      parscale<-rep(1,length(ipar))
      parscale[3:4]<-1e3
      parscale[5]<-1e-3
      locfit<-trust(dev,ipar,10,100,parscale,mterm=sqrt(.Machine$double.eps), minimize=TRUE, xri=xri,yri=yri,xls=xls,yls=yls,mm=1.05,useA3=useA3, lognorm.ls=lognorm.ls, monodispersed=monodispersed)
      #locfit<-trust(olddev,ipar,10,100,rep(1,length(ipar)),mterm=sqrt(.Machine$double.eps), minimize=TRUE, m=1, useA3=useA3, lognorm.ls=lognorm.ls, monodispersed=monodispersed)
      par<-locfit$argument
      #par[3:4]<-par[3:4]*1e-7  #Correct values
      out<-par
    }
    if(!useA3)
      names(out)<-c("A2","dndc","sigma2.RI","sigma2.LS","M")
    else
      names(out)<-c("A2","dndc","sigma2.RI","sigma2.LS","A3","M")
    out
  }
  #Perform inference, with bootstrap replication for SEs and such
  if(match.arg(method)=="Zimm"){  #Zimm's method(ish)
    #Create the data table for bootstrapping
    dtab<-cbind(conc[-1],ri[-1],ls[-1,])
    if(!is.null(RI.dat))
      dtab<-rbind(dtab,cbind(conc.add,ri.add,matrix(NA,length(conc.add),NCOL(ls))))
    dtab<-as.data.frame(dtab)
    #Perform estimation using a non-parametric bootstrap on the non-zero conc vals
    if(!is.null(RI.dat))
      fit<-boot(dtab,zimmest,bootstrap.replicates,strata=rep(1:2,times=c(length(conc)-1,length(ri.add))),conc.corr=use.RI.correction)
    else
      fit<-boot(dtab,zimmest,bootstrap.replicates,conc.corr=use.RI.correction)
  }else if(match.arg(method)=="inversion"){  #Direct conditional inversion
    #Create the data table for bootstrapping
    dtab<-cbind(conc[-1],ri[-1],ls[-1,])
    if(!is.null(RI.dat))
      dtab<-rbind(dtab,cbind(conc.add,ri.add,matrix(NA,length(conc.add),NCOL(ls))))
    dtab<-as.data.frame(dtab)
    #Perform estimation using a non-parametric bootstrap on the non-zero conc vals
    if(!is.null(RI.dat))
      fit<-boot(dtab,invest,bootstrap.replicates,strata=rep(1:2,times=c(length(conc)-1,length(ri.add))),conc.corr=use.RI.correction)
    else
      fit<-boot(dtab,invest,bootstrap.replicates,conc.corr=use.RI.correction)
  }else if(match.arg(method)=="MLE"){   #Maximum likelihood estimation
    #Create the data table for bootstrapping
    dtab<-cbind(conc[-1],ri[-1],ls[-1,])
    if(!is.null(RI.dat))
      dtab<-rbind(dtab,cbind(conc.add,ri.add,matrix(NA,length(conc.add),NCOL(ls))))
    dtab<-as.data.frame(dtab)
    #Perform estimation using a non-parametric bootstrap on the non-zero conc vals
    if(!is.null(RI.dat))
      fit<-boot(dtab,mlest,bootstrap.replicates, strata=rep(1:2,times=c(length(conc)-1,length(ri.add))),useA3=(match.arg(mle.est)=="a3"),lognorm.ls=LS.lognorm, conc.corr=use.RI.correction, monodispersed=monodispersed)
    else
      fit<-boot(dtab,mlest,bootstrap.replicates,useA3=(match.arg(mle.est)=="a3"),lognorm.ls=LS.lognorm,conc.corr=use.RI.correction, monodispersed=monodispersed)
  }
  #Compute residuals and such
  sel<-!is.na(dtab[,2])
  obsdri<-dtab[sel,2]                 #dRI observations (as passed to )
  predri<-dtab[sel,1]*fit$t0["dndc"]  #dRI predictions
  resdri<-obsdri-predri               #dRI residuals
  condri<-dtab[sel,1]                 #dRI concentrations
  sel<-!is.na(dtab[,3])
  obsdr<-unlist(dtab[sel,-(1:2)])
  condr<-rep(dtab[sel,1],NCOL(dtab)-2)
  if(match.arg(mle.est)=="a3")
    preddr<-dRfunA3(fit$t0["dndc"],fit$t0["A2"],fit$t0["A3"],fit$t0["M"],condr)
  else
    preddr<-dRfun(fit$t0["dndc"],fit$t0["A2"],fit$t0["M"],condr)
  resdr<-obsdr-preddr
  #Return the results
  out<-list(fit=fit,pred.RI=predri,pred.LS=preddr,resid.RI=resdri,resid.LS=resdr, obs.RI=obsdri,obs.LS=obsdr,con.RI=condri,con.LS=condr)
  out
}

