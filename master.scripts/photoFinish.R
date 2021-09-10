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


