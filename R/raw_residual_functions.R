resid_test<-function(object,...){
  z<-object$derived$psi_c$est
  psi<-object$real$psi$est

  occ.resid<-(z-psi)/sqrt((psi*(1-psi)))

  y<-unlist(object$data$det.data)
  p<-object$real$p$est

  det.resid<-(y-p)/sqrt(p*(1-p))
  det.resid<-z*det.resid
  det.resid[is.na(det.resid)]<-0

  return(list(occ=occ.resid,det=det.resid))
}


bootData<-function(model){
  psi<-model$real$psi$est
  p<-model$real$p$est
  z<-rbinom(length(psi),1,psi)
  y<-rbinom(length(p),1,p)
  det<-z*y
  is.na(det)<-(1:length(p))*is.na(model$data$det.data)
  return(det)
}

bootstrapResiduals<-function(model,nboot,pao.data){
  ## pull out model defn
  mod.defn<-list(as.formula(paste0("psi",deparse(model$model$psi,width.cutoff = 500))),
                 as.formula(paste0("p",deparse(model$model$p,width.cutoff = 500))))
  boot<-lapply(1:nboot,function(ii){
    bootResid(model,mod.defn,pao.data=pao.data,ii)
  })

  temp<-do.call(rbind,boot)
  occ.resid<-t(do.call(rbind,temp[,"occ"]))
  det.resid<-t(do.call(rbind,temp[,"det"]))
  num.det<-t(do.call(rbind,temp[,"num.det"]))

  return(list(occ=occ.resid,det=det.resid,num.det=num.det))
}

bootResid<-function(model,mod.defn,pao.data,count) {
  ok<-FALSE
  while(!ok){
    data<-bootData(model)
    det.data<-array(data,dim=dim(model$data$det.data))

    num.det<-rowSums(det.data,na.rm=TRUE)

    bootData<-pao.data
    bootData$det.data<-det.data

    paoname<-paste0("bootdata_",count,".pao")
    modname<-paste0("bootmod_",count,".out")


    fit<-occMod(model=mod.defn,data=bootData,type="so",
                 paoname=paoname,modname = modname)
    ok<-is.null(fit$warnings$VC)
  }


  psi<-fit$real$psi$est
  p<-fit$real$p$est
  is.na(p)<-(1:length(p))*is.na(data)
  pMat<-array(1-p,dim=dim(model$data$det.data))

  detected<-apply(det.data,1,function(row)return(sum(row,na.rm=TRUE)>0))

  p0<-apply(pMat,1,prod,na.rm=TRUE)
  z<-psi*p0/(1-psi*(1-p0))
  z[detected]<-1

  occ.resid<-(z-psi)/sqrt((psi*(1-psi)))

  y<-data
  det.resid<-(y-p)/sqrt(p*(1-p))
  det.resid<-z*det.resid
  det.resid[is.na(det.resid)]<-0

  return(list(occ=occ.resid,det=det.resid,num.det=num.det))
}

cusumPlot<-function(x,resid,boot.resid,jitter=TRUE,plot.series=NULL,...){

  if(class(x)=="character")x<-as.factor(x)

  factor.levels<-NULL
  if(is.factor(x)) {
    factor.levels<-levels(x)
    x <-as.numeric(x)
  }

  if(jitter){
    ## determine smallest increment of unique x values
    uni<-sort(unique(x))
    incr<-min(uni[-1]-uni[-length(uni)])
    jit<-runif(length(x),-0.25,0.25)*incr
    order<-order(x,jit)
    x<-x+jit
    x<-x[order]
    obs.cusum<-cumsum(resid[order])
    boot.cusum<-apply(boot.resid[order,],2,cumsum)
  } else {
    ## take means at tied values instead
    means<-tapply(resid,list(x),mean,na.rm=TRUE)
    obs.cusum<-cumsum(means)
    boot.cusum<-apply(boot.resid,2,function(col){
      temp<-tapply(col,list(x),mean)
      return(cumsum(temp))
    })
    x<-sort(unique(x))
  }

  quan.cusum<-t(apply(boot.cusum,1,quantile,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE))

  quan.col=c(grey(0.8),grey(0.6),grey(0.4),grey(0.2),grey(0.4),grey(0.6),grey(0.8))

#  plot(x,obs.cusum,type="n",ylim=range(cbind(obs.cusum,boot.cusum),na.rm=TRUE),...)
  plot(x,obs.cusum,type="n",...)
  for (ii in 1:7) points(x,quan.cusum[,ii],type="l",lwd=2,col=quan.col[ii]) # min

  if(!is.null(plot.series)) {
    for(ii in 1:plot.series) points(x,boot.cusum[,ii],type="l",col="blue",lwd=0.5)
  }
  points(x,obs.cusum,type="l",lwd=2,col="red")
}


movAvgPlot<-function(x,resid,boot.resid,jitter=TRUE,plot.series=NULL,window=5,...){

  if(class(x)=="character")x<-as.factor(x)

  factor.levels<-NULL
  if(is.factor(x)) {
    factor.levels<-levels(x)
    x <-as.numeric(x)
  }

  if(jitter){
    ## determine smallest increment of unique x values
    uni<-sort(unique(x))
    incr<-min(uni[-1]-uni[-length(uni)])
    jit<-runif(length(x),0,1)*incr
    order<-order(x,jit)
    x<-x+jit
  } else {
    order<-order(x)
  }

  x<-x[order]
  x<-filter(x,rep(1/window,window))
  obs.ma<-filter(resid[order],rep(1/window,window))
  boot.ma<-apply(boot.resid[order,],2,function(col) filter(col,rep(1/window,window)))
  quan.ma<-t(apply(boot.ma,1,quantile,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE))

  quan.col=c(grey(0.8),grey(0.6),grey(0.4),grey(0.2),grey(0.4),grey(0.6),grey(0.8))
  plot(x,obs.ma,type="n",ylim=range(cbind(obs.ma,boot.ma),na.rm=TRUE),...)
  for (ii in 1:7) points(x,quan.ma[,ii],type="l",lwd=2,col=quan.col[ii]) # min

  if(!is.null(plot.series)) {
    for(ii in 1:plot.series) points(x,boot.ma[,ii],type="l",col="blue",lwd=0.5)
  }
  points(x,obs.ma,type="l",lwd=2,col="red")
}


sim_occ_mod_so<-function(model,unitcov,survcov,psi.coeff,p.coeff){
  psi<-plogis(model.matrix(model$psi,data=unitcov)%*%as.matrix(psi.coeff))
  p<-plogis(model.matrix(model$p,data=cbind(unitcov,survcov))%*%as.matrix(p.coeff))

  s<-length(psi)
  k<-length(p)/s

  z<-rbinom(s,1,psi)
  y<-rbinom(s*k,1,z*p)

  det.data<-as.data.frame(array(y,dim=c(s,k)))
  return(det.data)
}



cusumNDPlot<-function(ND,resid,boot.resid,boot.ND,jitter=TRUE,plot.series=NULL,...){

  k<-length(resid)/length(ND)
  ND<-rep(ND,k)

  if(k>1) {
    temp<-boot.ND
    for (kk in 2:k){
      temp<-rbind(temp,boot,ND)
    }
    boot.ND<-temp
  }

  if(class(ND)=="character")ND<-as.factor(ND)

  factor.levels<-NULL
  if(is.factor(ND)) {
    factor.levels<-levels(ND)
    ND <-as.numeric(ND)
  }

  if(jitter){
    ## determine smallest increment of unique ND values
    uni<-sort(unique(ND))
    incr<-min(uni[-1]-uni[-length(uni)])
    jit<-runif(length(ND),-0.25,0.25)*incr
    order<-order(ND,jit)
    ND<-ND+jit
    ND<-ND[order]
    obs.cusum<-cumsum(resid[order])
    boot.cusum<-sapply(1:ncol(boot.resid),function(ii){
      order2<-order(boot.ND[,ii],jit)
      boot.ND[,ii]<-boot.ND[,ii]+jit
      boot.ND[,ii]<-boot.ND[,ii][order2]
      return(cumsum(boot.resid[order2,ii]))
    })
  } else {
    ## take means at tied values instead
    means<-tapply(resid,list(ND),mean,na.rm=TRUE)
    obs.cusum<-cumsum(means)
    boot.cusum<-apply(boot.resid,2,function(col){
      temp<-tapply(col,list(ND),mean)
      return(cumsum(temp))
    })
    ND<-sort(unique(ND))
  }

  quan.cusum<-t(apply(boot.cusum,1,quantile,probs=c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE))

  quan.col=c(grey(0.8),grey(0.6),grey(0.4),grey(0.2),grey(0.4),grey(0.6),grey(0.8))

  plot(ND,obs.cusum,type="n",ylim=range(cbind(obs.cusum,boot.cusum),na.rm=TRUE),...)
  for (ii in 1:7) points(ND,quan.cusum[,ii],type="l",lwd=2,col=quan.col[ii]) # min

  if(!is.null(plot.series)) {
    for(ii in 1:plot.series) points(ND,boot.cusum[,ii],type="l",col="blue",lwd=0.5)
  }
  points(ND,obs.cusum,type="l",lwd=2,col="red")
}

