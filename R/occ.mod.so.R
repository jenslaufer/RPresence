occMod_SO<-function(psi=call(),psi.cov=data$unitcov,
                     p=call(),p.cov=cbind(data$unitcov,data$survcov),
                     modname=NULL,paoname=NULL,outfile,model=100,fixed=NULL,initvals=NULL,data,conf,miscopts){
  #' Fit static occupancy/single season model
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="so"}.
  #'
  #'@param psi the right-hand side of the formula for the model to fit for occupancy probability.
  #'@param psi.cov a data frame containing the unit-specific covariates to use for the occupancy
  #'component of the model.
  #'@param p the right-hand side of the formula for the model to fit for detection probability.
  #'@param p.cov a data frame containing the survey-specific covariates to use for the detection
  #'component of the model.
  #'@param modname (optional) a string containing the model name
  #'@param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #'@param outfile name for output file (use outfile='modname') for outfile named via model name
  #'@param model the PRESENCE model code. DO NOT CHANGE.
  #'@param fixed a single-column matrix containing values for real parameters to be fixed at.
  #'\code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #'@param initvals initial values for the beta parameters at which PRESENCE begins the optimisation.
  #'The default values in PRESENCE is 0.
  #'@param data the \code{pao} data object containing the detection data and other information.
  #'@param conf level for confidence interval (may be vector valued).
  #'@param miscopts
  #'
  #'@return returns a list of class \code{"occMod"} and \code{"so"}.
  #'
  #'\code{occMod$beta} contains the objects:
  #'  #'\item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence.}
  #'\item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence.}
  #'\item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'\item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'\item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'
  #'\code{occMod$real} contains the objects:
  #'\item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard
  #'errors and limits of 95\% confidence interval.}
  #'\item{p}{estimated probabilities of detection for each survey, along with standard errors and
  #'limits of 95\% confidence interval.}
  #'
  #'\code{occMod$derived} contains the objects:
  #'\item{psi_c}{estimated probabilities of occurrence given the detection history for each sampling
  #'unit, along with standard errors and limits of 95\% confidence interval. Will be \code{=1} for
  #'any unit where the species was detected at least once.}
  #'
  #'@author Darryl MacKenzie
  #'
  #'@seealso \code{\link{occMod}}, \code{pao}

  options("na.action"="na.pass")
  if (is.null(psi.cov)) psi_mat=model.matrix(psi,data.frame(rep(1,data$nunits))) else psi_mat=model.matrix(psi,psi.cov)
  if (is.null(p.cov)) p_mat=model.matrix(p,data.frame(SURVEY=as.factor(rep(1:data$nsurveys,each=data$nunits)))) else
                      p_mat<-model.matrix(p,p.cov)

  ## create and output temporary pao file
  temp.pao=data
  temp.pao$nunitcov=ncol(psi_mat); temp.pao$unitcov=as.data.frame(psi_mat)
  temp.pao$nsurvcov=ncol(p_mat)
  temp.pao$survcov=as.data.frame(p_mat)
  temp.pao$paoname=ifelse(is.null(paoname),"paodata.pao",paoname)
  #  temp.pao$paoname="test.pao"

  colnames(temp.pao$unitcov)=gsub("[(]Intercept[)]","int",colnames(temp.pao$unitcov))
  colnames(temp.pao$unitcov)=paste("psi.",colnames(temp.pao$unitcov),sep="")
  colnames(temp.pao$survcov)=gsub("[(]Intercept[)]","int",colnames(temp.pao$survcov))
  colnames(temp.pao$survcov)=paste("p.",colnames(temp.pao$survcov),sep="")

  #writePao(temp.pao)

  ## create design matrices file
  psi.dm=matrix(colnames(temp.pao$unitcov),1,temp.pao$nunitcov)
  psi.dm.col=ncol(psi.dm)
  rownames(psi.dm)="psi"
  colnames(psi.dm)=paste("a",1:temp.pao$nunitcov,sep="")

  p.dm=matrix(rep(colnames(temp.pao$survcov),temp.pao$nsurveys),
               nrow=temp.pao$nsurveys,ncol=temp.pao$nsurvcov,byrow=TRUE);
  p.dm.col=ncol(p.dm)
  rownames(p.dm)=paste0("p",1:temp.pao$nsurveys)
  colnames(p.dm)=paste0("b",1:ncol(p.dm))
  npar=psi.dm.col+p.dm.col
  psi.dm=simplify_dm(psi.dm,temp.pao);  p.dm=simplify_dm(p.dm,temp.pao)

  if(!is.null(fixed)){
    rnames=unlist(lapply(list(psi.dm,p.dm),rownames))
    fixed$idx<-match(fixed$param,rnames)
    if (sum(is.na(fixed$idx))>0) cat('fixed names must be in ',rnames,'\n')
  }
  rv=runPresence(temp.pao,list(psi.dm,p.dm,NULL,NULL,NULL,NULL),
                model,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)

  ### extract results from output file
  v=readLines(outfile)

  warn.conv=check_conv_warn(v); warn.VC=check_VC_warn(v)   ##### check for warnings

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  v=v[-1:-i]

  ### get Untransformed coefficients (Beta's)
  i=grep("^A.+psi.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
  j=grep("^B.+ : ",v); p.coeff=getbetas(v[j],p.dm)

  ### get_VC matrix
  ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)
  instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC=as.matrix(VC[1:length(i), 1:length(i)])
  p.VC  =as.matrix(VC[length(i)+1:length(j), length(i)+1:length(j)])

  psi.est=p.est=psi_c.est=NULL
  VCoutopt=miscopts[1]
  if (! (VCoutopt %in% c("nose","betavc","noreal"))) {  ## return real parameters
    ## Calc real parameter estimates
    psi.est=calc_real(cov=psi_mat,coeff=psi.coeff$est,VC=psi.VC,conf=conf,rownames=paste0('psi_',temp.pao$unitnames))
    s=paste0('p',rep(1:temp.pao$nsurveys,each=temp.pao$nunits),'_',rep(temp.pao$unitnames,temp.pao$nsurveys))
    p.est=calc_real(cov=p_mat,coeff=p.coeff$est,VC=p.VC,conf=conf,rownames=s,fixed=fixed)

    psi_c.est=calc_psi_c(psi.cov=psi_mat,psi.coeff=psi.coeff$est,p.cov=p_mat,p.coeff=p.coeff$est,
                  VC=as.matrix(VC),det.data=temp.pao$det.data,
                  rownames=temp.pao$unitnames,conf=conf)
  } # end if(!noReal)
  #  get gof results if desired
  gof=list(chat=NA,TS=NA,tbl=NULL,bs=NA)
  if (length(grep('boot2=',miscopts))>0) {
    i=grep('boot2=',miscopts); modfitboot=as.numeric(gsub('.+=','',miscopts[i]))
    i=grep("  Chi-square",v); j=grep("Test Statistic =",v); n=j-i-1
    v2=gsub(' +',' ',v[i+1:n])
    v5=lapply(1:n,function(ii){
      v3=unlist(strsplit(v2[ii],')',fixed=T));
      v4=c(paste0(v3[1],')'),as.numeric(unlist(strsplit(v3[2],' '))[-1]))
    })
    Prb=as.numeric(gsub('.+= +','',v[grep("stic >= obs",v)]))
    TS=as.numeric(gsub(' min.+','',gsub('.+= +','',v[grep("tic .data",v)])))
    TS_low=as.numeric(gsub('.+= +','',v[grep("Lowest",v)]))
    TS_hi=as.numeric(gsub('.+= +','',v[grep("Highest",v)]))
    TS=c(TS,TS_low,TS_hi); names(TS)=c('TestStat','TS_low','TS_High')
    j=grep('of c-hat',v); chat=as.numeric(gsub('.+= +','',gsub(' +[(]=.+','',v[j])))
    tbl=matrix(unlist(v5),ncol=4,byrow=T); colnames(tbl)=c('History','Observed','Expected','Chi-square')
    gof=list(chat=chat,TS=TS,Prb=Prb,tbl=tbl,bootstraps=modfitboot)
  }


  result=list(modname=modname, model=list(psi=psi,p=p), dmat=list(psi=psi.dm,p=p.dm),
               data=temp.pao, outfile=outfile, neg2loglike=neg2loglike, npar=npar, aic=aic,
               beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                         p=as.data.frame(p.coeff),p.VC=p.VC, VC=VC),
               real=list(psi=as.data.frame(psi.est), p=as.data.frame(p.est)),
               derived=list(psi_c=as.data.frame(psi_c.est)), gof=gof,
               warnings=list(conv=warn.conv,VC=warn.VC),
               version=list(PRESENCE=version,RPresence=packageVersion("RPresence")))
  class(result)=c("occMod","so")
  return(result)
}
