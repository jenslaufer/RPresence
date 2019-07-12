occMod_DO_fp<-function(psi=call(),psi.cov=data$unitcov,
                       gamma=call(),gamma.cov=data$unitcov,
                       epsilon=call(),epsilon.cov=data$unitcov,
                        p11=call(),p11.cov=cbind(data$unitcov,data$survcov),
                        p10=call(),p10.cov=cbind(data$unitcov,data$survcov),
                        b=call(),b.cov=cbind(data$unitcov,data$survcov),
                        modname=NULL,paoname=NULL,outfile,
                     model=100,fixed=NULL,initvals=NULL,data,conf,miscopts=''){
  #' Fit static occupancy/single season model
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="do.fp"}.
  #'
  #'@param psi the right-hand side of the formula for the model to fit for occupancy probability.
  #'@param psi.cov a data frame containing the unit-specific covariates to use for the occupancy
  #'component of the model.
  #' @param gamma the right-hand side of the formula for the model to fit for colonization probability in each season.
  #' @param gamma.cov a data frame containing the unit-specific covariates to use for the colonization component of the model, with number of rows = \code{data$nunits*(data$nseasons-1)}.
  #' @param epsilon the right-hand side of the formula for the model to fit for extinction probability in each season.
  #' @param epsilon.cov a data frame containing the unit-specific covariates to use for the extinction component of the model, with number of rows = \code{data$nunits*(data$nseasons-1)}.
  #'@param p11 the right-hand side of the formula for the model to fit for "sure" detection probability.
  #'@param p11.cov a data frame containing the survey-specific covariates to use for the detection
  #'component of the model.
  #'@param p10 the right-hand side of the formula for the model to fit for false-positive detection probability.
  #'@param p10.cov a data frame containing the survey-specific covariates to use for the detection
  #'component of the model.
  #'@param b the right-hand side of the formula for the model to fit for Pr(detection is "sure")
  #'@param b.cov a data frame containing the survey-specific covariates to use for the b
  #'@param modname (optional) a string containing the model name
  #'@param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #'@param model the PRESENCE model code. DO NOT CHANGE.
  #'@param fixed a single-column matrix containing values for real parameters to be fixed at.
  #'\code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #'@param initvals initial values for the beta parameters at which PRESENCE begins the optimisation.
  #'The default values in PRESENCE is 0.
  #'@param data the \code{pao} data object containing the detection data and other information.
  #'@param conf level for confidence interval (may be vector valued).
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param miscopts see \code{\link{occMod}}
  #'
  #'@return returns a list of class \code{"occMod"} and \code{"soFp"}.
  #'
  #'\code{occMod$beta} contains the objects:
  #'\item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence.}
  #'\item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence.}
  #'\item{gamma}{estimated logistic regression coefficients and standard errors for probability of colonization.}
  #'\item{gamma.VC}{variance-covariance matrix for logistic regression coefficients for probability of colonization.}
  #'\item{epsilon}{estimated logistic regression coefficients and standard errors for probability of extinction.}
  #'\item{epsilon.VC}{variance-covariance matrix for logistic regression coefficients for probability of extinction.}
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

  SEASON=SRVY=NULL;
  for (i in 1:data$nseasons) {
    SEASON=c(SEASON,rep(i,data$nunits*data$nsurveyseason[i]))
    SRVY  =c(SRVY  ,rep(1:data$nsurveyseason[i],each=data$nunits))
  }
  SEASON=as.factor(SEASON); SRVY=as.factor(SRVY)
  #contrasts(SEASON)=contr.treatment(data$nseasons,base=data$nseasons)
  if (is.null(p11.cov)) p11.cov=data.frame(SEASON=SEASON,SRVY=SRVY) else { p11.cov$SEASON=SEASON; p11.cov$SRVY=SRVY }
  # delete last season for gam/eps covs
  SEASON=as.factor(rep(1:(data$nseasons-1),each=data$nunits));
  #contrasts(SEASON)=contr.treatment(data$nseasons-1,base=data$nseasons-1)
  if(is.null(gamma.cov)) gamma.cov=data.frame(SEASON=SEASON) else gamma.cov$SEASON=SEASON
  if(is.null(epsilon.cov)) epsilon.cov=data.frame(SEASON=SEASON) else epsilon.cov$SEASON=SEASON

  if (is.null(psi.cov))
    psi_mat<-model.matrix(psi,data.frame(rep(1,data$nunits))) else
    psi_mat<-model.matrix(psi,psi.cov)
  colnames(psi_mat)=paste0("psi.",gsub("(Intercept)","int",colnames(psi_mat),fixed=T))
  if (is.null(p11.cov))
    p11_mat<-model.matrix(p11,data.frame(rep(1,data$nunits*data$nsurveys))) else
    p11_mat<-model.matrix(p11,p11.cov)
  colnames(p11_mat)=paste0("p11.",gsub("(Intercept)","int",colnames(p11_mat),fixed=T))
  if (is.null(p10.cov))
    p10_mat<-model.matrix(p10,data.frame(rep(1,data$nunits*data$nsurveys))) else
    p10_mat<-model.matrix(p10,p11.cov)
  colnames(p10_mat)=paste0("p10.",gsub("(Intercept)","int",colnames(p10_mat),fixed=T))
  if (is.null(b.cov))
    b_mat<-model.matrix(b,data.frame(rep(1,data$nunits*data$nsurveys))) else
    b_mat<-model.matrix(b,b.cov)
  colnames(b_mat)=paste0("b.",gsub("(Intercept)","int",colnames(b_mat),fixed=T))
  if (is.null(gamma.cov))
    gamma_mat<-model.matrix(gamma,data.frame(rep(1,(data$nunits*data$nseasons-1)))) else
      try(gamma_mat <-  model.matrix(gamma,gamma.cov))
  if (is.null(epsilon.cov))
    epsilon_mat<-model.matrix(epsilon,data.frame(rep(1,(data$nunits*data$nseasons-1)))) else
      epsilon_mat=  model.matrix(epsilon,epsilon.cov)

  ## create and output temporary pao file
  temp.pao<-data

  temp.cov=psi_mat
  for(ii in 1:(data$nseasons-1)){
    temp.cov=cbind(temp.cov, gamma_mat[  (ii-1)*data$nunits+1:data$nunits,],
                   epsilon_mat[(ii-1)*data$nunits+1:data$nunits,])
  }
  s=c(colnames(psi_mat),
      paste0(rep(c(paste0("gamma.",  colnames(gamma_mat)),
                   paste0("epsilon.",colnames(epsilon_mat))),data$nseasons-1),
             rep(1:(data$nseasons-1),each=(ncol(gamma_mat)+ncol(epsilon_mat)))))
  colnames(temp.cov)=gsub("(Intercept)","int",s,fixed=T)

  temp.pao$nunitcov<-ncol(temp.cov)
  temp.pao$unitcov<-as.data.frame(temp.cov)
  temp.pao$nsurvcov<-ncol(p11_mat)+ncol(p10_mat)+ncol(b_mat)
  temp.pao$survcov<-as.data.frame(cbind(p11_mat,p10_mat,b_mat))
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)

  #writePao(temp.pao)

  ## create design matrices file
  v=colnames(temp.pao$unitcov); i=grep("^psi",v)
  psi.dm<-matrix(v[i],1,length(i))
  psi.dm.col<-ncol(psi.dm)
  rownames(psi.dm)="psi"; colnames(psi.dm)<-paste0("a",1:ncol(psi.dm))

  index=grep("gamma.",colnames(temp.pao$unitcov),fixed=TRUE)
  gamma.dm=matrix(colnames(temp.pao$unitcov)[index], (temp.pao$nseasons-1),ncol(gamma_mat),byrow=TRUE)
  rownames(gamma.dm)=paste0("gamma",1:(temp.pao$nseasons-1))
  colnames(gamma.dm)=paste0("b",1:ncol(gamma_mat))

  index=grep("epsilon.",colnames(temp.pao$unitcov),fixed=TRUE)
  epsilon.dm=matrix(colnames(temp.pao$unitcov)[index], (temp.pao$nseasons-1),ncol(epsilon_mat),byrow=TRUE)
  rownames(epsilon.dm)=paste0("epsilon",1:(temp.pao$nseasons-1))
  colnames(epsilon.dm)=paste0("c",1:ncol(epsilon_mat))

  v=colnames(temp.pao$survcov);
  i=grep("^p11[.]",v)
  p11.dm<-matrix(rep(v[i],temp.pao$nsurveys),nrow=temp.pao$nsurveys,ncol=length(i), byrow=TRUE)
  i=grep("^p10[.]",v)
  p10.dm<-matrix(rep(v[i],temp.pao$nsurveys),nrow=temp.pao$nsurveys,ncol=length(i), byrow=TRUE)
  i=grep("^b[.]",v)
  b.dm<-matrix(rep(v[i],temp.pao$nsurveys),nrow=temp.pao$nsurveys,ncol=length(i), byrow=TRUE)
  p.dm=diagbind(diagbind(p11.dm,p10.dm),b.dm)
  p.dm.col<-ncol(p.dm)
  v<-paste0(rep(c("p11(","p10(","b("),each=temp.pao$nsurveys),1:temp.pao$nsurveys,")")
  rownames(p.dm)=v; colnames(p.dm)<-paste0("b",1:ncol(p.dm))
  npar<-psi.dm.col+p.dm.col
  psi.dm=simplify_dm(psi.dm,temp.pao); gamma.dm=simplify_dm(gamma.dm,temp.pao);
  epsilon.dm=simplify_dm(epsilon.dm,temp.pao); p.dm=simplify_dm(p.dm,temp.pao)

  if(!is.null(fixed)){
    fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,gamma.dm,epsilon.dm,p.dm),rownames)))
  }

  rv<-runPresence(temp.pao,list(psi.dm,gamma.dm,epsilon.dm,p.dm,NULL,NULL),
                  model,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)

  ### extract results from output file
  v=readLines(outfile)

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  if (length(i)>0) {
    v=v[-1:-i]

    ### get Untransformed coefficients (Beta's)
    i=grep("^A.+psi.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
    j=grep("^B.+ : ",v); gam.coeff=getbetas(v[j],gamma.dm)
    k=grep("^C.+ : ",v); eps.coeff=getbetas(v[k],epsilon.dm)
    l=grep("^D.+ : ",v); p.coeff=getbetas(v[l],p.dm)

    ### get_VC matrix
    ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)
    instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
    VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

    psi.VC<-VC[1:length(i), 1:length(i)];
    gam.VC=VC[length(i)+1:length(j), length(i)+1:length(j)]
    eps.VC=VC[length(i)+length(j)+1:length(k), length(i)+length(j)+1:length(k)]
    p.VC=VC[length(i)+length(j)+length(k)+1:length(l), length(i)+length(j)+length(k)+1:length(l)]

    unitnames=rownames(temp.pao$det.data); surveynames=1:temp.pao$nsurveys
    ## Get real parameter estimates from output

    psi.est=p.est=gam.est=eps.est=psi_c.est=NULL

    #  if (!noReal){  ## return real parameters
    ## Calc real parameter estimates
    psi.est=calc_real(cov=psi_mat,coeff=psi.coeff$est,VC=psi.VC,conf=conf,rownames=paste0('psi_',temp.pao$unitnames))
    s=paste0('gamma',rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits),'_',rep(temp.pao$unitnames,temp.pao$nseasons-1))
    gam.est=calc_real(cov=gamma_mat,coeff=gam.coeff$est,VC=gam.VC,conf=conf,rownames=s)
    s=gsub('gamma','epsilon',s)
    eps.est=calc_real(cov=epsilon_mat,coeff=eps.coeff$est,VC=eps.VC,conf=conf,rownames=s)

    s=paste0('p11(',rep(1:temp.pao$nsurveys,each=temp.pao$nunits),')_',rep(temp.pao$unitnames,temp.pao$nsurveys))
    i=grep("p11",rownames(p.coeff)); p11.VC=p.VC[i,i];
    p11.est<-calc_real(cov=p11_mat,coeff=p.coeff$est[i],VC=as.matrix(p11.VC),conf=conf,rownames=s)
    i=grep("p10",rownames(p.coeff)); p10.VC=p.VC[i,i];
    p10.est<-calc_real(cov=p10_mat,coeff=p.coeff$est[i],VC=as.matrix(p10.VC),conf=conf,rownames=gsub("^p11","p10",s))
    i=grep("^b",rownames(p.coeff)); b.VC=p.VC[i,i]
    b.est<-calc_real(cov=b_mat,coeff=p.coeff$est[i],VC=as.matrix(b.VC),conf=conf,rownames=gsub("^p11","b",s))
  }
    ##### check for warnings
  warn.conv<-check_conv_warn(v)
  warn.VC<-check_VC_warn(v)

  result<-list(modname=modname,
               model=list(psi=psi,p11=p11,p10=p10,b=b),dmat=list(psi=psi.dm,p=p.dm),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike,
               npar=npar, aic=aic,
               beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                         gamma=as.data.frame(gam.coeff),gamma.VC=gam.VC,
                         epsilon=as.data.frame(eps.coeff),epsilon.VC=eps.VC,
                         p=as.data.frame(p.coeff),p.VC=p.VC, VC=VC),
               real=list(psi=as.data.frame(psi.est),
                         gamma=as.data.frame(gam.est),epsilon=as.data.frame(eps.est),
                         p11=as.data.frame(p11.est),
                         p10=as.data.frame(p10.est),
                         b=as.data.frame(b.est)),
               warnings=list(conv=warn.conv,VC=warn.VC),
               version=list(PRESENCE=version,RPresence=packageVersion("RPresence")))
  class(result)<-c("occMod","soFp")
  return(result)
}
