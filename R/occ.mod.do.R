###   occMod.do1cd.R    ######################################
occMod_DO<-function(psi=call(),psi.cov=data$unitcov,
                      gamma=call(),gamma.cov=data$unitcov,
                      epsilon=call(),epsilon.cov=data$unitcov,
                      p=call(),p.cov=cbind(data$unitcov,data$survcov),
                      theta=NULL, theta.cov=data$unitcov,th0pi=NULL, th0pi.cov=data$unitcov,
                      modname=NULL, paoname=NULL, outfile, model=200,fixed=NULL,initvals=NULL,
                      data=data,conf=0.95,miscopts=''){

#'    Fit dynamic (multi-season) occupancy model using the first parameterisation in PRESENCE.
#'
#'        This is not intended for direct use, but instead the \code{\link{occMod}}
#'        function should be used with \code{type="do.1"}.
  #' @param psi the right-hand side of the formula for the model to fit for occupancy probability in the first season.
  #' @param psi.cov a data frame containing the unit-specific covariates to use for the occupancy component of the model, with number of rows = \code{data$nunits}.
  #' @param gamma the right-hand side of the formula for the model to fit for colonization probability in each season.
  #' @param gamma.cov a data frame containing the unit-specific covariates to use for the colonization component of the model, with number of rows = \code{data$nunits*(data$nseasons-1)}.
  #' @param epsilon the right-hand side of the formula for the model to fit for extinction probability in each season.
  #' @param epsilon.cov a data frame containing the unit-specific covariates to use for the extinction component of the model, with number of rows = \code{data$nunits*(data$nseasons-1)}.
  #' @param p the right-hand side of the formula for the model to fit for detection probability.
  #' @param p.cov a data frame containing the survey-specific covariates to use for the detection component of the model.
  #' @param theta (optional) the right-hand side of the formula for the model to fit for local occupancy parameter (th0,th1). Set equal to NULL for standard dynamic model.
  #' @param theta.cov a data frame containing the survey-specific covariates to use for the local occupancy component of the model.
  #' @param th0pi the right-hand side of the formula for the model to fit for initial local occupancy probability.
  #' @param th0pi.cov a data frame containing the site-specific covariates to use for th0pi.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #' @param model the PRESENCE model code. DO NOT CHANGE.
  #' @param fixed a single-column matrix containing values for real parameters to be fixed at. \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals initial values for the beta parameters at which PRESENCE begins the optimisation. The default values in PRESENCE is 0.
  #' @param data the \code{pao} data object containing the detection data and other information.
  #' @param conf limits for confidence intervals as a proportion (defalut=0.95 for 95\% conf. interval limits)
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param miscopts (see \code{\link{occMod}})
  #' @details Pre-defined covariates:
  #' \itemize{
  #' \item{SURVEY} {categorical covariate indicating survey (for detection parameters)}
  #' \item{SEASON} {categorical covariate indicating season (for detection or col/ext parameters)}
  #' }
  #' @return returns a list of class \code{"occMod"} and \code{"do1"}.
  #'
  #' \code{occMod$beta} contains the objects:
  #'  \item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence in the first year.}
  #'  \item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence in the first year.}
  #'  \item{gamma}{estimated logistic regression coefficients and standard errors for probability of colonization.}
  #'  \item{gamma.VC}{variance-covariance matrix for logistic regression coefficients for probability of colonization.}
  #'  \item{epsilon}{estimated logistic regression coefficients and standard errors for probability of extinction.}
  #'  \item{epsilon.VC}{variance-covariance matrix for logistic regression coefficients for probability of extinction.}
  #'  \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'  \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'  \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'
  #'  \code{occMod$real} contains the objects:
  #'  \item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard errors and limits of 95\% confidence interval. Estimates are provided for the first season (calculated directly from the estimated \code{beta} parameters). Estimates for later seasons are provided as derived parameters. The season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(psi)}).}
  #'  \item{gamma}{estimated probabilities of colonization for each sampling unit, along with standard errors and limits of 95\% confidence interval. The beginning season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(gamma)}).}
  #'  \item{epsilon}{estimated probabilities of extinction for each sampling unit, along with standard errors and limits of 95\% confidence interval. The beginning season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(epsilon)}).}
  #'  \item{phi}{estimated probabilities of persistence (i.e., \code{(1-epsilon)}) for each sampling unit, along with standard errors and limits of 95\% confidence interval. The beginning season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(phi)}).}
  #'  \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits of 95\% confidence interval.}
  #'
  #'  \code{occMod$derived} contains the objects:
  #'  \item{psi}{estimated probabilities of occurrence for each sampling unit for second season onwards, along with standard errors and limits of 95\% confidence interval. The season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(psi)}).}

  #' @author Darryl MacKenzie
  #' @seealso \code{\link{occMod}}

  options("na.action"="na.pass")
  #### create season-specific covariate table for gamma and epsilon
  if(is.null(psi.cov)) { psi.cov=data.frame(rep(1,data$nunits)) }
  SEASON=SRVY=NULL;
  for (i in 1:data$nseasons) {
    SEASON=c(SEASON,rep(i,data$nunits*data$nsurveyseason[i]))
    SRVY  =c(SRVY  ,rep(1:data$nsurveyseason[i],each=data$nunits))
  }
  SEASON=as.factor(SEASON); SRVY=as.factor(SRVY)
  #contrasts(SEASON)=contr.treatment(data$nseasons,base=data$nseasons)
  if (is.null(p.cov)) p.cov=data.frame(SEASON=SEASON,SRVY=SRVY) else { p.cov$SEASON=SEASON; p.cov$SRVY=SRVY }
            # delete last season for gam/eps covs
  SEASON=as.factor(rep(1:(data$nseasons-1),each=data$nunits));
  #contrasts(SEASON)=contr.treatment(data$nseasons-1,base=data$nseasons-1)
  if(is.null(gamma.cov))
    gamma.cov=  data.frame(SEASON=SEASON)
  else {
      s=names(gamma.cov); gamma.cov=data.frame(gamma.cov[1:(data$nunits*(data$nseasons-1)),]);
      names(gamma.cov)=s; gamma.cov$SEASON=SEASON
  }
  if(is.null(epsilon.cov))
    epsilon.cov=data.frame(SEASON=SEASON)
  else {
    s=names(epsilon.cov); epsilon.cov=data.frame(epsilon.cov[1:(data$nunits*(data$nseasons-1)),]);
    names(epsilon.cov)=s; epsilon.cov$SEASON=SEASON
  }
  if(is.null(th0pi.cov)) { th0pi.cov=data.frame(rep(1,data$nunits)) }
  ##################
  psi_mat=    model.matrix(psi,psi.cov)
  gamma_mat=  model.matrix(gamma,gamma.cov[1:(data$nunits*(data$nseasons-1)),])
  epsilon_mat=model.matrix(epsilon,epsilon.cov[1:(data$nunits*(data$nseasons-1)),])
  p_mat=      model.matrix(p,p.cov)
  th_mat=th0pi_mat=NULL
  if (!is.null(theta)) {
    thSRVY=rep(rep(1:data$nsurveys,each=data$nunits),2); thPRIME=rep(1:2,each=data$nunits*data$nsurveys)
    thetacov1=data.frame(SURVEY=as.factor(thSRVY),PRIME=as.factor(thPRIME))
    if (is.null(theta.cov)) theta.cov=data.frame(thetacov1) else theta.cov=data.frame(theta.cov,thetacov1)
    th_mat = model.matrix(theta,theta.cov)
    if (is.null(th0pi.cov)) {
      th0piSEASN=rep(1:data$nseasons,each=data$nunits)
      th0pi.cov=data.frame(SEASON=th0piSEASN)
    }
    th0pi_mat = model.matrix(th0pi,th0pi.cov)
  }

  ## create and output temporary pao file
  ### bind each block of psi covariates to unitcov
  temp.cov=psi_mat
  for(ii in 1:(data$nseasons-1)){
    temp.cov=cbind(temp.cov, gamma_mat[  (ii-1)*data$nunits+1:data$nunits,],
                              epsilon_mat[(ii-1)*data$nunits+1:data$nunits,])
  }

  s=c(paste0("psi.",colnames(psi_mat)),
                        paste0(rep(c(paste0("gamma.",  colnames(gamma_mat)),
                                     paste0("epsilon.",colnames(epsilon_mat))),data$nseasons-1),
                        rep(1:(data$nseasons-1),each=(ncol(gamma_mat)+ncol(epsilon_mat)))))

  if (!is.null(theta)) {
    temp.cov=cbind(temp.cov,th0pi_mat)
    s=c(s,paste0("th0pi.",colnames(th0pi_mat)))
  }
  colnames(temp.cov)=gsub("(Intercept)","int",s,fixed=T)
  temp.pao=data
  temp.pao$unitcov=as.data.frame(temp.cov); temp.pao$nunitcov=ncol(temp.cov)

  temp.cov=p_mat; l=nrow(p_mat); s=paste0('p.', colnames(p_mat))
  if (!is.null(theta)) {
    for (i in 0:1) temp.cov=cbind(temp.cov,th_mat[i*l+1:l,])
    s=c(s,paste0('th0.',colnames(th_mat)),paste0('th1.',colnames(th_mat)))
  }
  colnames(temp.cov)=gsub("(Intercept)","int",s,fixed=T)
  temp.pao$survcov=as.data.frame(temp.cov);    temp.pao$nsurvcov=ncol(temp.cov)
  temp.pao$paoname=ifelse(is.null(paoname),"paodata.pao",paoname)

  #writePao(temp.pao)

  ## create design matrices and run PRESENCE
  psi.dm=matrix(colnames(temp.pao$unitcov)[1:ncol(psi_mat)],1,ncol(psi_mat),byrow=TRUE)
  rownames(psi.dm)=c("psi1");   colnames(psi.dm)=paste("a",1:ncol(psi_mat),sep="")

  index=grep("gamma.",colnames(temp.pao$unitcov),fixed=TRUE)
  gamma.dm=matrix(colnames(temp.pao$unitcov)[index], (temp.pao$nseasons-1),ncol(gamma_mat),byrow=TRUE)
  rownames(gamma.dm)=paste("gamma",1:(temp.pao$nseasons-1),sep="")
  colnames(gamma.dm)=paste("b",1:ncol(gamma_mat),sep="")

  index=grep("epsilon.",colnames(temp.pao$unitcov),fixed=TRUE)
  epsilon.dm=matrix(colnames(temp.pao$unitcov)[index], (temp.pao$nseasons-1),ncol(epsilon_mat),byrow=TRUE)
  rownames(epsilon.dm)=paste("epsilon",1:(temp.pao$nseasons-1),sep="")
  colnames(epsilon.dm)=paste("c",1:ncol(epsilon_mat),sep="")

  index=grep("p.",colnames(temp.pao$survcov),fixed=TRUE)
  p.dm=matrix(rep(colnames(temp.pao$survcov)[index],temp.pao$nsurveys),nrow=temp.pao$nsurveys,  byrow=T)
  s=colnames(temp.pao$survcov); nsrvys=temp.pao$nsurveys; theta.dm=th0pi.dm=NULL
  if (!is.null(theta)) {
    index=grep("th0pi.",colnames(temp.pao$unitcov))
    th0pi.dm=matrix(colnames(temp.pao$unitcov)[index],nrow=temp.pao$nseasons)
    index=grep("th0.",s); theta.dm=matrix(rep(s[index],nsrvys),nrow=nsrvys,byrow=T)
    index=grep("th1.",s); theta.dm=rbind(theta.dm,matrix(rep(s[index],nsrvys),nrow=nsrvys,byrow=T))
    rownames(theta.dm)=c(paste0("th0(",1:nsrvys,')'),paste0('th1(',1:nsrvys,')'))
    colnames(theta.dm)=paste0('a',1:ncol(theta.dm)+1)
  }
  psi.dm=simplify_dm(psi.dm,temp.pao);       p.dm=simplify_dm(p.dm,temp.pao)
  gamma.dm=simplify_dm(gamma.dm,temp.pao);   epsilon.dm=simplify_dm(epsilon.dm,temp.pao)

  if (!is.null(theta)) {
    theta.dm=simplify_dm(theta.dm,temp.pao);  th0pi.dm=simplify_dm(th0pi.dm,temp.pao)
    rownames(th0pi.dm)=paste('th0pi(',1:temp.pao$nseasons,')'); colnames(th0pi.dm)=paste0('e',1:ncol(th0pi.dm))
    psi.dm=diagbind(psi.dm,theta.dm); rownames(psi.dm)=c('psi',rownames(theta.dm))
    colnames(psi.dm)=paste0('a',1:ncol(psi.dm))
  }
  rownames(p.dm)=paste0("P[",
                        unlist(lapply(1:temp.pao$nseasons, function(xx){rep(xx,temp.pao$nsurveyseason[xx])})),
                        "-",
                        unlist(lapply(1:temp.pao$nseasons, function(xx){1:temp.pao$nsurveyseason[xx]})),"]")
  colnames(p.dm)=paste0("d",1:ncol(p.dm))

  if(!is.null(fixed)){
    fixed$idx=match(fixed$param,unlist(lapply(list(psi.dm,gamma.dm,epsilon.dm,p.dm),rownames)))
  }
  rv=runPresence(temp.pao,list(psi.dm,gamma.dm,epsilon.dm,p.dm,th0pi.dm,NULL),
                model=model,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  ################################################################
  #### Extract results from output file
  ################################################################
  v=readLines(outfile)

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  v=v[-1:-i]

  ### get Untransformed coefficients (Beta's)
  i=grep("^A.+psi.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
  j=grep("^B.+ : ",v); gamma.coeff=getbetas(v[j],gamma.dm)
  k=grep("^C.+ : ",v); epsilon.coeff=getbetas(v[k],epsilon.dm)
  l=grep("^D.+ : ",v); p.coeff=getbetas(v[l],p.dm)
  theta.coeff=th0pi.coeff=NULL

  ### get_VC matrix
  ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)
  instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  ii=1:length(i); psi.VC=VC[ii,ii]
  if (is.null(theta)) {
    gamma.VC=VC[length(i)+1:length(j), length(i)+1:length(j)]
    epsilon.VC=VC[length(c(i,j))+1:length(k), length(c(i,j))+1:length(k)]
    p.VC=VC[length(c(i,j,k))+1:length(l), length(c(i,j,k))+1:length(l)]
    ijk=1:length(c(i,j,k))
  } else {
    i1=grep("^A.+th[01].+:",v); getbetas(v[i1],theta.dm) #theta.coeff=as.numeric(gsub(" +\\S+$","",gsub(".+: +","",v[i1])))
    i2=grep("^E.+th0pi.+:",v); th0pi.coeff=as.numeric(gsub(" +\\S+$","",gsub(".+: +","",v[i2])))
    theta.VC=VC[length(c(i))+1:length(i1), length(c(i))+1:length(i1)]
    gamma.VC=VC[length(c(i,i1))+1:length(j), length(c(i,i1))+1:length(j)]
    epsilon.VC=VC[length(c(i,i1,j))+1:length(k), length(c(i,i1,j))+1:length(k)]
    p.VC=VC[length(c(i,i1,j,k))+1:length(l), length(c(i,i1,j,k))+1:length(l)]
    th0pi.VC=VC[length(c(i,i1,j,k,l))+1:length(i2), length(c(i,i1,j,k,l))+1:length(i2)]
    ijk=c(1:length(i),length(c(i,i1))+1:length(c(j,k)))
  }
  psi.est=gamma.est=epsilon.est=p.est=theta.est=th0pi.est=NULL
  VCoutopt=miscopts[1]
  if (! (VCoutopt %in% c("nose","betavc","noreal"))) {  ## return real parameters
    #     ## Calculate real parameter estimates
    psi.est=calc_real(cov=psi_mat,coeff=as.numeric(psi.coeff$est),VC=psi.VC,conf=conf,rownames=paste0('psi_',temp.pao$unitnames))
    s=paste0('gamma',rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits),'_',rep(temp.pao$unitnames,temp.pao$nseasons-1))
    gamma.est=calc_real(cov=gamma_mat,coeff=as.numeric(gamma.coeff$est), VC=as.matrix(gamma.VC),conf=conf,rownames=s)
    s=gsub('gamma','epsilon',s)
    epsilon.est=calc_real(cov=epsilon_mat,coeff=as.numeric(epsilon.coeff$est), VC=as.matrix(epsilon.VC),conf=conf,rownames=s)
    s=paste0('p',rep(1:temp.pao$nsurveys,each=temp.pao$nunits),'_',rep(temp.pao$unitnames,temp.pao$nsurveys))
    p.est=calc_real(cov=p_mat,coeff=as.numeric(p.coeff$est), VC=as.matrix(p.VC),conf=conf,rownames=s)
    theta.est=th0pi.est=NULL
    if (!is.null(theta)) {
      i1=rep(1:temp.pao$nunits,temp.pao$nsurveys)
      i2=rep(1:temp.pao$nsurveys,each=temp.pao$nunits)
      rnames=paste0('th0(',i1,')unit',i2)
      rnames=c(rnames,gsub('th0','th1',rnames))
      theta.est=calc_real(cov=th_mat,coeff=as.numeric(theta.coeff), VC=as.matrix(theta.VC),
                        conf=conf,rownames=rnames)
      th0pi.est=calc_real(cov=th0pi_mat,coeff=as.numeric(th0pi.coeff), VC=as.matrix(th0pi.VC),
                        conf=conf,rownames=paste0('th0pi_',temp.pao$unitnames))
    }
    temp.psi=temp.SE=rep(NA,temp.pao$nunits*(temp.pao$nseasons-1))

    for(ii in 1:temp.pao$nunits){                                  #   compute derivitive of real parm wrt betas
      index=seq(ii,by=temp.pao$nunits,length.out=(temp.pao$nseasons-1))
      deriv=diagbind(
              diagbind(matrix(psi_mat[ii,]*psi.est$est[ii]*(1-psi.est$est[ii]),nrow=1),
                       matrix(gamma_mat[index,]*gamma.est$est[index]*(1-gamma.est$est[index]),nrow=temp.pao$nseasons-1)),
              matrix(epsilon_mat[index,]*epsilon.est$est[index]*(1-epsilon.est$est[index]),nrow=temp.pao$nseasons-1))
      temp.VC=deriv%*%VC[ijk,ijk]%*%t(deriv)  #  ikj=indices if psi,gam,eps in VC matrix
      psi.temp=psi.est$est[ii]
      psi.deriv=array(0,dim=c(temp.pao$nseasons,2*(temp.pao$nseasons-1)+1)); psi.deriv[1,1]=1

      for(jj in 1:(temp.pao$nseasons-1)){            #  compute seasonal psi and variance using delta method
        temp.psi[index[jj]]=psi.temp*(1-epsilon.est$est[index[jj]])+(1-psi.temp)*gamma.est$est[index[jj]]
        temp.ind=c(1,2:(2+jj-1),(temp.pao$nseasons+1):(temp.pao$nseasons+1+jj-1))
        psi.deriv[jj+1,temp.ind]=psi.deriv[jj,temp.ind]*((1-epsilon.est$est[index[jj]])-gamma.est$est[index[jj]])
        psi.deriv[jj+1,1+jj]=(1-psi.temp)
        psi.deriv[jj+1,temp.pao$nseasons+jj]= -psi.temp
        psi.temp=temp.psi[index[jj]]
      }
      temp.psi.VC=psi.deriv%*%temp.VC%*%t(psi.deriv); temp.SE[index]=sqrt(diag(temp.psi.VC)[-1])
    }
    logit.psi=qlogis(temp.psi); logit.se=temp.SE/(temp.psi*(1-temp.psi))

    alpha=(1-conf)/2; z=-qnorm(alpha)
    lower=plogis(sapply(z,function(zz) logit.psi-zz*logit.se))
    upper=plogis(sapply(z,function(zz) logit.psi+zz*logit.se))

    temp=data.frame(temp.psi,temp.SE,lower,upper); colnames(temp)=colnames(psi.est)

    psi.est=rbind(psi.est,temp)
    rownames(psi.est)=paste(rep(temp.pao$unitnames,temp.pao$nseasons),
                          rep(1:temp.pao$nseasons,each=temp.pao$nunits),sep="_")
  }  #  end if compute real estimates...
  ##### get gof results
  i=grep('c-hat',v); chat=as.numeric(gsub(' [(].+','',gsub('.+c-hat = ','',v[i])))
  i=grep('Test Statistic',v); tstat=as.numeric(gsub('.+= +','',v[i]))
  i=grep('Lowest',v); tslow=as.numeric(gsub('.+= +','',v[i]))
  i=grep('Highest',v); tshi=as.numeric(gsub('.+= +','',v[i]))
  i=grep('tic >= obs',v); tsprob=as.numeric(gsub('.+= +','',v[i]))

  ##### check for warnings
  warn.conv=check_conv_warn(v); warn.VC=check_VC_warn(v)

  beta=list(psi=psi.coeff,psi.VC=psi.VC,
            gamma=gamma.coeff,gamma.VC=gamma.VC,
            epsilon=epsilon.coeff,epsilon.VC=epsilon.VC,
            p=p.coeff,p.VC=p.VC,
            VC=VC)
  real=list(psi=as.data.frame(psi.est)[(1:temp.pao$nunits),],
            gamma=as.data.frame(gamma.est),
            epsilon=as.data.frame(epsilon.est),
            p=as.data.frame(p.est))
  if (!is.null(theta)) {
    beta[["theta"]]=as.data.frame(theta.coeff)
    beta[["th0pi"]]=as.data.frame(th0pi.coeff)
    real[["theta"]]=as.data.frame(theta.est)
    real[["th0pi"]]=as.data.frame(th0pi.est)
  }
  result=list(modname=modname,
               model=list(psi=psi,gamma=gamma,epsilon=epsilon,p=p),
               dmat=list(psi=psi.dm,gamma=gamma.dm,epsilon=epsilon.dm,p=p.dm,th0pi=th0pi.dm),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike, aic=aic,
               npar=(ncol(psi.dm)+ncol(gamma.dm)+ncol(epsilon.dm)+ncol(p.dm)),

               beta=beta,real=real,derived=list(psi=as.data.frame(psi.est)[-(1:temp.pao$nunits),]),
               gof=list(TS=tstat, Prob=tsprob, chat=chat),
               warnings=list(conv=warn.conv,VC=warn.VC),
               version=list(PRESENCE=version,RPresence=packageVersion("RPresence")))

  class(result)=c("occMod","do1")
  return(result)
}
