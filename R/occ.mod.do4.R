occMod_DO4<-function(psi=call(),psi.cov=data$unitcov,
                      p=call(),p.cov=cbind(data$unitcov,data$survcov),
                      modname=NULL, paoname=NULL, outfile, model=240,  fixed=NULL,initvals=NULL,data,conf,miscopts=''){
  #' Fit dynamic (multi-season) occupancy model using the fourth parameterisation
  #' (random changes in occupancy, epsilon=1-gamma) in PRESENCE.
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should
  #' be used with \code{type="do.4"}.
  #'
  #' @param psi the right-hand side of the formula for the model to fit for occupancy probability
  #' in each season.
  #' @param psi.cov a data frame containing the unit-specific covariates to use for the occupancy
  #' component of the model, with number of rows = \code{data$nunits} or \code{data$nunits*data$nseasons}.
  #' If the shorter version of the data frame is supplied the rows are recycled to the longer length.
  #' @param p the right-hand side of the formula for the model to fit for detection probability.
  #' @param p.cov a data frame containing the survey-specific covariates to use for the detection component
  #' of the model.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #' @param fixed a single-column matrix containing values for real parameters to be fixed at.
  #' \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals initial values for the beta parameters at which PRESENCE begins the optimisation.
  #' The default values in PRESENCE is 0.
  #' @param data the \code{pao} data object containing the detection data and other information.
  #' @param conf limits for confidence intervals as a proportion (defalut=0.95 for 95\% conf. interval limits)
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param model The PRESENCE model code. DO NOT CHANGE.
  #' @param miscopts see \code{\link{occMod}}
  #'
  #' @return list of class \code{"occMod"} and \code{"do4"}.
  #'
  #' \code{occMod$beta} contains the objects:
  #'
  #'  \item{psi}{estimated logistic regression coefficients and standard errors for probabilities of occurrence in each season.}
  #'  \item{psi.VC}{variance-covariance matrix for psi.}
  #'  \item{p}{estimated logistic regression coefficients and standard errors for probabilities of detection in each survey.}
  #'  \item{p.VC}{variance-covariance matrix for p.}
  #'  \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'
  #'\code{occMod$real} contains the objects:
  #'  \item{psi}{estimated probabilities of occurrence for each sampling unit, in each season,
  #'  along with standard errors and limits of 95\% confidence interval. The season for which an
  #'  estimate applies to can be identified from the final number of the row names
  #'  (\code{rownames(psi)}).}
  #'  \item{p}{estimated probabilities of detection for each survey, along with standard errors
  #'  and limits of 95\% confidence intervals.}
  #'
  #'\code{occMod$derived} contains the objects:
  #'  \item{gamma}{estimated probabilities of colonization for each sampling unit, along with
  #'  standard errors and limits of 95\% confidence intervals. The beginning season for which an
  #'  estimate applies to can be identified from the final number of the row names
  #'  (\code{rownames(gamma)}).}
  #'  \item{epsilon}{estimated probabilities of extinction for each sampling unit, along with
  #'  standard errors and limits of 95\% confidence intervals. The beginning season for which an
  #'  estimate applies to can be identified from the final number of the row names
  #'  (\code{rownames(epsilon)}). Under this parameterisation, extinction probabilities have the
  #'  same value as 1-colonization probabilties}
  #'
  #'Note: For this parameterization, \eqn{\gamma_t = \psi_{t+1}}{gamma(t)=psi(t+1)}, due to the
  #'constraint, \eqn{\epsilon=1-\gamma}{epsilon=1-gamma}
  #'            \preformatted{
  #'            ie., psi(t+1)=psi(t)*(1-epsilon(t))+(1-psi(t))*gamma(t)
  #'                 psi(t+1)=psi(t)*gamma(t)+(1-psi(t))*gamma(t)
  #'                 psi(t+1)=psi(t)*gamma(t)+gamma(t)-psi(t)*gamma(t)
  #'                 psi(t+1)=gamma(t)}
  #'
  #' @author Darryl MacKenzie
  #'
  #' @seealso \code{\link{occMod}}
  #'
  options("na.action"="na.pass")
  nunits=data$nunits; nseasons=data$nseasons; SEASN.cov=as.factor(rep(1:nseasons,each=nunits))
  #### create season-specific covariate table for psi

  if(is.null(psi.cov)) psi.cov=data.frame(INT=rep(1,nunits),SEASON=SEASN.cov) else {
    if (nrow(psi.cov)>(nunits*nseasons)) psi.cov=psi.cov[1:(nunits*nseasons),]
    psi.cov=data.frame(psi.cov,INT=rep(1,nunits),SEASON=SEASN.cov)
  }

  ## combine unit and survey specific covariates

  temp.cov<-data.frame();
  for (ii in 1:nseasons) for (jj in 1:data$nsurveyseason[ii]) temp.cov<-rbind(temp.cov,psi.cov[(ii-1)*nunits+1:nunits,])
  if (is.null(p.cov)) p.cov<-data.frame(temp.cov) else p.cov<-data.frame(p.cov,temp.cov)

  ##################
  psi_mat<-model.matrix(psi,psi.cov); p_mat<-model.matrix(p,p.cov)
  colnames(psi_mat)=gsub('[(]Intercept[)]','int',colnames(psi_mat))
  colnames(p_mat)=gsub('[(]Intercept[)]','int',colnames(p_mat))

  ## create and output temporary pao file
  ### bind each block of psi covariates to unitcov
  temp.cov=matrix(psi_mat,nrow=nunits)
  s=NULL; for (i in 1:ncol(psi_mat)) s=c(s,paste0('psi',1:nseasons,'.',colnames(psi_mat)[i]))
  colnames(temp.cov)=s

  temp.pao<-data
  temp.pao$unitcov<-as.data.frame(temp.cov); temp.pao$nunitcov<-ncol(temp.cov)
  temp.pao$survcov<-as.data.frame(p_mat);   temp.pao$nsurvcov<-ncol(p_mat)
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)

  colnames(temp.pao$survcov)<-paste0("p.",colnames(temp.pao$survcov))

  #writePao(temp.pao)

  ## create design matrices and run PRESENCE
  psi.dm=matrix(s,nrow=nseasons);
  rownames(psi.dm)=paste0("gam",1:nseasons-1); rownames(psi.dm)[1]="psi1"
  colnames(psi.dm)<-paste0("a",1:ncol(psi.dm))
  p.dm<-matrix(rep(colnames(temp.pao$survcov),temp.pao$nsurveys),nrow=temp.pao$nsurveys,byrow=T)
  s=NULL; for (i in 1:nseasons) for (j in 1:data$nsurveyseason[i]) s=c(s,paste0('P[',i,'-',j,']'))
  rownames(p.dm)<-s; colnames(p.dm)<-paste0("d",1:temp.pao$nsurvcov)
  psi.dm=simplify_dm(psi.dm,temp.pao)
  p.dm=simplify_dm(p.dm,temp.pao)
  if(!is.null(fixed)) fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,NULL,NULL,p.dm),rownames)))
  rv<-runPresence(temp.pao,list(psi.dm,NULL,NULL,p.dm,NULL,NULL),
                            model=model,modname=modname,
                            fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)

  ################################################################
  #### Extract results from output file
  ################################################################

  v=readLines(outfile); l=length(v);  npar<-ncol(psi_mat)+ncol(p_mat)

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))

  ### get_coefficients
  v=v[-1:-i]; vv=v[grep(" : ",v)][1:npar]

  ### get Untransformed coefficients (Beta's)
  i=grep("^A.+psi.+ : ",vv); psi.coeff=getbetas(vv[i],psi.dm)
  l=grep("^D.+ : ",vv); p.coeff=getbetas(vv[l],p.dm)

  ### get_VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  names<-c(rownames(psi.coeff),rownames(p.coeff))
  VC<-get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  ii=1:length(i); psi.VC<-VC[ii,ii]
  jj=length(i)+1:length(l); p.VC<-VC[jj,jj]
  VCoutopt=miscopts[1]
  if (! (VCoutopt %in% c("nose","betavc","noreal"))) {  ## return real parameters
    ### calc_real estimates
    psi.est<-calc_real(cov=psi_mat,coeff=psi.coeff[,"est"],VC=as.matrix(psi.VC),conf=conf,
                     rownames=paste0(rep(temp.pao$unitnames,nseasons),"_",
                                     rep(1:nseasons,each=temp.pao$nunits)))
    p.est<-calc_real(cov=p_mat,coeff=p.coeff[,"est"],VC=as.matrix(p.VC),conf=conf,
                   rownames=paste0(rep(temp.pao$unitnames,temp.pao$nsurveys),"_",
                                  rep(temp.pao$surveynames,each=temp.pao$nunits)))
    ### derived estimates
    gam.est<-psi.est[-1:-nunits,]; rownames(gam.est)=rownames(psi.est)[1:((nseasons-1)*nunits)]
    epsilon.est=1-gam.est; epsilon.est[,2]=gam.est[,2];
    #  since eps=1-gam, swap conf. interval limits
    tmp=epsilon.est[,3]; epsilon.est[,3]=epsilon.est[,4]; epsilon.est[,4]=tmp
  }  #  end if compute real parms
  ##### check for warnings
  warn.conv<-check_conv_warn(v)
  warn.VC<-check_VC_warn(v)

  result<-list(modname=modname,
               model=list(psi=psi,gamma=gamma,p=p),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike,
               npar=npar, aic=aic,
               beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                         p=as.data.frame(p.coeff),p.VC=p.VC,VC=VC),
               real=list(psi=as.data.frame(psi.est),p=as.data.frame(p.est)),
               derived=list(gamma=as.data.frame(gam.est),epsilon=as.data.frame(epsilon.est)),
               warnings=list(conv=warn.conv,VC=warn.VC),
               version=list(PRESENCE=version,RPresence=packageVersion("RPresence")))

  class(result)<-c("occMod","do4")
  return(result)

}
