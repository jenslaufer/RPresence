#########################################
occMod_DO_ms1<-function(psi=call(),psi.cov=data$unitcov,
                        phi=call(), phi.cov=data$unitcov,
                      p=call(),p.cov=data$survcov,
                      modname=NULL, paoname=NULL, outfile,
                      model=6000,fixed=NULL,initvals=NULL,data,conf,miscopts=''){
#' Fit multi-state dynamic (multi-season) occupancy model using the 1st parameterisation (\eqn{\psi/\phi}).
#'
#' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="do.ms.1"}.
#'
#' @param psi the right-hand side of the formula for the model to fit for occupancy probability.
#' @param psi.cov a data frame containing the unit-specific covariates to use for the occupancy component of the model.
#' @param phi the right-hand side of the formula for the model to fit for occupancy transition probabilities,.
#' @param phi.cov a data frame containing the unit-specific covariates to use for the occupancy component of the model.
#' @param p the right-hand side of the formula for the model to fit for detection probability.
#' @param p.cov a data frame containing the survey-specific covariates to use for the detection component of the model.
#' @param modname (optional) a string containing the model name
#' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #'@param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param model the PRESENCE model code. DO NOT CHANGE.
#' @param fixed a 2-column matrix containing indices and values for real parameters to be fixed at.
#' @param initvals initial values for the beta parameters at which PRESENCE begins the optimisation. Default = 0.
#' @param data the \code{pao} data object containing the detection data and other information.
  #'@param conf level for confidence interval (may be vector valued).
  #' @param miscopts see \code{\link{occMod}}
  #'
#' @note
#' NOTE THAT THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED SO EXACT IMPLEMENTATION MAY CHANGE.
#'
#' psi built-in covariates:
#' \itemize{
#'    \item{STATE   - for model where psi (initial occupancy) depends on occupancy state.}
#' }
#' phi built-in covariates:
#' \itemize{
#'    \item{FROM - for model where psi depends on state species is transitioning from.}
#'    \item{TO   - for model where psi depends on state species is transitioning to.}
#'    \item{FROM_TO - for model where psi depends on "from" state and "to" state (equalivent to "FROM*TO").}
#'    \item{SEASON - for model where psi changes from season to season.}
#' }
#' P built-in covariates:
#' \itemize{
#'    \item{OSTATE - detection depends on observed state of species.}
#'    \item{OTSTATE - detection depends on both observed and true state of species.}
#'    \item{SURVEY  - detection varies by survey.}
#'    \item{SEASON  - detection varies by season (constant within season).}
#' }
#' Note: P21 (prob detect as state=2, when true state=1) must be fixed to 0, or model will not converge
#' Also: setting outfile="modname" causes RPresence to save the output file using the model name
#'                         (with ".out" extension). Otherwise, you can name it something else.
#'                         Or, set outfile=NULL to not save the output file.
#'
#' @return returns a list of class \code{"occMod"} and \code{"do.1"}.
#'
#' \code{occMod$beta} contains the objects:
#'  \item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence in the first year.}
#'  \item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence in the first year.}
#'  \item{phi}{estimated logistic regression coefficients and standard errors for probability of state transition between years.}
#'  \item{phi.VC}{variance-covariance matrix for logistic regression coefficients for probability of transition between years.}
#'  \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
#'  \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
#'  \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
#'
#'  \code{occMod$real} contains the objects:
#'  \item{psi}{estimated probabilities of occurrence for each sampling unit and state, along with standard errors and limits
#'  of 95\% confidence interval.}
#'  \item{phi}{estimated probabilities of transition between occupancy states for each sampling unit and state, along with standard errors and limits
#'  of 95\% confidence interval.}
#'  \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits of 95\% confidence interval.}
#'
#' @author Jim Hines
#' @seealso \code{\link{occMod}}, \code{\link{occMod_DO_ms2}}

  if(!("pao"%in%class(data))){ print("data must be a pao object!!"); return(NULL) }
  options("na.action"="na.pass")
  nstates=max(data$det.data,na.rm=TRUE)+1; ns1=nstates-1; NU=data$nunits; NINT=data$nseasons-1
  #### create season-specific covariate table for psi
  psi_lbl=paste0('psi0',1:ns1,'(0)')
  phi_lbl=paste0('phi',rep(0:ns1,each=(data$nseasons-1)*ns1),rep(1:ns1,each=NINT),paste0('(',1:NINT,')'))
  p_lbl=paste0('p',rep(1:ns1,each=data$nsurveys),rep(1:ns1,each=data$nsurveys*ns1),paste0('(',1:(data$nsurveys),')'))

  STATE=1:2
  psicov=data.frame(STATE=as.factor(rep(STATE,each=data$nunits)))

  SEASON=rep(1:NINT,nstates*ns1)
  FROM=rep(-1+1:nstates,each=ns1*NINT)
  FROM_TO=rep(1:(nstates*ns1),each=NINT)
  TO=rep(rep(-1+2:nstates,each=NINT),nstates)

  phicov<-data.frame(SEASON=as.factor(rep(SEASON,each=data$nunits)),
                     FROM=as.factor(rep(FROM,each=data$nunits)),
                     TO=as.factor(rep(TO,each=data$nunits)),
                     FROM_TO=as.factor(rep(FROM_TO,each=data$nunits)))
  if (!is.null(psi.cov)) psi.cov=data.frame(psi.cov,psicov) else psi.cov=psicov
  if (!is.null(phi.cov)) phi.cov=data.frame(phi.cov,phicov) else phi.cov=phicov
  rm(FROM,FROM_TO,SEASON,TO)
  #  if(!is.null(data$unitcov)) psi.cov<-data.frame(psi.cov,data$unitcov)
  ## combine unit and survey specific covariates

  SURVEY=rep(1:data$nsurveys,ns1^2)
  SEASON=rep(unlist(sapply(1:data$nseasons,function(xx){rep(xx,data$nsurveyseason[xx])})),ns1^2)
  OSTATE=substr(p_lbl,2,2); OTSTATE=substr(p_lbl,2,3)

  pcov<-data.frame(SURVEY=as.factor(rep(SURVEY,each=data$nunits)),
                    SEASON=as.factor(rep(SEASON,each=data$nunits)),
                    OSTATE=as.factor(rep(OSTATE,each=data$nunits)),
                    OTSTATE=as.factor(rep(OTSTATE,each=data$nunits))     )
  rm(SEASON,SURVEY,OSTATE,OTSTATE)
  if (!is.null(p.cov)) p.cov=data.frame(p.cov,pcov) else p.cov=pcov
  ##################
  psi_mat<-model.matrix(psi,psi.cov)
  phi_mat=model.matrix(phi,phi.cov)
  p_mat<-model.matrix(p,p.cov)

  ## create and output temporary pao file
  ### bind each block of psi covariates to unitcov
  temp.unitcov=matrix(psi_mat,nrow=data$nunits)
  colnames(temp.unitcov)=unlist(sapply(colnames(psi_mat),function(x)paste(psi_lbl,x,sep='.')))
  temp2.unitcov=matrix(phi_mat,nrow=data$nunits)
  colnames(temp2.unitcov)=unlist(sapply(colnames(phi_mat),function(x)paste(phi_lbl,x,sep='.')))
  temp.unitcov=cbind(temp.unitcov,temp2.unitcov)
  colnames(temp.unitcov)=gsub("(Intercept)","int",colnames(temp.unitcov),fixed=TRUE)

  temp.pao<-data; temp.pao$unitcov<-temp.unitcov; temp.pao$nunitcov<-ncol(temp.unitcov)

  ### bind each block of p covariates to cov
  temp.survcov<-matrix(p_mat,nrow=data$nunits)
  colnames(temp.survcov)=unlist(sapply(colnames(p_mat),function(x)paste(p_lbl,x,sep='.')))
  colnames(temp.survcov)=gsub(")","",gsub("(","_",gsub("(Intercept)","int",colnames(temp.survcov),fixed=TRUE),fixed=T),fixed=T)

  temp.cov<-data.frame(temp.survcov)

  temp.pao$survcov<-temp.cov
  temp.pao$nsurvcov<-ncol(temp.cov)
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)

  #writePao(temp.pao)

  ## create design matrices and run PRESENCE
  s=colnames(temp.pao$unitcov); i=grep("^psi",s); psi.dm<-matrix(s[i],nrow=ns1,byrow=FALSE)
  rownames(psi.dm)=psi_lbl; colnames(psi.dm)=paste0("a",1:ncol(psi_mat))
  i=grep("^phi",s); phi.dm<-matrix(s[i],nrow=NINT*ns1*nstates,byrow=FALSE)
  rownames(phi.dm)=phi_lbl; colnames(phi.dm)=paste0("a",1:ncol(phi_mat))

  p.dm<-matrix(colnames(temp.pao$survcov),ncol=ncol(p_mat),byrow=FALSE)
  rownames(p.dm)=p_lbl; colnames(p.dm)=paste0("d",1:ncol(p_mat))

  i=unlist(sapply(fixed$param, function(x) which(x==rownames(p.dm)))); if (length(i)>0) p.dm[i,]=0
  j=unlist(sapply(fixed$param, function(x) which(x==rownames(psi.dm)))); if (length(j)>0) psi.dm[j,]=0

  tmp=new_simplify_dm(psi.dm,temp.pao); psi.dm=tmp$dm; psi.nzcol=tmp$nzcol
  tmp=new_simplify_dm(phi.dm,temp.pao); phi.dm=tmp$dm; phi.nzcol=tmp$nzcol
  tmp=new_simplify_dm(p.dm,temp.pao); p.dm=tmp$dm; p.nzcol=tmp$nzcol

  remove_unused_covs <- function(dm,pao) { #  get rid of site covariates not in design matrices
    s1=as.character(colnames(pao$unitcov))
    s2=as.character(dm)
    newunitcov=NULL; for (i in 1:length(s1)) if (s1[i] %in% s2) newunitcov=cbind(newunitcov,pao$unitcov[,i])
    if (is.null(newunitcov)) newunitcov=nunitcov=0 else nunitcov=ncol(newunitcov)
    pao$unitcov=newunitcov; pao$nunitcov=nunitcov
    s1=colnames(pao$survcov)
    newsurvcov=NULL; for (i in 1:length(s1)) if (s1[i] %in% s2) newsurvcov=cbind(newsurvcov,pao$survcov[,i])
    if (is.null(newsurvcov)) newsurvcov=nsurvcov=0 else nsurvcov=ncol(newsurvcov)
    pao$survcov=newsurvcov; pao$nsurvcov=nsurvcov
    return(pao)
  }
  temp.pao=remove_unused_covs(c(psi.dm,phi.dm,p.dm),temp.pao)
  psi.dm=diagbind(psi.dm,phi.dm); colnames(psi.dm)=paste0('a',1:ncol(psi.dm))
  if(!is.null(fixed)) fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm),rownames)))
  if (data$nseasons>1) {
    rv<-runPresence(temp.pao,list(psi.dm,NULL,NULL,NULL,p.dm,NULL),
                              model=model,modname=modname,
                              fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  } else {
    colnames(p.dm)=paste0("b",1:ncol(p_mat))
    rv<-runPresence(temp.pao,list(psi.dm,p.dm,NULL,NULL,NULL,NULL),
                              model=model,modname=modname,
                              fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  }

  ################################################################
  #### Extract results from output file
  ################################################################

  v=readLines(outfile); npar<-(ncol(psi.dm)+ncol(p.dm))

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  v=v[-1:-i]

  ### get Untransformed coefficients (Beta's)
  i=grep("^A.+p[sh]i.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
  j=grep("^E.+ : ",v); p.coeff=getbetas(v[j],p.dm)

  ### get_VC matrix
  ii=grep("^[AE].+ : ",v); names=c(rownames(psi.coeff),rownames(p.coeff)); npar=length(names)
  instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC=as.matrix(VC[1:length(i), 1:length(i)])
  p.VC  =as.matrix(VC[length(i)+1:length(j), length(i)+1:length(j)])

  i1=grep('phi',rownames(psi.coeff)); phi.coeff=psi.coeff[i1,]; phi.VC=psi.VC[i1,i1]
  psi.coeff=psi.coeff[-i1,]; psi.VC=psi.VC[-i1,-i1]

  psi.est=phi.est=p.est=NULL; NU=temp.pao$nunits
  VCoutopt=miscopts[1]
  if (! (VCoutopt %in% c("nose","betavc","noreal"))) {  ## return real parameters
    ## Calc real parameter estimates
    compute_ms_parm <- function(p_mat,est,VC,ind) {
      dlt=.1e-8; ind1=ind[1]; ind2=ind[2]
      x=p_mat %*% est; i1=which(ind1==fixed$idx); i2=which(ind2==fixed$idx);
      if (length(i1)>0) p=c(fixed$value[i1],plogis(x[2])) else if (length(i2)>0) p=c(plogis(x[1]),fixed$value[i2]) else p=exp(x)/(1+sum(exp(x)))
      grd=matrix(0,2,length(est))
      for (j in 1:length(est)) {
        est[j]=est[j]+dlt; x=p_mat %*% est;
        if (length(i1)>0) p1=c(fixed$value[i1],plogis(x[2])) else if (length(i2)>0) p1=c(plogis(x[1]),fixed$value[i2]) else p1=exp(x)/(1+sum(exp(x)))
        est[j]=est[j]-dlt
        grd[,j]<- (p1-p)/dlt
      }
      vc=grd %*% VC %*% t(grd)
      return(c(p,sqrt(diag(vc))))
    }
    new_compute_ms_parm <- function(p_mat,est,VC,ind) {
      dlt=.1e-8; x=p_mat %*% est; i1=NULL
      for (i in 1:length(fixed$idx)) i1=c(i1,which(ind==fixed$idx[i]))
      if (length(i1)>0) { p=exp(x)/(1+sum(exp(x[-i1]))); p[i1]=fixed$value[i1] } else p=exp(x)/(1+sum(exp(x)))
      grd=matrix(0,length(ind),length(est));
      for (j in 1:length(est)) {
        est[j]=est[j]+dlt; x=p_mat %*% est;
        if (length(i1)>0) { p1=exp(x)/(1+sum(exp(x[-i1]))); p1[i1]=fixed$idx[i1] } else p1=exp(x)/(1+sum(exp(x)))
        est[j]=est[j]-dlt
        grd[,j]=(p1-p)/dlt
      }
      vc=grd %*% VC %*% t(grd)
      return(c(p,sqrt(diag(vc))))
    }
    ipar=1:ns1; psi=matrix(NA,ns1*NU,2)
    psix=t(sapply(1:NU,function(i) new_compute_ms_parm(psi_mat[(ipar-1)*NU+i,psi.nzcol],psi.coeff$est,psi.VC,ipar)))
    for (l in 1:ns1) psi[(l-1)*NU+1:NU,]=psix[,c(l,ns1+l)]
    psi=cbind(psi,sapply(1:nrow(psi),function(i) psi[i,1]-1.96*psi[i,2]))
    psi=cbind(psi,sapply(1:nrow(psi),function(i) psi[i,1]+1.96*psi[i,2]))
    psilbl=rep(rownames(psi.dm)[grep('psi',rownames(psi.dm))],each=NU)
    colnames(psi)=c('est','se','lower_95','upper_95'); rownames(psi)=psilbl
    psi[psi<0]=0; psi[psi>1]=1; psi.est=psi;

    phi.est=matrix(NA,nstates*ns1*NINT*NU,2)
    for (k in 0:ns1) {
      for (yr in 1:NINT) {
        ipar=k*NINT*ns1+(1:ns1-1)*NINT+yr
        phix=t(sapply(1:NU, function(i)
          new_compute_ms_parm(phi_mat[(ipar-1)*NU+i,phi.nzcol],phi.coeff$est,phi.VC,ipar+2)
          ))
        for (l in 1:ns1) phi.est[k*ns1*NINT*NU+(l-1)*NINT*NU+(yr-1)*NU+1:NU,]=phix[,c(l,ns1+l)]
      }
    }
    phi.est=cbind(phi.est,sapply(1:nrow(phi.est),function(i) phi.est[i,1]-1.96*phi.est[i,2]))
    phi.est=cbind(phi.est,sapply(1:nrow(phi.est),function(i) phi.est[i,1]+1.96*phi.est[i,2]))
    phi.est[phi.est<0]=0; phi.est[phi.est>1]=1;
    philbl=rep(rownames(psi.dm)[grep('phi',rownames(psi.dm))],each=NU)
    philbl=paste0(philbl,rep(temp.pao$unitnames,length(rownames(phi.dm))))
    rownames(phi.est)=philbl; colnames(phi.est)=colnames(psi.est)

    nsrvys=temp.pao$nsurveys; p.est=matrix(NA,ns1*ns1*nsrvys*NU,2);
    k=which(rowSums(p.dm!="0")==0)
    if (length(k)>0) for (i in 1:length(k)) p_mat[(k[i]-1)*NU+1:NU,]=0
    for (k in 1:2) {
      p1=p2=NULL
      for (yr in 1:nsrvys) {
        ipar=(k-1)*nsrvys*ns1+(1:ns1-1)*nsrvys+yr
        px=t(sapply(1:NU,function(i)
          new_compute_ms_parm(p_mat[(ipar-1)*NU+i,p.nzcol],p.coeff$est,p.VC,ipar+nrow(psi.dm))
          ))
        for (l in 1:ns1) p.est[(k-1)*ns1*nsrvys*NU+(l-1)*nsrvys*NU+(yr-1)*NU+1:NU,]=px[,c(l,ns1+l)]
      }
    }
  }
  p.est=cbind(p.est,sapply(1:nrow(p.est),function(i) p.est[i,1]-1.96*p.est[i,2]))
  p.est=cbind(p.est,sapply(1:nrow(p.est),function(i) p.est[i,1]+1.96*p.est[i,2]))
  plbl=rep(rownames(p.dm),each=NU)
  plbl=paste0(plbl,rep(temp.pao$unitnames,length(rownames(p.dm))))
  rownames(p.est)=plbl; colnames(p.est)=colnames(psi.est)
  ##### check for warnings
  warn.conv<-check_conv_warn(v); warn.VC<-check_VC_warn(v)

  result<-list(modname=modname,
               model=list(psi=psi,p=p),dmat=list(psi=psi.dm,p=p.dm),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike,npar=npar, aic=aic,
               beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,phi=as.data.frame(phi.coeff),phi.VC=phi.VC,
                         p=as.data.frame(p.coeff),p.VC=p.VC,VC=VC),
               real=list(psi=as.data.frame(psi.est),phi=as.data.frame(phi.est),p=as.data.frame(p.est)),
               warnings=list(conv=warn.conv,VC=warn.VC))

  class(result)<-c("occMod","doMs1")
  return(result)
}
