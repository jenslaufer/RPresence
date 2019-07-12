#########################################
occMod_DO_ms2<-function(psi=call(),psi.cov=data$unitcov,
                      r=call(),r.cov=data$unitcov,
                      p=call(),p.cov=data$survcov,
                      delta=call(),delta.cov=data$survcov,
                      modname=NULL, paoname=NULL, outfile,
                      model=6100,fixed=NULL,initvals=NULL,data,conf,miscopts){
#' Fit multi-state dynamic occupancy (multi-season) occupancy model using the second (conditional binomial) parameterisation in PRESENCE.
#'
#' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="do.ms.2"}.
#'
#' @param psi the right-hand side of the formula for the model to fit for occupancy probability, both first year and dynamic occupancy. The terms \code{SEASON}, \code{DYN} and \code{PREV_STATE} can each be used without having been defined in \code{psi.cov}. \code{SEASON} allows for a seasonal effect, \code{DYN} enables dynamic occupancy probabilities to be different from first year occupancy, and \code{PREV_STATE} allows the probability to be different depending on the state in the previous season.
#' @param psi.cov a data frame containing the unit-specific covariates to use for the occupancy component of the model.
#' @param r the right-hand side of the formula for the model to fit for the probability of being in the second state, conditional on the unit being occupied. The terms \code{SEASON}, \code{DYN} and \code{PREV_STATE} can be used as above with the same effect.
#' @param r.cov a data frame containing the unit-specific covariates to use for the probability of being in the second state, conditional on the unit being occupied part of the model.
#' @param p the right-hand side of the formula for the model to fit for detection probability. The terms \code{SEASON} and \code{STATE} can each be used without having been defined in \code{p.cov}. \code{SEASON} allows for a seasonal effect, and \code{STATE} allows the probability of detection to be different depending on the state in the current season.
#' @param p.cov a data frame containing the survey-specific covariates to use for the detection component of the model.
#' @param delta the right-hand side of the formula for the model to fit for the probability of detecting the second state in a survey, conditional on the species being detected in the survey. The terms \code{SEASON} and \code{STATE} can each be used as above to the same effect.
#' @param delta.cov a data frame containing the survey-specific covariates to use for the probability of detecting the second state in a survey, conditional on the species being detected in the survey.
#' @param modname (optional) a string containing the model name
#' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #'@param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param model the PRESENCE model code. DO NOT CHANGE.
#' @param fixed a single-column matrix containing values for real parameters to be fixed at. \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
#' @param initvals initial values for the beta parameters at which PRESENCE begins the optimisation. The default values in PRESENCE is 0.
#' @param data the \code{pao} data object containing the detection data and other information.
  #'@param conf level for confidence interval (may be vector valued).
  #'@param miscopts see \code{\link{occMod}}
#'
#' NOTE THAT THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED SO EXACT IMPLEMENTATION MAY CHANGE.
#' @return returns a list of class \code{"occMod"} and \code{"do.1"}.
#'
#' \code{occMod$beta} contains the objects:
#'  \item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence in the first year.}
#'  \item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence in the first year.}
#'  \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
#'  \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
#'  \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
#'
#'  \code{occMod$real} contains the objects:
#'  \item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard errors and limits of 95\% confidence interval. Estimates are provided for the first season (calculated directly from the estimated \code{beta} parameters). Estimates for later seasons are provided as derived parameters. The season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(psi)}).}
#'  \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits of 95\% confidence interval.}
#'
#'  \code{occMod$derived} contains the objects:
#'  \item{psi}{estimated probabilities of occurrence for each sampling unit for second season onwards, along with standard errors and limits of 95\% confidence interval. The season for which an estimate applies to can be identified from the final number of the row names (\code{rownames(psi)}).}

#' @author Darryl MacKenzie
#' @seealso \code{\link{occMod}}

  if(!("pao"%in%class(data))) { print("data must be a pao object!!");  return(NULL)  }
  options("na.action"="na.pass")

  #### create season-specific covariate table for psi and r
  SEASON=0; if (data$nseasons>1) SEASON=c(SEASON,rep(1:(data$nseasons-1),each=3)); NU=data$nunits
  PREV_STATE=c(-1,rep(0:2,data$nseasons-1))
  DYN=c(0,rep(1,each=3*(data$nseasons-1)))
  psicov=data.frame(SEASON    =as.factor(rep(SEASON,    each=data$nunits)),
                         PREV_STATE=as.factor(rep(PREV_STATE,each=data$nunits)),
                    DYN=as.factor(rep(DYN,each=data$nunits)))
  if (!is.null(psi.cov)) psi.cov=data.frame(psi.cov,psicov) else psi.cov=psicov
  if (!is.null(r.cov)) r.cov=data.frame(r.cov,psicov) else r.cov=psicov
                  ## combine unit and survey specific covariates
  SURVEY=rep(1:data$nsurveys,each=data$nunits)
  SEASON=NULL; for (i in 1:data$nseasons) SEASON=c(SEASON,rep(i,data$nsurveyseason[i]))
  SEASON=rep(SEASON,each=data$nunits)
  STATE=rep(1,data$nunits*data$nsurveys)
  ALLDIFF=rep(1:data$nsurveys,each=data$nunits)
  temp.cov=cbind(SURVEY,SEASON,STATE,ALLDIFF)

  n=nrow(temp.cov)
  p.cov=data.frame(SURVEY=as.factor(rep(SURVEY,2)),SEASON=as.factor(rep(SEASON,2)),
                   STATE=as.factor(c(STATE,STATE+1)),ALLDIFF=as.factor(rep(ALLDIFF,2)),p.cov)

  deltacov=p.cov[1:n,]
  if (!is.null(delta.cov)) delta.cov=data.frame(delta.cov,deltacov) else delta.cov=deltacov

  ##################
  psi_mat<-model.matrix(psi,psi.cov)
  r_mat<-model.matrix(r,r.cov)
  p_mat<-model.matrix(p,p.cov)
  delta_mat<-model.matrix(delta,delta.cov)

  ## create and output temporary pao file
  ### bind each block of psi covariates to unitcov
  temp.unitcov=matrix(c(psi_mat,r_mat),nrow=data$nunits);
  s=NULL
  for (j in 1:ncol(psi_mat)) {
    s=c(s,paste0('psi0.',colnames(psi_mat)[j]))
    if (data$nseasons>1)
      for (i in 2:data$nseasons) s=c(s,paste0(paste0('Cpsi',0:2,'(',i-1,')'),colnames(psi_mat)[j]))
  }
  for (j in 1:ncol(r_mat)) {
    s=c(s,paste0('R0.',colnames(r_mat)[j]))
    if (data$nseasons>1)
      for (i in 2:data$nseasons) s=c(s,paste0(paste0('CR',0:2,'(',i-1,')'),colnames(r_mat)[j]))
  }
  colnames(temp.unitcov)=gsub("(Intercept)","int",s,fixed=TRUE)

  temp.pao<-data
  temp.pao$unitcov<-temp.unitcov;   temp.pao$nunitcov<-ncol(temp.unitcov)

  n2=nrow(p_mat)/2; temp.survcov<-matrix(c(delta_mat,p_mat[1:n2,],p_mat[n2+1:n2,]),nrow=data$nunits*data$nsurveys)
  s=NULL
  for (j in 1:ncol(delta_mat)) s=c(s,paste0(paste0('delta.'),colnames(delta_mat)[j]))
  for (ii in 1:2) {
    for (j in 1:ncol(p_mat)) s=c(s,paste0(paste0('p',ii,'.'),colnames(p_mat)[j]))
  }
  colnames(temp.survcov)=gsub("(Intercept)","int",s,fixed=TRUE)

  temp.cov<-data.frame(temp.unitcov,temp.survcov)

  ### bind each block of p covariates to cov
  temp.unitcov<-as.matrix(p_mat[1:(data$nunits*data$nsurveys),])
  x<-gsub("(Intercept)","int",colnames(p_mat),fixed=TRUE)
  colnames(temp.unitcov)<-paste("p1.",x,sep="")

  temp.survcov<-as.matrix(p_mat[(1+data$nunits*data$nsurveys):(2*data$nunits*data$nsurveys),])
  x<-gsub("(Intercept)","int",colnames(p_mat),fixed=TRUE)
  colnames(temp.survcov)<-paste("p2.",x,sep="")

  temp.pao$survcov<-temp.cov; temp.pao$nsurvcov<-ncol(temp.cov)
  if (is.null(paoname)) temp.pao$paoname="paodata.pao" else temp.pao$paoname=paoname

  #writePao(temp.pao)

  ## create design matrices and run PRESENCE
  temp<-grep("^psi|^Cpsi",colnames(temp.pao$unitcov),value=T)
  psi.dm<-matrix(temp,ncol=ncol(psi_mat),byrow=F)
  temp<-grep("^R0|^CR",colnames(temp.pao$unitcov),value=T)
  r.dm<-matrix(temp, ncol=ncol(r_mat),byrow=F)
  psi.dm=diagbind(psi.dm,r.dm)
#  psi.dm<-cbind(rbind(psi.dm,matrix(0,(1+3*(data$nseasons-1)),ncol(psi_mat))),
#                rbind(matrix(0,(1+3*(data$nseasons-1)),ncol(r_mat)),r.dm))
  v="psi0"; if (data$nseasons>1) v=c(v,paste0("Cpsi",rep(0:2,data$nseasons-1),"(",rep(1:(data$nseasons-1),each=3),")"))
  rownames(psi.dm)=c(v,gsub("psi","R",v)); colnames(psi.dm)<-paste0("a",1:ncol(psi.dm))

  temp<-grep("delta.",colnames(temp.pao$survcov),fixed=TRUE,value=T)
  delta.dm<-matrix(rep(temp,temp.pao$nsurveys),nrow=temp.pao$nsurveys,ncol=ncol(delta_mat), byrow=TRUE)
  temp<-grep("p1.",colnames(temp.pao$survcov),fixed=TRUE,value=T)
  p.dm<-matrix(rep(temp,temp.pao$nsurveys),nrow=temp.pao$nsurveys,ncol=ncol(p_mat),byrow=TRUE)
  temp<-grep("p2.",colnames(temp.pao$survcov),fixed=TRUE,value=T)
  p.dm<-rbind(p.dm,matrix(rep(temp,temp.pao$nsurveys),nrow=temp.pao$nsurveys,ncol=ncol(p_mat),byrow=T))

  v=paste0("p1(",unlist(lapply(1:temp.pao$nseasons, function(xx){rep(xx,temp.pao$nsurveyseason[xx])})),
           "-",  unlist(lapply(1:temp.pao$nseasons, function(xx){1:temp.pao$nsurveyseason[xx]})),")")
  if (data$nseasons>1) {

    p.dm=diagbind(delta.dm,p.dm)
    rownames(p.dm)=c(gsub('p1','delta',v),v,gsub('p1','p2',v))
    colnames(p.dm)<-paste0("e",1:ncol(p.dm))

    psi.dm=simplify_dm(psi.dm,temp.pao);  p.dm=simplify_dm(p.dm,temp.pao)
    temp.pao<-clean_covs_from_pao(c(psi.dm,p.dm),temp.pao)        #  get rid of site/survey covariates not in design matrices
    #s1=as.character(colnames(temp.pao$unitcov))
    #s2=unique(as.character(c(psi.dm,p.dm)))
    #s3=NULL; for (i in 1:length(s1)) if (s1[i] %in% s2) s3=c(s3,s1[i])
    #if (!is.null(s3)) {
    #  temp.pao$unitcov=temp.pao$unitcov[,s3,drop=F]; temp.pao$nunitcov=ncol(temp.pao$unitcov)
    #} else temp.pao$unitcov=temp.pao$nunitcov=0

    #                  #  get rid of survey covariates not in design matrix
    #s1=as.character(colnames(temp.pao$survcov))
    #s2=unique(as.character(p.dm))
    #s3=NULL; for (i in 1:length(s1)) if (s1[i] %in% s2) s3=c(s3,s1[i])
    #if (!is.null(s3)) { temp.pao$survcov=temp.pao$survcov[,s3,drop=F]; temp.pao$nsurvcov=ncol(temp.pao$survcov)
    #} else temp.pao$nsurvcov=temp.pao$survcov=0

    if(!is.null(fixed)) fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,NULL,NULL,NULL,p.dm),rownames)))

    rv<-runPresence(temp.pao,list(psi.dm,NULL,NULL,NULL,p.dm,NULL),
                              model=model,modname=modname,
                              fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  } else {
    rownames(p.dm)=c(v,gsub('p1','p2',v));   colnames(p.dm)=paste0("b",1:ncol(p_mat))
    rownames(delta.dm)=gsub('p1','delta',v); colnames(delta.dm)=paste0("c",1:ncol(delta.dm))

    psi.dm=simplify_dm(psi.dm,temp.pao);  p.dm=simplify_dm(p.dm,temp.pao); delta.dm=simplify_dm(delta.dm,temp.pao)
    temp.pao<-clean_covs_from_pao(c(psi.dm,p.dm,delta.dm),temp.pao)        #  get rid of site/survey covariates not in design matrices

    if(!is.null(fixed)){
      fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm,delta.dm),rownames)))
    }

    rv<-runPresence(temp.pao,list(psi.dm,p.dm,delta.dm,NULL,NULL,NULL),
                              model=model,modname=modname,
                              fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  }

  ################################################################
  #### Extract results from output file
  ################################################################

  npar<-(ncol(psi.dm)+ncol(p.dm)); if (data$nseasons==1) npar=ncol(psi.dm)+ncol(p.dm)+ncol(delta.dm)
  v=readLines(outfile); l=length(v);

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  v=v[-1:-i]

  ### get Untransformed coefficients (Beta's)
  i=grep("^A.+psi.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
  j=grep("^A.+R[012].+ : ",v); r.coeff=getbetas(v[j],r.dm)
  k=grep("^[B-F].+ p.+ : ",v); p.coeff=getbetas(v[k],p.dm)
  l=grep("^[A-F].+delta.+ : ",v); delta.coeff=getbetas(v[l],delta.dm)

  ### get_VC matrix
  ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)
  instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC=as.matrix(VC[1:length(i), 1:length(i)])
  r.VC  =as.matrix(VC[length(i)+1:length(j), length(i)+1:length(j)])
  p.VC  =as.matrix(VC[length(c(i,j))+1:length(k), length(c(i,j))+1:length(k)])
  delta.VC  =as.matrix(VC[length(c(i,j,k))+1:length(l), length(c(i,j,k))+1:length(l)])

  psi.est=p.est=psi_c.est=NULL
  VCoutopt=miscopts[1]
  if (! (VCoutopt %in% c("nose","betavc","noreal"))) {  ## return real parameters
    ## Calc real parameter estimates
    psi.est=calc_real(cov=psi_mat,coeff=psi.coeff$est,VC=psi.VC,conf=conf,rownames=paste0('psi_',temp.pao$unitnames))
    r.est  =calc_real(cov=r_mat,coeff=r.coeff$est,VC=r.VC,conf=conf,rownames=paste0('R_',temp.pao$unitnames))
    plbl=paste0(rep(rownames(p.dm),each=NU),rep(temp.pao$unitnames,nrow(p.dm)))
    p.est=calc_real(cov=p_mat,coeff=p.coeff$est,VC=p.VC,conf=conf,rownames=plbl,fixed=fixed)
    dlbl=paste0(rep(rownames(delta.dm),each=NU),rep(temp.pao$unitnames,nrow(delta.dm)))
    delta.est=calc_real(cov=delta_mat,coeff=delta.coeff$est,VC=delta.VC,conf=conf,rownames=dlbl,fixed=fixed)

    Cpsi0.est<-data.frame(); CR0.est<-data.frame()
    Cpsi1.est<-data.frame(); CR1.est<-data.frame()
    Cpsi2.est<-data.frame(); CR2.est<-data.frame()
    if (data$nseasons>1) {
      for (ii in 1:(data$nseasons-1)){
        Cpsi0.est<-#calc_real(cov=psi_mat,coeff=psi.coeff$est,VC=psi.VC,conf=conf,rownames=paste0('psi_',temp.pao$unitnames))
                    rbind(Cpsi0.est,get_real(real=paste("<Cpsi0(",ii,")>",sep=""),v=v,con.offset=1,
                    n1=1,nunits=temp.pao$nunits,
                    row.names=paste(temp.pao$unitnames,"_",ii,sep=""),
                    fixed=fixed,real.offset=1+(ii-1)*3))
        Cpsi1.est<-rbind(Cpsi1.est,get_real(real=paste("<Cpsi1(",ii,")>",sep=""),v=v,con.offset=1,
                                        n1=1,nunits=temp.pao$nunits,
                                        row.names=paste(temp.pao$unitnames,"_",ii,sep=""),
                                        fixed=fixed,real.offset=2+(ii-1)*3))
        Cpsi2.est<-rbind(Cpsi2.est,get_real(real=paste("<Cpsi2(",ii,")>",sep=""),v=v,con.offset=1,
                                        n1=1,nunits=temp.pao$nunits,
                                        row.names=paste(temp.pao$unitnames,"_",ii,sep=""),
                                        fixed=fixed,real.offset=3+(ii-1)*3))
        CR0.est<-rbind(CR0.est,get_real(real=paste("<CR0(",ii,")>",sep=""),v=v,con.offset=1,
                                      n1=1,nunits=temp.pao$nunits,
                                      row.names=paste(temp.pao$unitnames,"_",ii,sep=""),
                                      fixed=fixed,real.offset=2+3*(data$nseasons-1+ii-1)))
        CR1.est<-rbind(CR1.est,get_real(real=paste("<CR1(",ii,")>",sep=""),v=v,con.offset=1,
                                      n1=1,nunits=temp.pao$nunits,
                                      row.names=paste(temp.pao$unitnames,"_",ii,sep=""),
                                      fixed=fixed,real.offset=3+3*(data$nseasons-1+ii-1)))
        CR2.est<-rbind(CR2.est,get_real(real=paste("<CR2(",ii,")>",sep=""),v=v,con.offset=1,
                                      n1=1,nunits=temp.pao$nunits,
                                      row.names=paste(temp.pao$unitnames,"_",ii,sep=""),
                                      fixed=fixed,real.offset=4+3*(data$nseasons-1+ii-1)))
      }
    }
    getp <- function(x) {
      i=grep(x,v); vv=NULL
      for (k in 1:length(i)) {
        j=2
        while (nchar(v[i[k]+j])>10) { vv=c(vv,v[i[k]+j]); j=j+1 }
#       for (j in 1:length(i)) {
#         if (nchar(v[i[j]+3])<10) vv=c(vv,rep(v[i[j]+2],temp.pao$nunits)) else vv=c(vv,v[i[j]+1+1:temp.pao$nunits])
#       }
      }
      v2=gsub('.+: +','',vv); v3=gsub(',fixed','',gsub('-,','',gsub(' +',',',v2)))
      p1.est=matrix(suppressWarnings(as.numeric(unlist(strsplit(v3,',')))),ncol=4,byrow=T)
      colnames(p1.est)=c('est','se','lower_0.95','upper_0.95');
      v4=paste0(rep(temp.pao$unitnames,temp.pao$nsurveys),"_",rep(temp.pao$surveynames,each=temp.pao$nunits))
      rownames(p1.est)=gsub(' +','',gsub(':.+','',vv))
      return(p1.est)
    }

    p1.est<-getp("estimates of .p1")
    p2.est<-getp("estimates of .p2")
    delta.est<-getp("estimates of .delta")
  }
  ##### check for warnings
  warn.conv<-check_conv_warn(v,quiet)
  warn.VC<-check_VC_warn(v)

  result<-list(modname=modname,
               model=list(psi=psi,r=r,delta=delta,p=p),dmat=list(psi=psi.dm,p=p.dm),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike,npar=npar, aic=aic,
               beta=list(psi=  as.data.frame(psi.coeff),  psi.VC=psi.VC,
                         r  =  as.data.frame(r.coeff),    r.VC=  r.VC,
                         delta=as.data.frame(delta.coeff),delta.VC=delta.VC,
                         p  =  as.data.frame(p.coeff),    p.VC=p.VC,VC=VC),
               real=list(psi=as.data.frame(psi.est),
                         Cpsi0=as.data.frame(Cpsi0.est),
                         Cpsi1=as.data.frame(Cpsi1.est),
                         Cpsi2=as.data.frame(Cpsi2.est),
                         r=as.data.frame(r.est),
                         CR0=as.data.frame(CR0.est),
                         CR1=as.data.frame(CR1.est),
                         CR2=as.data.frame(CR2.est),
                         delta=as.data.frame(delta.est),
                         p1=as.data.frame(p1.est),
                         p2=as.data.frame(p2.est)),
               warnings=list(conv=warn.conv,VC=warn.VC))
  if (data$nseasons<=1) result$dmat$delta=delta.dm

  class(result)<-c("occMod","doMs2")
  return(result)

}
