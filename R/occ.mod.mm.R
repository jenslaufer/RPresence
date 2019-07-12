occMod_mm<-function(psi=call(),psi.cov=data$unitcov,
                     theta=call(),theta.cov=data$unitcov,
                     p=call(),p.cov=data$survcov,
                     gamma=call(), gamma.cov=data$survcov,
                     epsilon=call(), epsilon.cov=data$survcov,
                     modname=NULL, paoname=NULL, outfile,
                     model=103,fixed=NULL,initvals=NULL,data,miscopts=''){
  #' Fit multi-method (or multi-scale) static occupancy/single season model
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="so.mm"}.
  #'
  #' @param psi the right-hand side of the formula for the model to fit for occupancy probability.
  #' @param psi.cov a data frame containing the unit-specific covariates to use for the occupancy component of the model.
  #' @param theta the right-hand side of the formula for the model to fit for secondary-scale occupancy probability.
  #' @param theta.cov a data frame containing the secondary-scale specific covariates to use for the secondary-scaleoccupancy component of the model.
  #' @param p the right-hand side of the formula for the model to fit for detection probability.
  #' @param p.cov a data frame containing the survey-specific covariates to use for the detection component of the model.
  #' @param gamma the right-hand side of the formula for the model to fit for colonization.
  #' @param gamma.cov a data frame containing the site-specific covariates to use for the colonization component of the model.
  #' @param epsilon the right-hand side of the formula for the model to fit for extinction.
  #' @param epsilon.cov a data frame containing the site-specific covariates to use for the extinction component of the model.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param model the PRESENCE model code. DO NOT CHANGE.
  #' @param fixed a single-column matrix containing values for real parameters to be fixed at. \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals initial values for the beta parameters at which PRESENCE begins the optimisation. The default values in PRESENCE is 0.
  #' @param data the \code{pao} data object containing the detection data and other information.
  #' @param miscopts see \code{\link{occMod}}
  #'
  #' @return list of class \code{"occMod"} and \code{"soMm"}.
  #'
  #' \code{occMod$beta} contains the objects:
  #'   \item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence.}
  #'   \item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence.}
  #'   \item{theta}{estimated logistic regression coefficients and standard errors for probability of occurrence at the secondary scale.}
  #'   \item{theta.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence at the secondary scale.}
  #'   \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'   \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'   \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'
  #' \code{occMod$real} contains the objects:
  #'   \item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard errors and limits of 95\% confidence interval.}
  #'   \item{theta}{estimated probabilities of occurrence at the secondary scale for each sampling unit, along with standard errors and limits of 95\% confidence interval.}
  #'   \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits of 95\% confidence interval.}
  #'   \item{psi_c}{estimated probabilities of occurrence given the detection history for each sampling unit, along with standard errors and limits of 95\% confidence interval. Will be \code{=1} for any unit where the species was detected at least once.}
  #'@examples
  #' #
  #' #  sim_so_mm.R - Single-season, multi-method example for RPresence
  #' #                   where theta is modelled as a function of a survey covariate
  #' #
  #' rm(list=ls()); library(RPresence); setwd('~')
  #' #
  #' #  simulate data...
  #' #
  #' N=1000; K=5; M=2             #  N=number of sites, K=number of surveys, M=number of methods
  #' psi=.75;                     #  Prob of occupancy (constant for all sites)
  #' X=matrix(round(rnorm(N*K),4),N,K)   #  Generate site and survey-specific covariate for local occupancy (theta)
  #' b0=0.5; b1=-1;               #  b0=intercept, b1=effect of covariate X on theta
  #' theta=plogis(b0+b1*X)        #  logit(local occupancy) = b0 + b1*X
  #' p=c(.9,.8)                   #  p[1]=detection prob for method1, p[2]=detection prob for method 2
  #'
  #' v=hist(plogis(b0+b1*X),xlab='local occupancy/use') # plot range of true local occupancy...
  #'
  #' occupied_sites=0+(runif(N)<psi)    #  site occupied if random number < Prob of occ
  #'
  #' locally_occ_sites=NULL     #  site locally occ if random number < theta and occupied
  #' for (i in 1:K) locally_occ_sites=cbind(locally_occ_sites,(0+(runif(N)<theta[,i]))*occupied_sites)
  #'
  #' detected=NULL    #  simulate detection...
  #' for (i in 1:K){    #  for each survey...
  #'   for (m in 1:M)   #    for each method...
  #'                    #      detection=1 if rand number<p(method) and locally occupied in survey
  #'     detected=cbind(detected,(0+runif(N)<p[m])*locally_occ_sites[,i])
  #' }
  #' Xcov=NULL   #  For the multi-method model, we need to repeat each survey-specific covariate
  #'             #  for each method.  For example, we have 5 surveys and 2 methods, so we need
  #'             #  10 columns for the survey covariate (X1,X1,X2,X2,X3,X3,X4,X4,X5,X5)
  #' for (i in 1:K){    #  for each survey...
  #'   for (m in 1:M)
  #'     Xcov=cbind(Xcov,X[,i])  #  save covariate (repeated for each method) as survey covariate
  #' }
  #' #
  #' #       generate Presence input object...
  #' #
  #' pao=createPao(data=detected, survcov=data.frame(X=as.numeric(Xcov)),nmethods=2)
  #' #
  #' #       run multi-method model...
  #' #
  #' mmod1=occMod(model=list(psi~1,theta~X,p~DEVICE),cov.list=list(X),data=pao,type='so.mm',outfile='tmp.out1')
  #' print(unique(mmod1$real$psi))  ## print real estimates of psi
  #' print(mmod1$beta$theta)        ## print beta estimates of theta
  #' print(unique(mmod1$real$p))    ## print real estiamtes of p

  options("na.action"="na.pass")
  #### create covariate table for theta
  nsites=data$nunits; nsrvys=data$nsurveys; nmeth=data$nmethods; Msurveys=nsrvys/nmeth; nseasns=length(data$nsurveyseason)
  psicov<-data.frame(rep(1,nsites))
  if(is.null(psi.cov)) psi.cov<-psicov else psi.cov=data.frame(psi.cov,psicov)
  thetacov1=rep(1:Msurveys,each=nsites*nmeth); thetacov2=as.factor(thetacov1)
  thetacov1=data.frame(TREND=thetacov1,PERIOD=thetacov2)
  if (is.null(theta.cov)) theta.cov=data.frame(thetacov1,psi.cov) else theta.cov=data.frame(theta.cov,thetacov1,psi.cov)

  ## add period and device covariates for detection
  pper=NULL; for (i in 1:Msurveys) pper=c(pper,rep(i,nsites*nmeth));
  pdev=NULL; for (i in 1:Msurveys) pdev=c(pdev,rep(1:nmeth,each=nsites))
  pper=as.factor(pper); pdev=as.factor(pdev)
  if (is.null(p.cov)) p.cov=data.frame(PERIOD=pper,DEVICE=pdev) else p.cov=data.frame(p.cov,PERIOD=pper,DEVICE=pdev)

  ##################
  psi_mat<-model.matrix(psi,psi.cov);
  colnames(psi_mat)=gsub('\\(Intercept\\)','int',colnames(psi_mat))
  colnames(psi_mat)<-paste0("psi.",colnames(psi_mat)) # coeffs. extracted later by matching '.psi.'

  theta_mat<-model.matrix(theta,theta.cov);
  colnames(theta_mat)=gsub('\\(Intercept\\)','int',colnames(theta_mat))
  colnames(theta_mat)<-paste0("theta.",colnames(theta_mat)) # coeffs. extracted later by matching '.theta.'

  p_mat<-model.matrix(p,p.cov);
  colnames(p_mat)=gsub('\\(Intercept\\)','int',colnames(p_mat))
  colnames(p_mat)<-paste0("p.",colnames(p_mat)) # coeffs. extracted later by matching '.p.'

  temp.cov=psi_mat; s=colnames(psi_mat)

  if (!is.null(gamma)) {
    if (is.null(gamma.cov))
      gamma_mat<-model.matrix(gamma,data.frame(rep(1,(data$nunits*data$nseasons-1)))) else
      try(gamma_mat <-  model.matrix(gamma,gamma.cov))
    if (is.null(epsilon.cov))
      epsilon_mat<-model.matrix(epsilon,data.frame(rep(1,(data$nunits*data$nseasons-1)))) else
      epsilon_mat=  model.matrix(epsilon,epsilon.cov)
    for(ii in 1:(data$nseasons-1)){
    temp.cov=cbind(temp.cov, gamma_mat[  (ii-1)*data$nunits+1:data$nunits,],
                   epsilon_mat[(ii-1)*data$nunits+1:data$nunits,])
    }
    s=c(colnames(psi_mat),
      paste0(rep(c(paste0("gamma.",  colnames(gamma_mat)),
                   paste0("epsilon.",colnames(epsilon_mat))),data$nseasons-1),
             rep(1:(data$nseasons-1),each=(ncol(gamma_mat)+ncol(epsilon_mat)))))
  }
  colnames(temp.cov)=gsub("(Intercept)","int",s,fixed=T)

  theta_p_mat=cbind(theta_mat,p_mat); rownames(theta_p_mat)=rownames(p_mat)
  colnames(theta_p_mat)=gsub("(Intercept)","int",colnames(theta_p_mat),fixed=T)

  ## create and output temporary pao file
  ### bind each block of theta covariates to unitcov

  temp.pao=data   ## create and output temporary pao file
  temp.pao$unitcov=as.data.frame(temp.cov); temp.pao$nunitcov=ncol(temp.cov)
  temp.pao$nunitcov<-ncol(temp.cov)
  temp.pao$survcov<-as.data.frame(theta_p_mat)
  temp.pao$nsurvcov<-ncol(theta_p_mat)
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)
  if (sum(temp.pao$nsurveyseason)<temp.pao$nsurveys) temp.pao$nsurveyseason=rep(temp.pao$nsurveyseason,nmeth)
  nseasns=length(temp.pao$nsurveyseason)
  #writePao(temp.pao)

  ## create design matrices and run PRESENCE
  psi.dm<-matrix(colnames(psi_mat),nrow=1); rownames(psi.dm)="psi";
  theta.dm<-matrix(rep(colnames(theta_mat),Msurveys),nrow=Msurveys,byrow=TRUE)
  theta.dm=gsub("(Intercept)","int",theta.dm,fixed=T)
  rownames(theta.dm)=paste0("theta",1:Msurveys)
  psi.dm<-diagbind(psi.dm,theta.dm); psi.dm=gsub("(Intercept)","int",psi.dm,fixed=T)
  colnames(psi.dm)<-paste0("a",1:ncol(psi.dm))

  p.dm<-matrix(rep(colnames(p_mat),nsrvys),nrow=nsrvys,byrow=TRUE)
  rownames(p.dm)<-paste0("p",1:nmeth,'(',rep(1:(nseasns/nmeth),each=Msurveys),'-',rep(1:Msurveys,each=nmeth),')')
  colnames(p.dm)<-paste0("b",1:ncol(p.dm))

  gamma.dm=epsilon.dm=NULL
  if (!is.null(gamma)) {
    index=grep("gamma.",colnames(temp.pao$unitcov),fixed=TRUE)
    gamma.dm=matrix(colnames(temp.pao$unitcov)[index], (temp.pao$nseasons-1),ncol(gamma_mat),byrow=TRUE)
    rownames(gamma.dm)=paste0("gamma",1:(temp.pao$nseasons-1)); colnames(gamma.dm)=paste0("b",1:ncol(gamma_mat))
    gamma.dm=simplify_dm(gamma.dm,temp.pao)

    index=grep("epsilon.",colnames(temp.pao$unitcov),fixed=TRUE)
    epsilon.dm=matrix(colnames(temp.pao$unitcov)[index], (temp.pao$nseasons-1),ncol(epsilon_mat),byrow=TRUE)
    rownames(epsilon.dm)=paste0("epsilon",1:(temp.pao$nseasons-1)); colnames(epsilon.dm)=paste0("c",1:ncol(epsilon_mat))
    epsilon.dm=simplify_dm(epsilon.dm,temp.pao)
  }
  psi.dm=simplify_dm(psi.dm,temp.pao);  p.dm=simplify_dm(p.dm,temp.pao); theta.dm=simplify_dm(theta.dm,temp.pao);
  if(!is.null(fixed)){
    fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm),rownames)))
  }
  rv<-runPresence(temp.pao,list(psi.dm,p.dm,gamma.dm,epsilon.dm,NULL,NULL),
                  modname=modname,model=2000,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)

  ################################################################
  #### Extract results from output file
  ################################################################
  npar<-ncol(psi.dm)+ncol(p.dm)
  v=readLines(outfile); l=length(v);
  ## find AIC value for each model
  instns=grep("-2log(likelihood)",v,fixed=TRUE)
  z<-get_line(v,instns[1])
  neg2loglike<-as.numeric(z[4])
  z<-get_line(v,instns[1]+1)
  aic<-as.numeric(z[3])

  ### get_coefficients
  instns<-grep("Untransformed",v,fixed=TRUE)
  extract<-v[(instns[1]+3):(instns[1]+2+npar)]
  index<-grep(" psi",extract,fixed=TRUE)
  psi.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(psi_mat),
                       cov.names=colnames(psi_mat),b.string="a",b.count=0)

  index<-grep("A.+ .theta.+ :",extract)
  theta.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(theta_mat),
                         cov.names=colnames(theta_mat),b.string="a",b.count=ncol(psi_mat))

  index<-grep("B.+ .p.",extract,fixed=FALSE)
  p.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(p_mat),
                     cov.names=colnames(p_mat),b.string="b",b.count=0)

  ### get_VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  names<-c(rownames(psi.coeff),rownames(theta.coeff),rownames(p.coeff))
  VC<-get_VC(v=v,offset=instns+1,n.coeff=ncol(psi.dm)+ncol(p.dm),cov.names=names)
  psi.VC<-VC[psi.coeff[,"index"],psi.coeff[,"index"]]
  theta.VC<-VC[theta.coeff[,"index"],theta.coeff[,"index"]]
  p.VC<-VC[p.coeff[,"index"],p.coeff[,"index"]]

  ## Get real parameter estimates from output
  psi.est<-get_real(real="<psi>",v=v,con.offset=1,n1=1,nunits=nsites,
                    row.names=temp.pao$unitnames,fixed=fixed,real.offset=0)
  theta.est<-get_real(real="<theta1>",v=v,con.offset=1,n1=Msurveys,nunits=nsites,
                      row.names=paste0(temp.pao$unitnames,"_",rep(1:Msurveys,each=nsites)),
                      fixed=fixed,real.offset=1)
  i=grep('^p[12]\\(.+ :',v); vi=v[i]; v1=gsub('-,','',gsub(' +',',',gsub('.+: +','',vi)))
  p.est=matrix(as.numeric(unlist(strsplit(v1,','))),ncol=4,byrow=T); colnames(p.est)=colnames(psi.est)
  rownames(p.est)=gsub(' +','_',gsub(' +:.+','',v[i]))

  psi_c.est<-get_real(real="Psi-conditional",v=v,con.offset=2,n1=1,nunits=nsites,
                      row.names=temp.pao$unitnames,fixed=fixed,real.offset=0)

  ##### check for warnings
  warn.conv<-check_conv_warn(v)
  warn.VC<-check_VC_warn(v)

  result<-list(modname=modname,
               model=list(psi=psi,theta=theta,p=p),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike,
               npar=(ncol(psi.dm)+ncol(p.dm)), aic=aic,
               beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                         theta=as.data.frame(theta.coeff),theta.VC=theta.VC,
                         p=as.data.frame(p.coeff),p.VC=p.VC,VC=VC),
               real=list(psi=as.data.frame(psi.est),
                         theta=as.data.frame(theta.est),
                         p=as.data.frame(p.est),
                         psi_c=as.data.frame(psi_c.est)),
               warnings=list(conv=warn.conv,VC=warn.VC))
  class(result)<-c("occMod","soMm")
  return(result)

}
