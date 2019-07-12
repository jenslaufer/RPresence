occMod_SO_cd<-function(psi  =call(),   psi.cov=data$unitcov,
                        theta=call(), theta.cov=data$unitcov,
                        p    =call(),     p.cov=data$survcov,
                        th0pi=call(), th0pi.cov=data$unitcov,
                     modname=NULL, paoname=NULL,outfile, model=103,fixed=NULL,initvals=NULL,data,conf,miscopts=''){
  #' Fit static occupancy, correlated detections model
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="so.cd"}.
  #'
  #' @param psi the right-hand side of the formula for the model to fit for occupancy probability (eg., ~1).
  #' @param psi.cov a data frame containing the unit-specific covariates to use for the occupancy component of the model.
  #' @param theta the right-hand side of the formula for the model to fit for secondary-scale occupancy probability.
  #' (eg., ~PRIME -> th0 != th1, ~SURVEY -> th0(1)=th1(1), th0(2)=th1(2),...)
  #' @param theta.cov a data frame containing the secondary-scale specific covariates to use for the secondary-scaleoccupancy component of the model.
  #' @param p the right-hand side of the formula for the model to fit for detection probability.
  #' @param p.cov a data frame containing the survey-specific covariates to use for the detection component of the model.
  #' @param th0pi the right-hand side of the formula for the model to fit for local use before 1st segment.
  #' @param th0pi.cov a data frame containing the unit-specific covariates to use for th0pi.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #'@param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param model the PRESENCE model code. DO NOT CHANGE.
  #' @param fixed a single-column matrix containing values for real parameters to be fixed at. \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals initial values for the beta parameters at which PRESENCE begins the optimisation. The default values in PRESENCE is 0.
  #' @param data the \code{pao} data object containing the detection data and other information.
  #'@param conf level for confidence interval (may be vector valued).
  #' @param miscopts see \code{\link{occMod}}
  #'
  #' @return list of class \code{"occMod"} and \code{"soCd"}.
  #'
  #' \code{occMod$beta} contains the objects:
  #'   \item{psi}{estimated logistic regression coefficients and standard errors for probability of occurrence.}
  #'   \item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence.}
  #'   \item{theta}{estimated logistic regression coefficients and standard errors for probability of occurrence at the secondary scale.}
  #'   \item{theta.VC}{variance-covariance matrix for logistic regression coefficients for probability of occurrence at the secondary scale.}
  #'   \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'   \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'   \item{th0pi}{estimated logistic regression coefficients and standard errors for probability of local-use before 1st segment.}
  #'   \item{th0pi.VC}{variance-covariance matrix for logistic regression coefficients for th0pi.}
  #'   \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'
  #' \code{occMod$real} contains the objects:
  #'   \item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard errors and limits of 95\% confidence interval.}
  #'   \item{theta}{estimated probabilities of occurrence at the secondary scale for each sampling unit, along with standard errors and limits of 95\% confidence interval.}
  #'   \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits of 95\% confidence interval.}
  #'   \item{psi_c}{estimated probabilities of occurrence given the detection history for each sampling unit, along with standard errors and limits of 95\% confidence interval. Will be \code{=1} for any unit where the species was detected at least once.}
  #'
  #' @examples
  #' #   example 1 - read previously created input file and run a model
  #'data=readPao('genpres.pao')  #  read simulated data from program GENPRES
  #'mod1=occMod(model=list(psi~1, theta~PRIME, p~1, th0pi~1), data=data, type="so.cd")
  #'
  #' #   example 2 - simulate data and run a mode
  #'sim_corr_det_data <- function(sites=500,surveys=9,psi=.9,th=c(.2,.8),p=.6,th0pi=0) {
  #'   # simulates single-season correlated detections data
  #'   h=matrix(0,sites,surveys); px=c(0,p) # px=detection prob for each state (local-unocc, local-occ)
  #'   occ=0+(runif(sites)<psi);            # randomly assign occupancy state of each site (0=unocc, 1=occ)
  #'   locc=occ*(runif(sites)<th0pi);       # randomly assign local-occ state before 1st segment of trail
  #'   for (i in 1:surveys) {
  #'      locc=occ*(runif(sites)<th[locc+1]);   # local-occ for survey i depends on previous local-occ
  #'      h[,i]=locc*(runif(sites)<px[locc+1])  # detection depends on local-occ state (p=0 if not local-occ)
  #'   }
  #'return(h)
  #'}
  #'x=sim_corr_det_data(sites=1000,surveys=15)
  #'data=createPao(x)
  #'mod1=occMod(model=list(psi~1,       # constant psi
  #'                      theta~PRIME, # theta not equal theta'
  #'                      p~1,         # constant detection
  #'                      th0pi~1      # constant prop unocc before 1st segment
  #'), data=data, type="so.cd")  #  run ss-cd model
  #'print_one_site_estimates(mod1,site=1)

  #'

  options("na.action"="na.pass")
  #### create covariate table for theta
  nsites=data$nunits; nsrvys=data$nsurveys; pst=function(x) paste(as.character(x),collapse='')
  modname=paste0('psi(',pst(psi),')theta(',pst(theta),')p(',pst(p),')th0pi(',pst(th0pi),')')
  het_flag=length(grep("HET",as.character(p))>0)

  thSRVY=rep(rep(1:nsrvys,each=nsites),2); thPRIME=rep(1:2,each=nsites*nsrvys)
  thetacov1=data.frame(SURVEY=as.factor(thSRVY),PRIME=as.factor(thPRIME))
  if (is.null(theta.cov)) theta.cov=data.frame(thetacov1) else theta.cov=data.frame(theta.cov,thetacov1)

  ## add period and device covariates for detection
  pper=rep(1:nsrvys,each=nsites); if (het_flag) { pper=c(pper,pper); p.cov=rbind(p.cov,p.cov)}
  l=length(pper); l2=length(pper)/2; pper=as.factor(pper); het=as.factor(0+((1:l)>l2))
  if (is.null(p.cov)) p.cov=data.frame(SURVEY=pper,HET=het) else p.cov=data.frame(p.cov,SURVEY=pper,HET=het)

  ##################
  psi_mat=matrix(1,nsites,1); rownames(psi_mat)=1:nsites; colnames(psi_mat)="psi.int"
  if (!is.null(psi.cov)) psi_mat=model.matrix(psi,psi.cov);
  theta_mat=model.matrix(theta,theta.cov); colnames(theta_mat)=gsub('.Intercept.','int',colnames(theta_mat))
  p_mat=model.matrix(p,p.cov); colnames(p_mat)=gsub('.Intercept.','p.int',colnames(p_mat))
  th0pi_mat=matrix(1,nsites,1); rownames(th0pi_mat)=1:nsites; colnames(th0pi_mat)="th0pi.int"
  if (!is.null(th0pi.cov)) th0pi_mat=model.matrix(th0pi,th0pi.cov);
  ## create and output temporary pao file
  ### bind each block of theta covariates to unitcov

  v1=paste0("th",rep(0:1,each=ncol(theta_mat)),".",rep(colnames(theta_mat),2)); v1
  v2=paste0(rep(colnames(p_mat),1+het_flag),"_",rep(c("p1","p2")[1:1+het_flag],each=ncol(p_mat)))

  temp.pao=data
  temp.pao$unitcov=as.data.frame(cbind(psi_mat,th0pi_mat))
  temp.pao$nunitcov=ncol(temp.pao$unitcov); l=nsites*nsrvys; i=seq(1,nrow(theta_mat),nsites)

  temp.th.cov=cbind(theta_mat[1:l,],theta_mat[l+1:l,])
  colnames(temp.th.cov)=v1

  if (het_flag) {temp.p.cov=cbind(p_mat[1:l,],p_mat[l+1:l,]); colnames(temp.p.cov)[1:(2*ncol(p_mat))]=v2} else {
    temp.p.cov=p_mat; colnames(temp.p.cov)=v2[1:ncol(p_mat)]
  }
  temp.pao$survcov=as.data.frame(cbind(temp.th.cov,temp.p.cov)); temp.pao$nsurvcov=ncol(temp.pao$survcov)

  temp.pao$paoname=ifelse(is.null(paoname),"paodata.pao",paoname)
  if (length(temp.pao$unitnames)<1) temp.pao$unitnames=paste0('unit',1:temp.pao$nunits)
  if (length(temp.pao$surveynames)<1) temp.pao$surveynames=1:temp.pao$nsurveys

  ## create design matrices and run PRESENCE
  psi.dm=matrix(colnames(psi_mat),nrow=1); rownames(psi.dm)="psi"
  theta.dm=matrix(rep(matrix(v1,2,byrow=T),each=nsrvys),nrow=2*nsrvys)
  rownames(theta.dm)=c(paste0("th0(",1:nsrvys,")"),paste0("th1(",1:nsrvys,")"))
  psi.dm=diagbind(psi.dm,theta.dm); colnames(psi.dm)=paste0("a",1:ncol(psi.dm))
  th0pi.dm=matrix(colnames(th0pi_mat),nrow=1); rownames(th0pi.dm)="th0pi"

  dm=matrix(rep(v2,each=nsrvys),nrow=nsrvys); l=ncol(dm)/2; p.dm=dm
  rownames(p.dm)=paste0("p",1:nrow(p.dm)); colnames(p.dm)=paste0("b",1:ncol(p.dm)); pi.dm=NULL
  if (het_flag) {
    p.dm=rbind(dm[,1:l],dm[,l+1:l]); rownames(p.dm)=c(paste0('p1_',1:nsrvys),paste0('p2_',1:nsrvys));
    pi.dm=matrix(1,1,1,dimnames=list("pi","e2"))
  }
  psi.dm=simplify_dm(psi.dm,temp.pao)
  p.dm=simplify_dm(p.dm,temp.pao)
  th0pi.dm=simplify_dm(th0pi.dm,temp.pao); colnames(th0pi.dm)=paste0('e',1:ncol(th0pi.dm))

  unitcov=NULL; temp.pao$nunitcov=0; lst=get_not01_dm(c(psi.dm,th0pi.dm),names(temp.pao$unitcov))
  if (sum(lst)>0) {
    for (i in 1:length(lst)) if (lst[i]>0) unitcov=cbind(unitcov,temp.pao$unitcov[,i])
    colnames(unitcov)=colnames(temp.pao$unitcov)[which(lst>0)]
  }
  temp.pao$unitcov=unitcov; if (!is.null(unitcov)) temp.pao$nunitcov=ncol(unitcov);

  survcov=NULL; temp.pao$nsurvcov=0; lst=get_not01_dm(c(psi.dm,p.dm),names(temp.pao$survcov))
  if (sum(lst)>0) {
    for (i in 1:length(lst)) if (lst[i]>0) survcov=cbind(survcov,temp.pao$survcov[,i])
    colnames(survcov)=colnames(temp.pao$survcov)[which(lst>0)]
  }
  temp.pao$survcov=survcov; if (!is.null(survcov)) temp.pao$nsurvcov=ncol(survcov)
  #writePao(temp.pao)

  if(!is.null(fixed)){
    fixed$idx=match(fixed$param,unlist(lapply(list(psi.dm,NULL,NULL,p.dm,th0pi.dm,pi.dm),rownames)))
  }
  rv=runPresence(temp.pao,list(psi.dm,NULL,NULL,p.dm,th0pi.dm,pi.dm),
                  model=103,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)

  ### extract results from output file
  v=readLines(outfile)

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  v=v[-1:-i]

  ### get Untransformed coefficients (Beta's)
  i=grep("^A.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
  theta.coeff=psi.coeff[grep("^th",rownames(psi.coeff)),]; psi.coeff=psi.coeff[grep("^p",rownames(psi.coeff)),]
  k=grep("^D.+ : ",v); p.coeff=getbetas(v[k],p.dm)
  l=grep("^E.+ : ",v); th0pi.coeff=getbetas(v[l],th0pi.dm)

  ### get_VC matrix
  ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)
  instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC=as.matrix(VC[1:nrow(psi.coeff), 1:nrow(psi.coeff)])
  theta.VC=as.matrix(VC[nrow(psi.coeff)+1:nrow(theta.coeff), nrow(psi.coeff)+1:nrow(theta.coeff)])
  p.VC  =as.matrix(VC[length(i)+1:length(k), length(i)+1:length(k)])
  th0pi.VC=VC[length(c(i,k))+1:length(l),length(c(i,k))+1:length(l)]

  psi.est=p.est=psi_c.est=NULL
  VCoutopt=miscopts[1]
  if (! (VCoutopt %in% c("nose","betavc","noreal"))) {  ## return real parameters
    ## Get real parameter estimates from output
    psi.est=calc_real(cov=psi_mat,coeff=psi.coeff$est,VC=psi.VC,conf=conf,rownames=paste0('psi_',temp.pao$unitnames))
    s=paste0('th0(',rep(1:temp.pao$nsurveys,each=temp.pao$nunits),')_',rep(temp.pao$unitnames,temp.pao$nsurveys))
    w=gsub('th0','th1',s)
    theta.est=calc_real(cov=theta_mat,coeff=theta.coeff$est,VC=theta.VC,conf=conf,rownames=c(s,w))
    p.est=calc_real(cov=p_mat,coeff=p.coeff$est,VC=p.VC,conf=conf,rownames=gsub('th0','p',s),fixed=fixed)
    th0pi.est=calc_real(cov=th0pi_mat,coeff=th0pi.coeff$est,VC=th0pi.VC,conf=conf,rownames=paste0('th0pi_',temp.pao$unitnames))
  }
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

  ##### check for warnings
  warn.conv=check_conv_warn(v); warn.VC=check_VC_warn(v)

  result=list(modname=modname,
               model=list(psi=psi,theta=theta,p=p),dmat=list(psi=psi.dm,p=p.dm,th0pi=th0pi.dm,pi=pi.dm),
               data=temp.pao,outfile=outfile,
               neg2loglike=neg2loglike,
               npar=(ncol(psi.dm)+ncol(p.dm)), aic=aic,
               beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                         theta=as.data.frame(theta.coeff),theta.VC=theta.VC,
                         p=as.data.frame(p.coeff),p.VC=p.VC,
                         th0pi=as.data.frame(th0pi.coeff),th0pi.VC=th0pi.VC,VC=VC),
               real=list(psi=as.data.frame(psi.est),
#                         theta=data.frame(th0=theta0.est,th1=theta1.est),
                         theta=as.data.frame(theta.est),
                         p=as.data.frame(p.est),th0pi=as.data.frame(th0pi.est)),
               gof=gof,  warnings=list(conv=warn.conv,VC=warn.VC))
  class(result)=c("occMod","soCd")
  return(result)
}
