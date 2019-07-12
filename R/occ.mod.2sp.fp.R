occMod_2SPfp<-function(psi=call(),psi.cov=data$unitcov,p=call(),p.cov=cbind(data$unitcov,data$survcov),
                         omeg=call(),omeg.cov, conf=call, conf.cov,
                         param="PsiBA", sp.contr=TRUE,modname=NULL,paoname=NULL,outfile,model=3000,
                         fixed=NULL,initvals=NULL,data,miscopts=''){

  #' Fit two species, static occupancy (single season) model with false positive detections
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used
  #' with \code{type="so.2sp.cd"}. NOTE THAT THERE MAY BE SOME CHANGES TO HOW THIS MODEL IS IMPLEMENTED
  #' IN THE NEAR FUTURE.
  #'
  #' @param psi The right-hand side of the formula for the model to fit for occupancy probability.
  #' The terms \code{SP} and \code{INT} can be used to define a species effect on occupancy and an
  #' occurrence-level interaction between species accordingly, without them being defined in \code{psi.cov}.
  #' @param psi.cov A data frame containing the unit-specific covariates to use for the occupancy component
  #' of the model.
  #' @param p The right-hand side of the formula for the model to fit for detection probability. The terms
  #' \code{SP}, \code{INT_o} and \code{INT_d} can be used to define a species effect on detection, a
  #' detection-level interaction where the occurrence of one species changes the detection probability of
  #' the other species and a second detection-level interaction where the detection of one species changes
  #' the detection probability of the other species in the same survey. These terms do not have to be
  #' defined as variables in \code{p.cov}.
  #' @param p.cov A data frame containing the survey-specific covariates to use for the detection component
  #' of the model.
  #' @param param The parameterisation to be used, either "\code{psiBA}" or "\code{nu}", which relate to
  #' \code{type="so.2sp.1"} and \code{type="so.2sp.2"} in the function \code{\link{occMod}} respectively.
  #' @param sp.contr \code{TRUE} (default) or \code{FALSE}. Specifies the type of contrast used for the
  #' \code{SP} term. It it used as the \code{contrasts} input in the function \code{\link{contr.treatment}}.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #' @param model The PRESENCE model code. DO NOT CHANGE.
  #' @param maxfn maximum number of function evaluations before aborting optimization.
  #' @param fixed A single-column matrix containing values for real parameters to be fixed at.
  #' \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals Initial values for the beta parameters at which PRESENCE begins the optimisation.
  #' The default values in PRESENCE is 0.
  #' @param conf limits for confidence intervals as a proportion (defalut=0.95 for 95\% conf. interval limits)
  #' @param data The \code{pao} data object containing the detection data and other information.
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param miscopts (see \code{\link{occMod}})
  #'
  #' @return returns a list of class \code{occMod} and \code{so2spCd}
  #'
  #' \code{occMod$beta} contains the objects:
  #'  \item{psi}{estimated logistic regression coefficients and std.err for prob. of occurrence.}
  #'  \item{psi.VC}{var-covar matrix for logistic regression coefficients for prob. of occurrence.}
  #'  \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'  \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'  \item{omeg}{estimated logistic regression coefficients and std.err for prob. of local occurrence before 1st survey.}
  #'  \item{omeg.VC}{var-covar matrix for omeg.}
  #'  \item{conf}{estimated logistic regression coefficients and std.err for prob. of sample is confirmed.}
  #'  \item{conf.VC}{var-covar matrix for conf.}
  #'  \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'  \code{occMod$real} contains the objects:
  #'  \item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard errors and limits of 95\% confidence interval.}
  #'  \item{theta}{estimated prob. of local occurrence in each survey w/ std.err and 95\% conf interval limits.}
  #'  \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits of 95\% confidence interval.}
  #'  \item{th0pi}{estimated prob. of local occurrence before 1st survey w/ std err and 95\% conf. interval limits.}
  #'
  #' @author Jim Hines and Darryl MacKenzie
  #'
  #' @seealso \code{\link{occMod}}
  #' @examples
  #' \dontrun{
  #'# load a csv data file
  #'filename<-system.file("extdata/twosp_exmpl.csv",package="RPresence")
  #'dethist<-read.csv(filename,as.is=T)
  #'
  #'nsites=nrow(dethist); nsrvys=ncol(dethist)         #  set number of sites,surveys from det. history data
  #'dethist=matrix(as.integer(unlist(dethist)),nrow=nsites) # replace missing values (-) with NA
  #'
  #'##          create input "pao" object, for use with occMod function
  #'data=createPao(dethist,unitcov=NULL,survcov=NULL,title="twosp corr.det. example")
  #'
  #'## fit some models
  #'
  #'mod1<-occMod(model=list(psi~SP,   # occupancy species-specific, no interaction, parameters: psiA, psiBA=psiBa
  #'                         theta~SP, # local occ. species-specific: parameters, thetaA0, thetaA1 (no corr. det. model)
  #'                         p~SP,     # detection sp. specific: parms: pA, pB
  #'                         th0pi~1),data=data,type="so.2sp.cd",param="PsiBA")
  #'
  #'mod2=occMod(model=list(psi~SP+INT,    # species and interaction effect (psiA != psiBA != psiBa)
  #'                        theta~SP*PRIME+BAa+INTth,   # thA != thA' != thBA != thBA' != thBa != thBa'
  #'                        p~SP+INT_o+INT_d+INT_so,    # pA != pB != rA != rBA !=rBa
  #'                        th0pi~1),                   # constant prop unocc before 1st segment
  #'             data=data, type="so.2sp.cd", fixed=NULL)
  #'#
  #'tbl=createAicTable(list(mod1,mod2)); print(tbl$table)
  #'
  #' sim_2sp_corr_det_data <- function(sites=100,surveys=12,psiA=.8,psiBA=.3,psiBa=.7,
  #'                                   thA=c(.4,.9),thBA=c(.3,.9),thBa=c(.3,.9),  #  1st is theta, 2nd is theta\'
  #'                                   pA=.6,pB=.66,rA=.7,rBA=.4,rBa=.5,
  #'                                   th0piA=0,th0piBA=0,th0piBa=0) {
  #'               #  simulates single-season 2-species correlated detections data
  #'  psiB=c(psiBa,psiBA); th0piB=c(th0piBa,th0piBA); h=matrix(0,sites,surveys)
  #'  occA=0+(runif(sites)<psiA);               # randomly assign occupancy state for sp. A of each site (0=unocc, 1=occ)
  #'  occB=0+(runif(sites)<psiB[occA+1])        # randomly assign occupancy state for sp. B; depends on occupancy of sp. A
  #'  loccA=occA*(runif(sites)<th0piA);         # randomly assign local-occ state before 1st segment of trail for sp. A
  #'  loccB=occB*(runif(sites)<th0piB[loccA+1]) # randomly assign local-occ state before 1st segment of trail for sp. B
  #'  thB=c(thBa,thBA); pAx=c(pA,rA); pBx=c(pB,rBa,NA,rBA)
  #'  for (i in 1:surveys) {
  #'    loccA=occA*(runif(sites)<thA[loccA+1]);       # local-occ for survey i depends on previous local-occ
  #'    loccB=occB*(runif(sites)<thB[loccB+1+2*occA])
  #'    detA=loccA*(runif(sites)<pAx[occB+1])         # detection depends on local-occ state (p=0 if not local-occ)
  #'    detB=loccB*(runif(sites)<pBx[occA+1+2*detA])
  #'    h[,i]=detA+2*detB  #  detection-history in compressed format: 0=neither sp. detected, 1=sp. A only, 2=sp. B only, 3=both
  #'  }
  #'  return(h)
  #'}
  #'
  #'x=sim_2sp_corr_det_data(1000)          # simulate some data
  #'data=createPao(x)                     # create_pao file from data
  #'#    run a model : psiA()psiBA(),psiBa,thA(),thA'(),thBA(),thBA'(),thBa(),thBa'(),pA(),pB(),rA(),rBA(),rBa()
  #'mod3=occMod(model=list(psi~SP+INT,    # species and interaction effect (psiA != psiBA != psiBa)
  #'           theta~SP*PRIME+BAa+INTth,   # thA != thA' != thBA != thBA' != thBa != thBa'
  #'           p~SP+INT_o+INT_d+INT_so,    # pA != pB != rA != rBA !=rBa
  #'           th0pi~1),                   # constant prop unocc before 1st segment
  #'            data=data, type="so.2sp.cd", fixed=NULL)
  #'print_one_site_estimates(mod3,site=1)
  #' }
  options("na.action"="na.pass")

  detdata=data$det.data
  if (max(data$det.data)<4) detdata=data$det.data+matrix(data$survcov$conf,data$nunits,data$nsurveys)*4

  samp.units<-data$nunits; data$unitnames=paste0('unit',1:samp.units)

  if (is.null(psi.cov)) psi.cov=data.frame(SP=as.factor(c(rep("A",samp.units),rep("B",2*samp.units)))) else {
    if (nrow(psi.cov)==samp.units) psi.cov=rbind(psi.cov,psi.cov,psi.cov)
    psi.cov$SP<-as.factor(c(rep("A",samp.units),rep("B",2*samp.units)))
  }
  contrasts(psi.cov$SP)<-contr.treatment(2,contrasts=sp.contr)
  psi.cov$INT<-as.factor(c(rep(1,2*samp.units),rep(2,samp.units)))  #  INT - psiBA,psiBa interaction
  psi_mat<-model.matrix(psi,psi.cov)

  len<-samp.units*data$nsurveys ## length of covariates in detection design matrix
  if(is.null(p.cov)) p.cov<-data.frame(SURVEY=rep(1:data$nsurveys,each=samp.units)) else
    p.cov$SURVEY=rep(1:data$nsurveys,each=samp.units)
  p.cov=data.frame(apply(p.cov,2,rep,5))
  p.cov$SP<-as.factor(c(rep("A",len), rep("B",len), rep("A",len), rep("B",2*len)))
  seasn=time=NULL; for (i in 1:data$nseasons) seasn=c(seasn,rep(i,samp.units*data$nsurveyseason[i]))
  for (i in 1:data$nsurveys) time=c(time,rep(i,samp.units))
  p.cov$SEASON<-as.factor(rep(seasn,5));
  contrasts(p.cov$SP)<-contr.treatment(2,contrasts=sp.contr)
  p.cov$INT_o<-as.factor(c(rep('1',2*len),rep('2',3*len)))  #  interaction: r != p
  p.cov$INT_d<-as.factor(c(rep(1,4*len),rep(2,len)))        #  interaction: rBa != rBA
  p.cov$INT_so<-as.factor(c(rep(1,3*len),rep(2,2*len)))        #  interaction: rBa != rBA
  p_mat<-model.matrix(p,p.cov);

  if(is.null(omeg.cov)) omeg.cov<-data.frame(SURVEY=rep(1:data$nsurveys,each=samp.units)) else
    omeg.cov$SURVEY=rep(1:data$nsurveys,each=samp.units)
  omeg.cov=data.frame(apply(omeg.cov,2,rep,4))
  omeg.cov$SP<-as.factor(c(rep("A",2*len), rep("B",2*len)))
  seasn=time=NULL; for (i in 1:data$nseasons) seasn=c(seasn,rep(i,samp.units*data$nsurveyseason[i]))
  for (i in 1:data$nsurveys) time=c(time,rep(i,samp.units))
  omeg.cov$SEASON<-as.factor(rep(seasn,4));
  omeg.cov$INT_d<-as.factor(c(rep(1,2*len),rep(2,2*len)))        #  interaction: oa != ob
  omeg.cov$INT_so<-as.factor(rep(c(rep(1,len),rep(2,len)),2))   #  interaction: oa != oA, ob != oB
  omeg_mat<-model.matrix(omeg,omeg.cov)

  if(is.null(conf.cov)) conf.cov<-data.frame(SURVEY=rep(1:data$nsurveys,each=samp.units)) else
    conf.cov$SURVEY=rep(1:data$nsurveys,each=samp.units)
  conf.cov=data.frame(apply(conf.cov,2,rep,2))
  conf.cov$SP<-as.factor(c(rep("A",len), rep("B",len)))
  seasn=time=NULL; for (i in 1:data$nseasons) seasn=c(seasn,rep(i,samp.units*data$nsurveyseason[i]))
  for (i in 1:data$nsurveys) time=c(time,rep(i,samp.units))
  conf.cov$SEASON<-as.factor(rep(seasn,2));
  conf.cov$INT_d<-as.factor(c(rep(1,len),rep(2,len)))        #  interaction: cA!= cB
  conf_mat<-model.matrix(conf,conf.cov)

  psi.label<-c("psiA","psiBA","psiBa"); p.label<-c("pA","pB","rA","rBA","rBa")
  omeg.label<-c('oa','oA','ob','oB'); conf.label<-c('cA','cB')

  ## create and output temporary pao file
  ### bind each block of occupancy covariates to unitcov
  temp.cov<-cbind(psi_mat[1:samp.units,],psi_mat[samp.units+1:samp.units,],psi_mat[2*samp.units+1:samp.units,])
  v1=c(paste0(rep(colnames(psi_mat),3),"_",rep(psi.label,each=ncol(psi_mat))))
  colnames(temp.cov)<-v1
  temp.pao<-data; temp.pao$nunitcov<-ncol(temp.cov); temp.pao$unitcov<-as.data.frame(temp.cov)
  temp.pao$det.data=detdata;
  rownames(temp.pao$unitcov)<-data$unitnames;
  if (length(temp.pao$surveynames)<1) temp.pao$surveynames=1:temp.pao$nsurveys

  ## bind each block of detection covariates
  temp.cov<-cbind(p_mat[1:len,],p_mat[len+1:len,],p_mat[2*len+1:len,],p_mat[3*len+1:len,],p_mat[4*len+1:len,],
                  omeg_mat[1:len,],omeg_mat[len+1:len,],omeg_mat[2*len+1:len,],omeg_mat[3*len+1:len,],
                  conf_mat[1:len,],conf_mat[len+1:len,])
  v2=c(paste0(rep(colnames(p_mat),5),"_",rep(p.label,each=ncol(p_mat))))
  v3=c(paste0(rep(colnames(omeg_mat),4),"_",rep(omeg.label,each=ncol(omeg_mat))))
  v4=c(paste0(rep(colnames(conf_mat),2),"_",rep(conf.label,each=ncol(conf_mat))))
  colnames(temp.cov)<-c(v2,v3,v4)

  temp.pao$nsurvcov<-ncol(temp.cov)
  temp.pao$survcov<-as.data.frame(temp.cov)
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)
  rownames(temp.pao$survcov)<-1:dim(temp.pao$survcov)[1]

  colnames(temp.pao$unitcov)<-gsub(".Intercept.","Int",colnames(temp.pao$unitcov))
  colnames(temp.pao$survcov)<-gsub(".Intercept.","Int",colnames(temp.pao$survcov))

  ## create design matrices file
  v=colnames(temp.pao$unitcov); i=grep("_psi",v)
  psi.dm<-matrix(v[i],3,ncol(psi_mat),byrow=TRUE)
  rownames(psi.dm)<-psi.label; colnames(psi.dm)<-paste0("a",1:ncol(psi_mat))

  psi.dm.col<-ncol(psi.dm)
  psi.dm=simplify_dm(psi.dm,pao=temp.pao); colnames(psi.dm)=paste0('a',1:ncol(psi.dm))

  v=colnames(temp.pao$survcov); i=grep("_[pr]",v)
  p.dm<-matrix(rep(matrix(v[i],nrow=5,byrow=T),each=temp.pao$nsurveys),nrow=5*temp.pao$nsurveys)
  p.dm.col<-ncol(p.dm)
  rownames(p.dm)<-paste0(rep(p.label,each=temp.pao$nsurveys),"[",1:temp.pao$nsurveys,"]")
  colnames(p.dm)<-paste0(ifelse(is.null(gamma),"b","d"),1:ncol(p_mat))
  p.dm.save=p.dm; p.dm=simplify_dm(p.dm.save,temp.pao)

  v=colnames(temp.pao$survcov); i=grep("_o[aAbB]",v)
  o.dm<-matrix(rep(matrix(v[i],nrow=4,byrow=T),each=temp.pao$nsurveys),nrow=4*temp.pao$nsurveys)
  o.dm.col<-ncol(o.dm)
  rownames(o.dm)<-paste0(rep(omeg.label,each=temp.pao$nsurveys),"[",1:temp.pao$nsurveys,"]")
  colnames(o.dm)<-paste0(ifelse(is.null(gamma),"b","e"),1:ncol(omeg_mat))
  o.dm.save=o.dm; o.dm=simplify_dm(o.dm.save,temp.pao)

  v=colnames(temp.pao$survcov); i=grep("_c[AB]",v)
  c.dm<-matrix(rep(matrix(v[i],nrow=2,byrow=T),each=temp.pao$nsurveys),nrow=2*temp.pao$nsurveys)
  c.dm.col<-ncol(c.dm)
  rownames(c.dm)<-paste0(rep(conf.label,each=temp.pao$nsurveys),"[",1:temp.pao$nsurveys,"]")
  colnames(c.dm)<-paste0(ifelse(is.null(gamma),"b","f"),1:ncol(conf_mat))
  c.dm.save=c.dm; c.dm=simplify_dm(c.dm.save,temp.pao)

  tmp.unitcov=NULL; temp.pao$nunitcov=0; lst=get_not01_dm(c(psi.dm),names(temp.pao$unitcov))
  if (sum(lst)>0) {
    for (i in 1:length(lst)) if (lst[i]>0) tmp.unitcov=cbind(tmp.unitcov,temp.pao$unitcov[,i])
    colnames(tmp.unitcov)=colnames(temp.pao$unitcov)[which(lst>0)]
  }
  temp.pao$unitcov=tmp.unitcov; if (!is.null(tmp.unitcov)) temp.pao$nunitcov=ncol(tmp.unitcov);

  tmp.survcov=NULL; temp.pao$nsurvcov=0; lst=get_not01_dm(c(psi.dm),names(temp.pao$survcov))
  if (sum(lst)>0) {
    for (i in 1:length(lst)) if (lst[i]>0) tmp.survcov=cbind(tmp.survcov,temp.pao$survcov[,i])
    colnames(tmp.survcov)=colnames(temp.pao$survcov)[which(lst>0)]
  }
  temp.pao$survcov=tmp.survcov; if (!is.null(tmp.survcov)) temp.pao$nsurvcov=ncol(tmp.survcov)

  #write_pao(temp.pao)
  modname=gsub(' ','',paste0('psi(',as.character(psi),')p(',as.character(p),
                             ')omeg(',as.character(omeg),')conf(',as.character(omeg),')')[2])

  if(!is.null(fixed)){
    fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm,o.dm,c.dm),rownames)))
  }
  rv<-runPresence(temp.pao,list(psi.dm,p.dm,NULL,NULL,o.dm,c.dm),
                  model,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  npar<-psi.dm.col+p.dm.col+o.dm.col+c.dm.col;

  ### extract results from output file
   v=readLines(outfile); l=length(v);
   ## find AIC value for each model
   instns=grep("-2log(likelihood)",v,fixed=TRUE)
   z<-get_line(v,instns[1])
   neg2loglike<-as.numeric(z[4])
   z<-get_line(v,instns[1]+1)
   aic<-as.numeric(z[3]);   cat('aic=',aic,'\n')

   ### get_coefficients
  instns<-grep("Untransformed",v,fixed=TRUE)
  extract<-v[(instns[1]+3):(instns[1]+2+npar)]
  index<-grep("^A.+ psi",extract)
  psi.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(psi_mat),
                       cov.names=gsub('Intercept','Int',colnames(psi_mat)),b.string="a",b.count=0)
  index<-grep("^B",extract);
  p.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(p_mat),
                     cov.names=gsub('Intercept','Int',colnames(p_mat)),b.string="b",b.count=0)
  index<-grep("^E",extract)
  omeg.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(omeg_mat),
                        cov.names=gsub('Intercept','Int',colnames(omeg_mat)),b.string="e",b.count=0)
  index<-grep("^F",extract)
  conf.coeff<-get_coeff(v=extract,index=index,n.coeff=ncol(conf_mat),
                        cov.names=gsub('Intercept','Int',colnames(conf_mat)),b.string="e",b.count=0)
  names<-c(rownames(psi.coeff),rownames(p.coeff),rownames(omeg.coeff),rownames(conf.coeff))

  ### get_VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC<-get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC<-VC[psi.coeff[,"index"],psi.coeff[,"index"]]
  p.VC<-VC[p.coeff[,"index"],p.coeff[,"index"]]
  omeg.VC=VC[omeg.coeff[,"index"],omeg.coeff[,"index"]]
  conf.VC=VC[conf.coeff[,"index"],conf.coeff[,"index"]]

  ## Get real parameter estimates from output
  i=grep("parameter +site +estimate",v)[1]; v=v[i:length(v)]
  parse_est <- function(parm,v) {
    i=grep(paste0('^',parm),v); vv=v[i];
    v2=gsub('D0-','D0 -',gsub('^ ','',gsub(' +',' ',gsub('.+:','',gsub(' - ','',vv)))))
    v3=matrix(as.numeric(unlist(strsplit(v2,' '))),ncol=4,byrow=T)
    rownames(v3)=gsub(' .+','',vv); colnames(v3)=c('est','se','lowCI','hiCI')
    return(v3)
  }
  psi.est=parse_est('psi',v)
  i=grep("parameter +site +estimate",v)[2]; v=v[i:length(v)]
  p.est<-parse_est('[pr][AB]',v)
  omeg.est=parse_est('o[aAbB]',v)
  conf.est=parse_est('c[AB]',v)

  ##### check for warnings
  warn.conv<-check_conv_warn(v)
  warn.VC<-check_VC_warn(v)

   result<-list(modname=modname,
                model=list(psi=psi,p=p),dmat=list(psi=psi.dm,p=p.dm,o=o.dm,c=c.dm),
                data=temp.pao,outfile=outfile,
                neg2loglike=neg2loglike,
                npar=npar, aic=aic,
                beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                          p=as.data.frame(p.coeff),p.VC=p.VC,
                          omeg=as.data.frame(omeg.coeff),omeg.VC=omeg.VC,VC=VC,
                          conf=as.data.frame(conf.coeff),conf.VC=conf.VC,VC=VC),
                real=list(psi=as.data.frame(psi.est),
                          p=as.data.frame(p.est),
                          omeg=as.data.frame(omeg.est),
                          conf=as.data.frame(conf.est)),
                warnings=list(conv=warn.conv,VC=warn.VC))
   class(result)<-c("occMod","so2spfp")
   return(result)
}
