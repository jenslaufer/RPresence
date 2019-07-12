occMod_SO_2SP_cd<-function(psi=call(),psi.cov=data$unitcov,
                         theta=call(), theta.cov=cbind(data$unitcov,data$survcov),
                         p=call(),p.cov=cbind(data$unitcov,data$survcov),
                         th0pi=call(), th0pi.cov=data$unitcov,
                         param="nu", sp.contr=TRUE,modname=NULL,paoname=NULL,outfile,model=3000,
                         fixed=NULL,initvals=NULL,data,miscopts=''){

  #' Fit two species, static occupancy (single season) model with correlated detectoins
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
  #' @param theta The right-hand side of the formula for the model to fit for "local-use" probability. The
  #' following terms are pre-defined covariates to allow constraints on the theta's (thA, thA\', thBA, thBA\',
  #' thBa, thBa\').
  #' \itemize{
  #' \item{\code{SP} - species effect on theta's (thA != thB), }
  #' \item{\code{PRIME} - "PRIME" effect (thA !=thA\', thB. !=thB.\') }
  #' \item{\code{BAa} - species interaction on theta (thBA != thBa, thBA\' != thBa\') }
  #' \item{\code{INTth} - interaction of \code{PRIME} and \code{BAa} (thBA != thBA\' and thBA != thBa)}
  #' }
  #' So, the most generl model (all theta\'s different) would be: theta(\code{~SP*PRIME+BAa+INTth)}
  #' @param theta.cov A data frame containing the survey-specific covariates to use for the theta component
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
  #' @param th0pi right-side of formula for model to fit for init proportion prev. unoccupied.
  #' @param th0pi.cov data frame containing the site-specific covariates to use for th0pi.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param model The PRESENCE model code. DO NOT CHANGE.
  #' @param fixed A single-column matrix containing values for real parameters to be fixed at.
  #' \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals Initial values for the beta parameters at which PRESENCE begins the optimisation.
  #' The default values in PRESENCE is 0.
  #' @param data The \code{pao} data object containing the detection data and other information.
  #' @param miscopts (see \code{\link{occMod}})

  #' NOTE THAT THERE MAY BE SOME CHANGES TO HOW THIS MODEL IS IMPLEMENTED IN THE NEAR FUTURE.
  #'
  #' @return returns a list of class \code{occMod} and \code{so2spCd}
  #'
  #' \code{occMod$beta} contains the objects:
  #'  \item{psi}{estimated logistic regression coefficients and std.err for prob. of occurrence.}
  #'  \item{psi.VC}{var-covar matrix for logistic regression coefficients for prob. of occurrence.}
  #'  \item{theta}{estimated logistic regression coefficients and std.err for prob. of local occurrence in each survey.}
  #'  \item{theta.VC}{var-covar matrix for theta.}
  #'  \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'  \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'  \item{th0pi}{estimated logistic regression coefficients and std.err for prob. of local occurrence before 1st survey.}
  #'  \item{th0pi.VC}{var-covar matrix for th0pi.}
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
  if(!(param%in%c("psiBA","nu"))){
    stop("Invalid parameterisation. Value of param must be either 'psiBA' or 'nu'!")
  }
  multiseason=F
  options("na.action"="na.pass")

  if (sum(data$det.data>1)==0) {  # if no 2's or 3's, then stacked input...
    data$nunits=data$nunits/2     #     convert to compressed input... 0=neither, 1=sp.A, 2=sp.B, 3=both
    data$det.data=data$det.data[1:data$nunits,]+2*data$det.data[data$nunits+1:data$nunits,]
  }
  samp.units<-data$nunits; data$unitnames=paste0('unit',1:samp.units)

  if (is.null(psi.cov)) psi.cov=data.frame(SP=as.factor(c(rep("A",samp.units),rep("B",2*samp.units)))) else {
    if (nrow(psi.cov)==samp.units) psi.cov=rbind(psi.cov,psi.cov,psi.cov)
    psi.cov$SP<-as.factor(c(rep("A",samp.units),rep("B",2*samp.units)))
  }
  contrasts(psi.cov$SP)<-contr.treatment(2,contrasts=sp.contr)
  psi.cov$INT<-as.factor(c(rep(1,2*samp.units),rep(2,samp.units)))  #  INT - psiBA,psiBa interaction
  psi_mat<-model.matrix(psi,psi.cov)


  ns=samp.units*data$nsurveys
  if(is.null(theta.cov)) theta.cov<-data.frame(SURVEY=rep(1:data$nsurveys,each=samp.units)) else
    theta.cov$SURVEY=rep(1:data$nsurveys,each=samp.units)
  theta.cov=data.frame(apply(theta.cov,2,rep,6))
  theta.cov$SP<-as.factor(c(rep("A",2*ns),rep("B",4*ns)))
  theta.cov$PRIME<-as.factor(rep(c(rep(0,ns),rep(1,ns)),3))
  theta.cov$BAa <- as.factor(c(rep(0,4*ns),rep(1,2*ns)));
  theta.cov$INTth=as.factor(c(rep(0,5*ns),rep(1,ns)))
  theta.cov$ALLDIFF=as.factor(rep(1:6,each=ns))
  th_mat<-model.matrix(theta,theta.cov);   # full model: th=~SP*PRIME+BAa;

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
  p.cov$INT_so<-as.factor(c(rep(1,3*len),rep(2,2*len)))        #  interaction: rB != rA
  p.cov$ALLDIFF<-as.factor(rep(1:5,each=len))
  p_mat<-model.matrix(p,p.cov);

  #   gamma and epsilon built-in covariates
  if (multiseason) {
    SP<-as.factor(c(rep("A",2*samp.units),rep("B",4*samp.units))); SP=rep(SP,data$nseasons-1)
    contrasts(SP)<-contr.treatment(2,contrasts=sp.contr)
    INT_Ba<-as.factor(c(rep(0,4*samp.units),rep(1,2*samp.units)))                  #  interaction: BA? != Ba?
    INT_B_a<-as.factor(c(rep(0,3*samp.units),rep(1,samp.units),rep(0,samp.units),rep(1,samp.units)))  #  interaction: B?A != B?a
    INT_Baa<-as.factor(c(rep(0,5*samp.units),rep(1,samp.units)))  #                             interaction: B?? != Baa
    INT_Ab<-as.factor(c(rep(0,samp.units),rep(1,samp.units),rep(0,4*samp.units)))  #     interaction: B?? != Baa
    INT_Ba=rep(INT_Ba,data$nseasons-1)
    INT_B_a=rep(INT_B_a,data$nseasons-1)
    INT_Baa=rep(INT_Baa,data$nseasons-1)
    INT_Ab=rep(INT_Ab,data$nseasons-1)
    SEASON=as.factor(rep(rep(1:(data$nseasons-1),each=samp.units),6))
    gamma.cov<-data.frame(SP,INT_Ba,INT_B_a,INT_Baa,INT_Ab,SEASON)
    gamma_mat<-model.matrix(gamma,gamma.cov)
    epsilon.cov<-data.frame(SP,INT_Ba,INT_B_a,INT_Baa,INT_Ab,SEASON)
    epsilon_mat<-model.matrix(epsilon,epsilon.cov)
  }
  psi.label<-c("psiA","psiBA","psiBa"); p.label<-c("pA","pB","rA","rBA","rBa")
  th.label<-c("thA","thA'","thBA","thBA'","thBa","thBa'")
  gamma.label=paste0("gam",c("AB","Ab","BAA","BAa","BaA","Baa"))
  epsilon.label=paste0("eps",c("AB","Ab","BAA","BAa","BaA","Baa"))

  if(param=="nu"){
    ## set terms to zero if using the nu parameterisation and fix nu/rho to 1
    c_index<-(1:ncol(psi_mat))[!grepl("INT",colnames(psi_mat),fixed=TRUE)]
    psi_mat[(2*samp.units+1):(3*samp.units),c_index]<-0

    c_index<-(1:ncol(p_mat))[!grepl("INT_d",colnames(p_mat),fixed=TRUE)]
    p_mat[(4*len+1):(5*len),c_index]<-0

    ## create labels
    psi.label<-c("psiA","psiBa","nu");    p.label<-c("pA","pB","rA","rBa","rho")
  }
  th0pi_mat<-model.matrix(th0pi,psi.cov); th0pi.label=c("th0piA","th0piBA","th0piBa")

  ## create and output temporary pao file
  ### bind each block of occupancy covariates to unitcov
  temp.cov<-cbind(psi_mat[1:samp.units,],psi_mat[samp.units+1:samp.units,],psi_mat[2*samp.units+1:samp.units,],
                th0pi_mat[1:samp.units,],th0pi_mat[samp.units+1:samp.units,],th0pi_mat[2*samp.units+1:samp.units,])
  v1=c(paste0(rep(colnames(psi_mat),3),"_",rep(psi.label,each=ncol(psi_mat))),
       paste0(rep(colnames(th0pi_mat),3),"_",rep(th0pi.label,each=ncol(th0pi_mat))))
  if (multiseason) {
    temp.cov=cbind(temp.cov,
                  gamma_mat[1:samp.units,],gamma_mat[samp.units+1:samp.units,],
                  gamma_mat[2*samp.units+1:samp.units,],gamma_mat[3*samp.units+1:samp.units,],
                  gamma_mat[4*samp.units+1:samp.units,],gamma_mat[5*samp.units+1:samp.units,],
                  epsilon_mat[1:samp.units,],              epsilon_mat[samp.units+1:samp.units,],
                  epsilon_mat[2*samp.units+1:samp.units,],epsilon_mat[3*samp.units+1:samp.units,],
                  epsilon_mat[4*samp.units+1:samp.units,],epsilon_mat[5*samp.units+1:samp.units,])
    v1=c(v1,paste0(rep(colnames(gamma_mat),6),"_",rep(gamma.label,each=ncol(gamma_mat))),
            paste0(rep(colnames(epsilon_mat),6),"_",rep(epsilon.label,each=ncol(epsilon_mat))))
  }
  colnames(temp.cov)<-v1
  temp.pao<-data; temp.pao$nunitcov<-ncol(temp.cov); temp.pao$unitcov<-as.data.frame(temp.cov)
  rownames(temp.pao$unitcov)<-data$unitnames; if (length(temp.pao$surveynames)<1) temp.pao$surveynames=1:temp.pao$nsurveys

  ## bind each block of theta covariates
  temp.cov<-cbind(th_mat[1:len,],th_mat[len+1:len,],th_mat[2*len+1:len,],th_mat[3*len+1:len,],
                  th_mat[4*len+1:len,],th_mat[5*len+1:len,])
  colnames(temp.cov)<-paste0(rep(colnames(th_mat),6),"_",rep(th.label,each=ncol(th_mat)))
  v2=paste0(rep(colnames(th_mat),6),"_",rep(th.label,each=ncol(th_mat)))

  ## bind each block of detection covariates
  temp.cov<-cbind(temp.cov,p_mat[1:len,],p_mat[len+1:len,],p_mat[2*len+1:len,],p_mat[3*len+1:len,],p_mat[4*len+1:len,])
  colnames(temp.cov)<-c(v2,paste0(rep(colnames(p_mat),5),"_",rep(p.label,each=ncol(p_mat))))

  temp.pao$nsurvcov<-ncol(temp.cov)
  temp.pao$survcov<-as.data.frame(temp.cov)
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)
  rownames(temp.pao$survcov)<-1:dim(temp.pao$survcov)[1]

  colnames(temp.pao$unitcov)<-gsub(".Intercept.","Int",colnames(temp.pao$unitcov))
  colnames(temp.pao$survcov)<-gsub(".Intercept.","Int",colnames(temp.pao$survcov))

  ## create design matrices file
  v=colnames(temp.pao$unitcov); i=grep("_psi",v)
  psi.dm<-matrix(v[i],3,ncol(psi_mat),byrow=TRUE)
  rownames(psi.dm)<-psi.label; colnames(psi.dm)<-paste0("a",1:ncol(psi_mat)); psi.dm.save=psi.dm

  i=grep("_th0pi",v)
  th0pi.dm<-simplify_dm(matrix(v[i],3,ncol(th0pi_mat),byrow=TRUE),temp.pao)
  th0pi.dm.col<-ncol(th0pi.dm); rownames(th0pi.dm)<-th0pi.label; colnames(th0pi.dm)<-paste0("c",1:ncol(th0pi_mat))

  v=colnames(temp.pao$survcov); i=grep("_th",v)
  th.dm<-matrix(rep(matrix(v[i],nrow=6,byrow=T),each=temp.pao$nsurveys),nrow=6*temp.pao$nsurveys)
  psi.dm.col<-ncol(psi.dm)+ncol(th.dm);
  rownames(th.dm)<-paste0(rep(th.label,each=temp.pao$nsurveys),"[",1:temp.pao$nsurveys,"]")
  colnames(th.dm)<-paste0("a",ncol(psi_mat)+1:ncol(th_mat))
  psi.dm=simplify_dm(diagbind(psi.dm,th.dm),temp.pao); colnames(psi.dm)=paste0('a',1:ncol(psi.dm))

  v=colnames(temp.pao$survcov); i=grep("_[pr]",v)
  p.dm<-matrix(rep(matrix(v[i],nrow=5,byrow=T),each=temp.pao$nsurveys),nrow=5*temp.pao$nsurveys)
  p.dm.col<-ncol(p.dm)
  rownames(p.dm)<-paste0(rep(p.label,each=temp.pao$nsurveys),"[",1:temp.pao$nsurveys,"]")
  colnames(p.dm)<-paste0(ifelse(is.null(gamma),"b","d"),1:ncol(p_mat))
  p.dm.save=p.dm; p.dm=simplify_dm(p.dm.save,temp.pao)

  rownames(th0pi.dm)=c("th0piA","th0piBA","th0piBa"); colnames(th0pi.dm)=paste0("c",1:ncol(th0pi.dm))

  gamma.dm.col=epsilon.dm.col=0; gamma.dm=epsilon.dm=NULL

  tmp.unitcov=NULL; temp.pao$nunitcov=0; lst=get_not01_dm(c(psi.dm,th0pi.dm),names(temp.pao$unitcov))
  if (sum(lst)>0) {
    for (i in 1:length(lst)) if (lst[i]>0) tmp.unitcov=cbind(tmp.unitcov,temp.pao$unitcov[,i])
    colnames(tmp.unitcov)=colnames(temp.pao$unitcov)[which(lst>0)]
  }
  temp.pao$unitcov=tmp.unitcov; if (!is.null(tmp.unitcov)) temp.pao$nunitcov=ncol(tmp.unitcov);

  tmp.survcov=NULL; temp.pao$nsurvcov=0; lst=get_not01_dm(c(psi.dm,th0pi.dm),names(temp.pao$survcov))
  if (sum(lst)>0) {
    for (i in 1:length(lst)) if (lst[i]>0) tmp.survcov=cbind(tmp.survcov,temp.pao$survcov[,i])
    colnames(tmp.survcov)=colnames(temp.pao$survcov)[which(lst>0)]
  }
  temp.pao$survcov=tmp.survcov; if (!is.null(tmp.survcov)) temp.pao$nsurvcov=ncol(tmp.survcov)

  #write_pao(temp.pao)
  modname=gsub(' ','',paste0('psi(',as.character(psi),')th(',as.character(theta),
                             ')p(',as.character(p),')th0pi(',as.character(th0pi),')')[2])

  if (multiseason) {
    for (j in 1:length(gamma.label)) {
      i=grep(gamma.label[j],v); gamma.dm<-rbind(gamma.dm,matrix(v[i],data$nseasons-1,ncol(gamma_mat),byrow=T))
    }
    gamma.dm.col<-ncol(gamma.dm);
    rownames(gamma.dm)<-paste0(rep(gamma.label,each=temp.pao$nseasons-1),"[",1:(temp.pao$nseasons-1),"]")
    colnames(gamma.dm)<-paste0("b",1:ncol(gamma_mat));

    for (j in 1:length(epsilon.label)) {
      i=grep(epsilon.label[j],v); epsilon.dm<-rbind(epsilon.dm,matrix(v[i],data$nseasons-1,ncol(epsilon_mat),byrow=T))
    }
    epsilon.dm.col<-ncol(epsilon.dm);
    rownames(epsilon.dm)<-paste0(rep(epsilon.label,each=temp.pao$nseasons-1),"[",1:temp.pao$nseasons,"]")
    colnames(epsilon.dm)<-paste0("c",1:ncol(epsilon_mat))

    if(!is.null(fixed)){
      fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,gamma.dm,epsilon.dm,p.dm),rownames)))
    }
    outfile<-runPresence(temp.pao,list(psi.dm,gamma.dm,epsilon.dm,p.dm,NULL,NULL),
                          model,modname=modname,fixed=fixed,initvals=initvals,miscopts=miscopts)
  } else {
    if(!is.null(fixed)){
      fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm,th0pi.dm),rownames)))
    }
    rv<-runPresence(temp.pao,list(psi.dm,p.dm,th0pi.dm,NULL,NULL,NULL),
                    model,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  }
  npar<-psi.dm.col+p.dm.col+gamma.dm.col+epsilon.dm.col+th0pi.dm.col;

  ################################################################
  #### Extract results from output file
  ################################################################
  v=readLines(outfile)

  ### get PRESENCE version
  i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

  ## find AIC value for each model
  i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))
  v=v[-1:-i]

  gam.coeff=eps.coeff=theta.coeff=th0pi.coeff=NULL
  ### get Untransformed coefficients (Beta's)
  i=grep("^A[0-9]+ +psi.+ : ",v); psi.coeff=getbetas(v[i],psi.dm)
  i2=grep("^A[0-9]+ +th.+ : ",v); th.coeff=getbetas(v[i2],th.dm)
  j=grep("^B[0-9]+ gam.+ : ",v); gam.coeff=getbetas(v[j],gamma.dm)
  k=grep("^C[0-9]+ eps.+ : ",v); eps.coeff=getbetas(v[k],epsilon.dm)
  l=grep("^[BD][0-9]+ +[pr].+ : ",v); p.coeff=getbetas(v[l],p.dm)
  m=grep("^[CE][0-9]+ +th.+ : ",v); th0pi.coeff=getbetas(v[m],th0pi.dm)

  ### get_VC matrix
  ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)
  instns=grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC=get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  ii=1:length(i); psi.VC=VC[ii,ii]
  ii=1:length(i2); n=length(i); th.VC=VC[n+ii,n+ii]
  ii=1:length(j); n=n+length(i2); if (length(j)>0) p.VC<-VC[n+ii,n+ii]
  ii=1:length(k); n=n+length(j); if (length(k)>0) p.VC<-VC[n+ii,n+ii]
  ii=1:length(l); n=n+length(k); p.VC=VC[n+ii,n+ii]
  ii=1:length(m); n=n+length(l); th0pi.VC=VC[n+ii,n+ii]

  ## Get real parameter estimates from output
  psiA.est<-get_real(real="<psiA>",v=v,con.offset=1,n1=1,nunits=samp.units,
                    row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
  psiBa.est<-get_real(real="<psiBa>",v=v,con.offset=1,n1=1,nunits=samp.units,
                   row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
  if(param=="nu") {
    nu.est<-get_real(real="<nu>",v=v,con.offset=1,n1=1,nunits=samp.units,
                    row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    ## derive psiBA
    psiBA.mat<-as.matrix(psi_mat[(samp.units+1):(2*samp.units),]+psi_mat[(2*samp.units+1):(3*samp.units),])
    logit.psiBA<-psiBA.mat%*%psi.coeff[,"est"]
    psiBA<-plogis(logit.psiBA)
    logit.var<-apply(psiBA.mat,1,function(xx) xx%*%psi.VC%*%xx)
    se<-sqrt(logit.var)*(psiBA*(1-psiBA))
    lower<-plogis(logit.psiBA-1.96*sqrt(logit.var))
    upper<-plogis(logit.psiBA+1.96*sqrt(logit.var))
    psiBA.est<-cbind(psiBA,se,lower,upper)
    colnames(psiBA.est)<-c("est","se","lower","upper")
    rownames(psiBA.est)<-temp.pao$unitnames[1:samp.units]

  } else {
    psiBA.est<-get_real(real="<psiBA>",v=v,con.offset=1,n1=1,nunits=samp.units,
                     row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    ## derive nu
    nu.mat<-as.matrix(psi_mat[(samp.units+1):(2*samp.units),]-psi_mat[(2*samp.units+1):(3*samp.units),])
    log.nu<-nu.mat%*%psi.coeff[,"est"]
    nu<-exp(log.nu)
    log.var<-apply(nu.mat,1,function(xx) xx%*%psi.VC%*%xx)
    se<-sqrt(log.var)*nu
    lower<-exp(log.nu-1.96*sqrt(log.var))
    upper<-exp(log.nu+1.96*sqrt(log.var))
    nu.est<-cbind(nu,se,lower,upper)
    colnames(nu.est)<-c("est","se","lower","upper")
    rownames(nu.est)<-temp.pao$unitnames[1:samp.units]
  }

  rnames=paste0(rep(temp.pao$unitnames[1:samp.units],temp.pao$nsurveys),"_",rep(temp.pao$surveynames,each=samp.units))
  thA.est<-get_real(real="<thA[1]>",    v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  thAp.est<-get_real(real="<thA'[1]>",  v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  thBA.est<-get_real (real="<thBA[1]>", v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  thBAp.est<-get_real(real="<thBA'[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  thBa.est<-get_real( real="<thBa[1]>", v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  thBap.est<-get_real(real="<thBa'[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  pA.est<-get_real(real="<pA[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  pB.est<-get_real(real="<pB[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  rA.est<-get_real(real="<rA[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  rBa.est<-get_real(real="<rBa[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,
                 row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
  if(param=="nu") {
    rho.est<-get_real(real="<rho[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,
                      row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
    ## derive rBA
    rBA.mat<-as.matrix(p_mat[(3*len+1):(4*len),]+p_mat[(4*len+1):(5*len),])
    logit.rBA<-rBA.mat%*%p.coeff[,"est"]
    rBA<-plogis(logit.rBA)
    p.VC<-VC[(psi.dm.col+1):npar,(psi.dm.col+1):npar]
    logit.var<-apply(rBA.mat,1,function(xx) xx%*%p.VC%*%xx)
    se<-sqrt(logit.var)*(rBA*(1-rBA))
    lower<-plogis(logit.rBA-1.96*sqrt(logit.var))
    upper<-plogis(logit.rBA+1.96*sqrt(logit.var))

    rBA.est<-cbind(rBA,se,lower,upper)
    colnames(rBA.est)<-c("est","se","lower","upper")
    rownames(rBA.est)<-rnames

  } else {
    rBA.est<-get_real(real="<rBA[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,
                      row.names=rnames,fixed=fixed,real.offset=temp.pao$nseasons+1)
    ## derive rho
    rho.mat<-as.matrix(p_mat[(3*len+1):(4*len),]-p_mat[(4*len+1):(5*len),])
    log.rho<-rho.mat%*%p.coeff[,"est"]
    rho<-exp(log.rho)
    i=grep("^B",rownames(VC)); if (multiseason) i=grep("^D",rownames(VC))
    p.VC<-VC[i,i]
    log.var<-apply(rho.mat,1,function(xx) xx%*%p.VC%*%xx)
    se<-sqrt(log.var)*rho
    lower<-exp(log.rho-1.96*sqrt(log.var))
    upper<-exp(log.rho+1.96*sqrt(log.var))

    rho.est<-cbind(rho,se,lower,upper)
    colnames(rho.est)<-c("est","se","lower","upper")
    rownames(rho.est)<-rnames
}
  v1=v[grep('th0pi.+:.+ - ',v)]
  v2=gsub('fixed','',gsub('^,','',gsub(' +',',',gsub(' - ','',gsub('.+:','',v1)))))
  th0pi.est=matrix(as.numeric(unlist(sapply(v2,strsplit,","))),ncol=4,byrow=T);
  rownames(th0pi.est)=gsub(' .+','',v1); colnames(th0pi.est)=c("est","se","lower","upper")

    if (multiseason) {
    index<-grep("^B",extract)
    gamAB.est<-get_real(real="<gamAb",v=v,con.offset=1,n1=1,nunits=samp.units,
                        row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    gamAb.est<-get_real(real="<gamAb",v=v,con.offset=1,n1=1,nunits=samp.units,
                        row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    gamBAA.est<-get_real(real="<gamBAA",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    gamBAa.est<-get_real(real="<gamBAa",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    gamBaA.est<-get_real(real="<gamBaA",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    gamBaa.est<-get_real(real="<gamBaa",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    index<-grep("^C",extract)
    epsAB.est<-get_real(real="<epsAb",v=v,con.offset=1,n1=1,nunits=samp.units,
                        row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    epsAb.est<-get_real(real="<epsAb",v=v,con.offset=1,n1=1,nunits=samp.units,
                        row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    epsBAA.est<-get_real(real="<epsBAA",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    epsBAa.est<-get_real(real="<epsBAa",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    epsBaA.est<-get_real(real="<epsBaA",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
    epsBaa.est<-get_real(real="<epsBaa",v=v,con.offset=1,n1=1,nunits=samp.units,
                         row.names=temp.pao$unitnames[1:samp.units],fixed=fixed,real.offset=0)
  }
  ##### check for warnings
  warn.conv<-check_conv_warn(v)
  warn.VC<-check_VC_warn(v)

   result<-list(modname=modname,
                model=list(psi=psi,p=p),dmat=list(psi=psi.dm,p=p.dm,th0pi=th0pi.dm),
                data=temp.pao,outfile=outfile,
                neg2loglike=neg2loglike,
                npar=npar, aic=aic,
                beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
                          theta=as.data.frame(th.coeff),th.VC=th.VC,
                          p=as.data.frame(p.coeff),p.VC=p.VC,
                          th0pi=as.data.frame(th0pi.coeff),th0pi.VC=th0pi.VC,VC=VC),
                real=list(psiA=as.data.frame(psiA.est),
                          psiBa=as.data.frame(psiBa.est), psiBA=as.data.frame(psiBA.est),
                          nu=as.data.frame(nu.est),
                          thA=as.data.frame(thA.est),thAp=as.data.frame(thAp.est),
                          thBA=as.data.frame(thBA.est),thBAp=as.data.frame(thBAp.est),
                          thBa=as.data.frame(thBa.est),thBap=as.data.frame(thBap.est),
                          pA=as.data.frame(pA.est), pB=as.data.frame(pB.est),
                          rA=as.data.frame(rA.est), rBa=as.data.frame(rBa.est),
                          rBA=as.data.frame(rBA.est), rho=as.data.frame(rho.est)),
                          th0pi=as.data.frame(th0pi.est),
                warnings=list(conv=warn.conv,VC=warn.VC))
   if (multiseason) {
     result$real$gammaAB=as.data.frame(gamAB.est)
     result$real$gammaAb=as.data.frame(gamAb.est)
     result$real$gammaBAA=as.data.frame(gamBAA.est)
     result$real$gammaBAa=as.data.frame(gamBAa.est)
     result$real$gammaBaA=as.data.frame(gamBaA.est)
     result$real$gammaBaa=as.data.frame(gamBaa.est)
     result$real$epsilonAB=as.data.frame(epsAB.est)
     result$real$epsilonAb=as.data.frame(epsAb.est)
     result$real$epsilonBAA=as.data.frame(epsBAA.est)
     result$real$epsilonBAa=as.data.frame(epsBAa.est)
     result$real$epsilonBaA=as.data.frame(epsBaA.est)
     result$real$epsilonBaa=as.data.frame(epsBaa.est)
     result$dmat$gamma=gamma.dm; result$dmat$epsilon=epsilon.dm
   }
   class(result)<-c("occMod","so2spCd")
   return(result)
}
