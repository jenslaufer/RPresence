occMod_2SP<-function(psi=call(),psi.cov=data$unitcov,
                         p=call(),p.cov=cbind(data$unitcov,data$survcov),
                         gamma=call(),gamma.cov=cbind(data$unitcov,data$survcov),
                         epsilon=call(),epsilon.cov=cbind(data$unitcov,data$survcov),
                         param="nu", sp.contr=TRUE, modname=NULL,paoname=NULL, outfile,
                  model=3000,fixed=NULL,initvals=NULL,data,miscopts=miscopts){

  #' Fit two species, static occupancy (single season), or dynamic occupancy (multi-season) model
  #'
  #' This is not intended for direct use, but instead the \code{\link{occMod}} function should be used with \code{type="so.2sp.1"} or \code{type="so.2sp.2"} depending on the parameterisation required. NOTE THAT THERE MAY BE SOME CHANGES TO HOW THIS MODEL IS IMPLEMENTED IN THE NEAR FUTURE.
  #'
  #'  Note:  This function assumes data are in "compressed" format (ie.,\cr
  #'        0=neither species detected,\cr
  #'        1=only species A detected,\cr
  #'        2=only species B detected,\cr
  #'        3=both species detected)\cr
  #'
  #'        "Stacked" format can easily be converted to "compressed" format by:\cr
  #'        \code{nsites=nrow(det.data)/2\cr
  #'        new.det.data=det.data[1:nsites,]+2*det.data[nsites+1:nsites,]}
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
  #' @param gamma The right-hand side of the formula for the model to fit for colonization probability. The
  #' following terms can be used to define effects on colonization:
  #' \itemize{
  #' \item{\code{SP}}{ species effect on colonization (gammaB != gammaA)}
  #' \item{\code{SEASON}}{ seasonal effect on colonization (gammaXX(t))}
  #' \item{\code{INT_Ba}}{ interaction of effect on colonization for species B when species A was, or was
  #' not present (gammaBA != gammaBa)}
  #' \item{\code{INT_B_a}}{ interaction of effect on colonization for species B when species A was present
  #' in the previous season and was, or was not present in succeeding season (gammaB?A != gammaB?a)}
  #' \item{\code{INT_Baa}}{ interaction of effect on colonization for species B when species A was not present
  #' in the previous season and was, or was not present in succeeding season (gammaBaA != gammaBaa)}
  #' \item{\code{INT_Ab}}{ interaction of effect on colonization for species A when species B was, or was
  #' not present in the previous season (gammaAB != gammaAb)}
  #' }
  #' These terms do not have to be
  #' defined as variables in \code{p.cov}. For single-season model, set gamma to NULL.
  #' @param gamma.cov A data frame containing the site-season-specific covariates to use for the colonization
  #' component of the model.  This should be a NxK data frame, where N=number of sites, and
  #' K=number of sites*number of seasons.
  #' @param epsilon The right-hand side of the formula for the model to fit for extinction probability.
  #' The terms \code{SP}, \code{INT_Ba}, \code{INT_B_a}, \code{INT_Baa}, \code{INT_Ab}, \code{SEASON}
  #' can be used to define effects on extinction, similar to those defined for colonization.
  #' These terms do not have to be defined as variables in \code{p.cov}.
  #' @param epsilon.cov A data frame containing the site-specific covariates to use for the extinction
  #' component of the model, similar to gaam.cov.
  #' @param param The parameterisation to be used, either "\code{psiBA}" or "\code{nu}", which relate to
  #' \code{type="so.2sp.1"} and \code{type="so.2sp.2"} in the function \code{\link{occMod}} respectively.
  #' @param sp.contr \code{TRUE} (default) or \code{FALSE}. Specifies the type of contrast used for the
  #' \code{SP} term. It it used as the \code{contrasts} input in the function \code{\link{contr.treatment}}.
  #' @param modname (optional) a string containing the model name
  #' @param paoname (optional) a string containing the filename for the temporary PRESENCE data file.
  #' @param model The PRESENCE model code. DO NOT CHANGE.
  #' @param fixed A single-column matrix containing values for real parameters to be fixed at.
  #' \code{rownnames(fixed)} should contain the index of the real parameters to be fixed.
  #' @param initvals Initial values for the beta parameters at which PRESENCE begins the optimisation.
  #' The default values in PRESENCE is 0.
  #' @param data The \code{pao} data object containing the detection data and other information.
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name

  #' @param VCoutopt option for PRESENCE computation/output of beta/real var-cov matrices
  #' @details choices for \code{VCoutopt} are:
  #' \itemize{
  #'  \item{"nose"}{only compute likelihood and beta estimates, no beta var-cov matrix or beta SE's (useful for faster model selection of big models or simulations)}
  #'  \item{"betavc"}{compute but don't print beta var-cov matrix so we get beta SE's, don't compute real params}
  #'  \item{"noreal"}{compute and print beta beta var-cov matrix, don't compute real params (I think this is what you want for RPresence)}
  #'  \item{"(default)"}{compute real params and real var-cov matrix, don't print either var-cov matrix.}
  #'  \item{"realbetavc"}{compute real params and real var-cov matrix, print only beta var-cov matrix.}
  #'  \item{"bothvc"}{print both var-cov matrices.}
  #'}
  #' @param noDerived if TRUE, doesn't print derived estimates from model
  #' @param randinit number of random init value vectors to try for optimization
  #'
  #' @return returns a list of class \code{"occMod"} and \code{"do.1"}.
  #'
  #' \code{occMod$beta} contains the objects:
  #'  \item{psi}{estimated logistic regression coefficients and standard errors for probability of
  #'  occurrence in the first year.}
  #'  \item{psi.VC}{variance-covariance matrix for logistic regression coefficients for probability
  #'  of occurrence in the first year.}
  #'  \item{p}{estimated logistic regression coefficients and standard errors for probability of detection.}
  #'  \item{p.VC}{variance-covariance matrix for logistic regression coefficients for probability of detection.}
  #'  \item{gamma}{estimated logistic regression coefficients and standard errors for probability of colonization.}
  #'  \item{gamma.VC}{var-cov matrix for logistic regression coefficients for probability of colonization.}
  #'  \item{epsilon}{estimated logistic regression coefficients and standard errors for probability of extinction.}
  #'  \item{epsilon.VC}{var-cov matrix for logistic regression coefficients for probability of extinction.}
  #'  \item{VC}{the full variance-covariance matrix for all logistic regression coefficients.}
  #'
  #'  \code{occMod$real} contains the objects:
  #'  \item{psi}{estimated probabilities of occurrence for each sampling unit, along with standard errors
  #'  and limits of 95\% confidence interval. Estimates are provided for the first season (calculated
  #'  directly from the estimated \code{beta} parameters). Estimates for later seasons are provided as
  #'  derived parameters. The season for which an estimate applies to can be identified from the final
  #'  number of the row names (\code{rownames(psi)}).}
  #'  \item{p}{estimated probabilities of detection for each survey, along with standard errors and limits
  #'  of 95\% confidence interval.}
  #'  \item{gamma}{estimated probabilities of colonization for each sampling unit and season, along with standard errors and limits
  #'  of 95\% confidence interval.}
  #'  \item{epsilon}{estimated probabilities of extinction for each sampling unit and season, along with standard errors and limits
  #'  of 95\% confidence interval.}
  #'
  #'  \code{occMod$derived} contains the object:
  #'  \item{psi}{estimated probabilities of occurrence for each sampling unit for second season onwards,
  #'  along with standard errors and limits of 95\% confidence interval. The season for which an estimate
  #'  applies to can be identified from the final number of the row names (\code{rownames(psi)}).}
  #' @author Darryl MacKenzie
  #' @seealso \code{\link{occMod}}
  #' @examples
  #'# load a csv data file
  #'filename<-system.file("extdata/twosp_exmpl.csv",package="RPresence")
  #'csv<-read.csv(filename,as.is=T,header=F)
  #'##
  #'nsites=nrow(dethist); nsrvys=ncol(dethist)         #  set number of sites,surveys from det. history data
  #'dethist=matrix(as.integer(unlist(dethist)),nrow=nsites) # replace missing values (-) with NA
  #'cov1=cov2=NULL
  #'##          create input "pao" object, for use with occMod function
  #'data=createPao(dethist,unitcov=cov1,survcov=cov2,title="twosp example")

  #'
  #'## fit some models
  #'##    occupancy: species-specific, no interaction, parameters: psiA, psiBA=psiBa
  #'##    detection: species-specific, p=r
  #'mod1<-occMod(model=list(psi~SP,p~SP),data=data,type="so.2sp.1",param="PsiBA")
  #'
  #'##    occupancy: species-specific, interaction, parameters: psiA, psiBA, psiBa
  #'##    detection: species-specific, p=r
  #'mod2<-occMod(model=list(psi~SP+INT,p~SP),data=data,type="so.2sp.1",param="PsiBA")
  #'
  #'##    occupancy: species-specific, interaction,
  #'##    detection: species=specific, interaction p,r (INT_o) and rBA,rBA (INT_d)
  #'##     Note: this is default 2 sp. model (psiA,psiBA,psiBa,pA,pB,rA,rBA,rBa)
  #'mod3<-occMod(model=list(psi~SP+INT,p~SP+INT_o+INT_d+SP:INT_o),data=data,type="so.2sp.1",param="PsiBA")
  #'#
  #' \dontrun{
  #' m1<-occMod(model=list(psi~SP,p~SP,data=data,type="so.2sp.1",param="PsiBA")
  #'
  #' ## fit some multi-season, 2-species models (using different data)
  #' mod0<-occMod(model=list(psi=psi~SP+INT,                                   #  PsiA,psiB,psiBa
  #'                          gamma =gamma~SP+INT_Ab+INT_Ba+INT_B_a+INT_Baa,    #  gamAB,gamAb,gamBAA,gamBAa,gamBaA,gamBaa
  #'                          epsilon=epsilon~SP+INT_Ab+INT_Ba+INT_B_a+INT_Baa, #  epsAB,epsAb,epsBAA,epsBAa,epsBaA,epsBaa
  #'                          p      =p~SP+INT_o+INT_d+SP:INT_o),               #  pA,pB,rA,rBA,rBa
  #'               data=data,type="do.2sp.1")
  #' mod1<-occMod(model=list(psi=psi~SP+INT,                                   #  PsiA,psiB,psiBa
  #'                          gamma =gamma~SP+INT_Ab+INT_Ba,                    #  gamAB,gamAb,gamBA.,gamBa.
  #'                          epsilon=epsilon~SP+INT_Ab+INT_Ba,                 #  epsAB,epsAb,epsBA.,epsBa.
  #'                          p      =p~SP+INT_o+INT_d+SP:INT_o),               #  pA,pB,rA,rBA,rBa
  #'               data=data,type="do.2sp.1")
  #' mod2<-occMod(model=list(psi=psi~SP+INT,                                   #  PsiA,psiB,psiBa
  #'                          gamma =gamma~SP,                                  #  gamA.,gamB..
  #'                          epsilon=epsilon~SP,                               #  epsA.,epsB..
  #'                          p      =p~SP+INT_o+INT_d+SP:INT_o),               #  pA,pB,rA,rBA,rBa
  #'               data=data,type="do.2sp.1")
  #'
  #' tbl=create_aic_table(list(mod0,mod1,mod2))
  #' print(tbl$table)
  #' }

  if(!(param%in%c("psiBA","nu"))){
    stop("Invalid parameterisation. Value of param must be either 'psiBA' or 'nu'!")
  }
  multiseason=(!is.null(gamma))
  #modname<-paste0(modname,",",param);  #  append parameterization (psiBA or nu) to modelname
  options("na.action"="na.pass")
  samp.units=data$nunits
  if (sum(data$det.data>1,na.rm=T)==0) {   #  if no 2's or 3's, convert "stacked" data to "compressed"
    data$nunits=data$nunits/2; samp.units<-data$nunits;
    data$unitnames=gsub("_[AB]","",data$unitnames[1:samp.units])
    ###  compress 2 sets of sites into 1...
    n=nrow(data$det.data)/2
    det1=data$det.data[1:n,]
    det2=data$det.data[n+1:n,]
    ii=(!is.na(det1) & !is.na(det2)); det1[ii]=det1[ii]+2*det2[ii]  # = spA + spB*2 = 0, 1, 2 or 3
    ii=(!is.na(det1) & is.na(det2)); det1[ii]=4+det1[ii]            # spB missing, det=4(spA=0) or 5(spA=1)
    ii=(is.na(det1) & !is.na(det2)); det1[ii]=6+det2[ii]            # spA missing, det=6(spB=0) or 7(spB=1)
    data$det.data=det1
    data$unitcov=data$unitcov[1:data$nunits,]
    i=nrow(data$survcov); data$survcov=data$survcov[1:(i/2),]
  }
  if(is.null(psi.cov)) { psi.cov<-data.frame(rep(1,samp.units)) }

  SP<-as.factor(c(rep("A",samp.units),rep("B",2*samp.units)))
  contrasts(SP)<-contr.treatment(2,contrasts=sp.contr)
  INT<-as.factor(c(rep(1,2*samp.units),rep(2,samp.units)))  #  INT - psiBA,psiBa interaction
  cov.names<-names(psi.cov) ## need to respecify names in psi.cov below if current psi.cov has only 1 column
  psi.cov<-data.frame(SP,INT,psi.cov[rep(1:samp.units,3),])
  names(psi.cov)[3:ncol(psi.cov)]<-cov.names
  psi_mat<-model.matrix(psi,psi.cov)

  len<-samp.units*data$nsurveys ## length of covariates in detection design matrix
  SP<-as.factor(c(rep("A",len), rep("B",len), rep("A",len), rep("B",2*len)))
  seasn=time=NULL; for (i in 1:data$nseasons) seasn=c(seasn,rep(i,samp.units*data$nsurveyseason[i]))
  for (i in 1:data$nsurveys) time=c(time,rep(i,samp.units))
  SEASON<-as.factor(rep(seasn,5)); SURVEY=as.factor(rep(time,5))
  contrasts(SP)<-contr.treatment(2,contrasts=sp.contr)
  INT_o<-as.factor(c(rep('1',2*len),rep('2',3*len)))  #  interaction: r != p
  INT_d<-as.factor(c(rep(1,4*len),rep(2,len)))        #  interaction: rBa != rBA
  if(is.null(p.cov)) p.cov<-data.frame(SP,INT_o,INT_d,SEASON,SURVEY) else {
    cov.names<-names(p.cov)
    if (nrow(p.cov)!=len) p.cov=data.frame(SP,INT_o,INT_d,SEASON,SURVEY,p.cov[rep(1:len,5),])
    else p.cov=data.frame(SP,INT_o,INT_d,SEASON,SURVEY,p.cov)
    names(p.cov)[6:ncol(p.cov)]<-cov.names
  }
  p_mat<-model.matrix(p,p.cov)

  #   gamma and epsilon built-in covariates
  len=data$nunits*(data$nseasons-1)
  if (multiseason) {
    SP<-as.factor(c(rep("A",2*len),rep("B",4*len)))
    contrasts(SP)<-contr.treatment(2,contrasts=sp.contr)
    INT_Ba<-as.factor(c(rep(0,4*len),rep(1,2*len)))                  #  interaction: BA? != Ba?
    INT_B_a<-as.factor(c(rep(0,3*len),rep(1,len),rep(0,len),rep(1,len)))  #  interaction: B?A != B?a
    INT_Baa<-as.factor(c(rep(0,5*len),rep(1,len)))  #                             interaction: B?? != Baa
    INT_Ab<-as.factor(c(rep(0,len),rep(1,len),rep(0,4*len)))  #     interaction: B?? != Baa
    SEASON=rep(1:(data$nseasons-1),each=data$nunits)
    SEASON=as.factor(rep(SEASON,6))  #  Replicate for the 6 parameters
    if(is.null(gamma.cov)) gamma.cov<-data.frame(SP,INT_Ba,INT_B_a,INT_Baa,INT_Ab,SEASON) else {
      cov.names<-names(gamma.cov)                              #  duplicate site covariates for the 6 parameters
      tmpcov=NULL; for (i in 1:(6*(data$nseasons-1))) tmpcov=rbind(tmpcov,gamma.cov)   #  and nseasons-1 seasons
      gamma.cov=data.frame(SP,INT_Ba,INT_B_a,INT_Baa,INT_Ab,SEASON,tmpcov)
      names(gamma.cov)[-1:-6]<-cov.names
    }
    if(is.null(epsilon.cov)) epsilon.cov<-data.frame(SP,INT_Ba,INT_B_a,INT_Baa,INT_Ab,SEASON) else {
      cov.names<-names(epsilon.cov)                              #  duplicate site covariates for the 6 parameters
      tmpcov=NULL; for (i in 1:(6*(data$nseasons-1))) tmpcov=rbind(tmpcov,epsilon.cov)   #  and nseasons-1 seasons
      epsilon.cov=data.frame(SP,INT_Ba,INT_B_a,INT_Baa,INT_Ab,SEASON,tmpcov)
      names(epsilon.cov)[-1:-6]<-cov.names
    }
    gamma_mat<-model.matrix(gamma,gamma.cov)
    epsilon_mat<-model.matrix(epsilon,epsilon.cov)
  }
  psi.label<-c("psiA","psiBA","psiBa"); p.label<-c("pA","pB","rA","rBA","rBa")
  gamma.label=paste0("gam",c("AB","Ab","BAA","BAa","BaA","Baa"))
  epsilon.label=paste0("eps",c("AB","Ab","BAA","BAa","BaA","Baa"))

  len<-samp.units*data$nsurveys ## length of covariates in detection design matrix
  if (is.null(fixed)) fixed=data.frame()
  if(param=="nu"){
    ## set terms to zero if using the nu parameterisation and fix nu/rho to 1
    ii=grepl("INT",colnames(psi_mat),fixed=TRUE)
    c_index<-(1:ncol(psi_mat))[!ii]
    psi_mat[2*samp.units+1:samp.units,c_index]<-0
    if (sum(ii)==0)
      if (sum(fixed$parm=="nu")==0) fixed=rbind(fixed,data.frame(param="nu",value=1))

    ii=grepl("INT_d",colnames(psi_mat),fixed=TRUE)
    c_index<-(1:ncol(p_mat))[!ii]
    p_mat[4*len+1:len,c_index]<-0
    if (sum(ii)==0)
      if (sum(grep("rho",fixed$param))==0)
        fixed=rbind(fixed,data.frame(param=paste0("rho[",1:data$nsurveys,']'),value=rep(1,data$nsurveys)))

    ## create labels
    psi.label<-c("psiA","psiBa","nu");    p.label<-c("pA","pB","rA","rBa","rho")
  }

  ## create and output temporary pao file
  ### bind each block of occupancy covariates to unitcov
  temp.cov<-cbind(psi_mat[1:samp.units,],
                  psi_mat[samp.units+1:samp.units,],
                  psi_mat[2*samp.units+1:samp.units,])
  v1=paste0(rep(colnames(psi_mat),3),"_",rep(psi.label,each=ncol(psi_mat)))

  if (multiseason) {
    k=0
    for (i in 0:5)
      for (j in 1:(data$nseasons-1)) {
           temp.cov=cbind(temp.cov,gamma_mat[k*samp.units+1:samp.units,]); k=k+1
           v1=c(v1,paste0(colnames(gamma_mat),"_",gamma.label[i+1],"_",j))
      }
    k=0
    for (i in 0:5)
      for (j in 1:(data$nseasons-1)) {
        temp.cov=cbind(temp.cov,epsilon_mat[k*samp.units+1:samp.units,]); k=k+1
        v1=c(v1,paste0(colnames(epsilon_mat),"_",epsilon.label[i+1],"_",j))
      }
  }
  colnames(temp.cov)<-v1
  temp.pao<-data
  temp.pao$nunitcov<-ncol(temp.cov)
  temp.pao$unitcov<-as.data.frame(temp.cov)  #  temp.pao$unitcov<-as.data.frame(rbind(temp.cov,temp.cov))
  rownames(temp.pao$unitcov)<-data$unitnames

  ## bind each block of detection covariates
  temp.cov=NULL; len=data$nunits*data$nsurveys; for (i in 0:4) temp.cov=cbind(temp.cov,p_mat[i*len+1:len,])
  colnames(temp.cov)<-paste0(rep(colnames(p_mat),5),"_",rep(p.label,each=ncol(p_mat)))

  temp.pao$nsurvcov<-ncol(temp.cov)
  temp.pao$survcov=as.data.frame(temp.cov)   # temp.pao$survcov<-as.data.frame(rbind(temp.cov,temp.cov))
  temp.pao$paoname<-ifelse(is.null(paoname),"paodata.pao",paoname)
  #  temp.pao$paoname<-"test.pao"
  rownames(temp.pao$survcov)<-1:dim(temp.pao$survcov)[1]

  colnames(temp.pao$unitcov)<-gsub(".Intercept.","Int",colnames(temp.pao$unitcov))
  #colnames(temp.pao$unitcov)<-paste0("psi.",colnames(temp.pao$unitcov))
  colnames(temp.pao$survcov)<-gsub(".Intercept.","Int",colnames(temp.pao$survcov))
  #colnames(temp.pao$survcov)<-paste0("p.",colnames(temp.pao$survcov))

  ## create design matrices file
  v=colnames(temp.pao$unitcov); i=c(grep("_psi",v),grep("_nu",v))
  psi.dm<-matrix(v[i],3,ncol(psi_mat),byrow=TRUE)
  psi.dm.col<-ncol(psi.dm); rownames(psi.dm)<-psi.label; colnames(psi.dm)<-paste0("a",1:ncol(psi_mat))

  temp.names<-unlist(lapply(1:5,function(xx){
    rep(colnames(temp.pao$survcov)[((xx-1)*ncol(p_mat)+1):(xx*ncol(p_mat))],temp.pao$nsurveys)
  }))

  p.dm<-matrix(temp.names,nrow=5*temp.pao$nsurveys,ncol=ncol(p_mat),byrow=TRUE)
  p.dm.col<-ncol(p.dm)
  rownames(p.dm)<-paste0(rep(p.label,each=temp.pao$nsurveys),"[",1:temp.pao$nsurveys,"]")
  colnames(p.dm)<-paste0(ifelse(is.null(gamma),"b","d"),1:ncol(p_mat))

  gamma.dm.col=epsilon.dm.col=0; gamma.dm=epsilon.dm=NULL
  if (multiseason) {
    for (j in 1:length(gamma.label)) {
      i=grep(gamma.label[j],v);
      if (length(i)>0) gamma.dm<-rbind(gamma.dm,matrix(v[i],data$nseasons-1,ncol(gamma_mat),byrow=T))
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
      #   get rid of covariates which are always 0.0 or always 1.0
    izero=which(colSums(temp.pao$unitcov)==0); ione=which(colSums(temp.pao$unitcov!=1)==0)
    jzero=which(colSums(temp.pao$survcov)==0); jone=which(colSums(temp.pao$survcov!=1)==0)
    szero=c(names(izero),names(jzero)); sone=c(names(ione),names(jone))
    psi.dm[psi.dm %in% szero]="0"; psi.dm[psi.dm %in% sone]="1"
    p.dm[p.dm %in% szero]="0"; p.dm[p.dm %in% sone]="1"
    gamma.dm[gamma.dm %in% szero]="0"; gamma.dm[gamma.dm %in% sone]="1"
    epsilon.dm[epsilon.dm %in% szero]="0"; epsilon.dm[epsilon.dm %in% sone]="1"

    s=c(psi.dm[nchar(psi.dm)>1],p.dm[nchar(p.dm)>1],gamma.dm[nchar(gamma.dm)>1],epsilon.dm[nchar(epsilon.dm)>1])
    temp.pao$unitcov=temp.pao$unitcov[,colnames(temp.pao$unitcov) %in% s];
    temp.pao$nunitcov=ncol(temp.pao$unitcov); if (temp.pao$nunitcov<1) temp.pao$unitcov=NULL
    temp.pao$survcov=temp.pao$survcov[,colnames(temp.pao$survcov) %in% s];
    temp.pao$nsurvcov=ncol(temp.pao$survcov); if (temp.pao$nsurvcov<1) temp.pao$survcov=NULL

    psi.dm=simplify_dm(psi.dm,temp.pao);  p.dm=simplify_dm(p.dm,temp.pao)
    npar<-psi.dm.col+p.dm.col+gamma.dm.col+epsilon.dm.col
    if(!is.null(fixed)){
      fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm,gamma.dm,epsilon.dm),rownames)))
    }
    #writePao(temp.pao)
    rv<-runPresence(temp.pao,list(psi.dm,gamma.dm,epsilon.dm,p.dm,NULL,NULL),
                    model,modname=modname,fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  } else {             #  not multi-season
    npar<-psi.dm.col+p.dm.col
    if(!is.null(fixed)){
      fixed$idx<-match(fixed$param,unlist(lapply(list(psi.dm,p.dm),rownames)))
    }
    #writePao(temp.pao)
    psi.dm=simplify_dm(psi.dm,temp.pao);  p.dm=simplify_dm(p.dm,temp.pao)
    rv<-runPresence(temp.pao,list(psi.dm,p.dm,NULL,NULL,NULL,NULL),model,modname=modname,
                         fixed=fixed,initvals=initvals,outfile=outfile,miscopts=miscopts)
  }

  ### extract results from output file
   v=readLines(outfile)

   ### get PRESENCE version
   i=grep("Version",v,fixed=TRUE)[1]; version=gsub(".+Version ","",v[i])

   ## find AIC value for each model
   i=grep("-2log.likelihood",v); neg2loglike=as.numeric(gsub(".+= ","",v[i])); aic=as.numeric(gsub(".+= ","",v[i+1]))

   ### get Untransformed coefficients (Beta's)

   ### get Untransformed coefficients (Beta's)
   i=grep("^A.+ : ",v); psi.coeff=matrix(as.numeric(unlist(strsplit(gsub(" +",",",gsub(".+: +","",v[i])),","))),ncol=2,byrow=T)
   rownames(psi.coeff)=gsub(" +","_",gsub(" +:.+","",v[i]))
   #i=grep("^A.+psi.+ : ",v); psi.coeff=as.numeric(gsub(" +\\S+$","",gsub(".+: +","",v[i])))
   j=grep("^B.+ : ",v); gam.coeff=matrix(as.numeric(unlist(strsplit(gsub(" +",",",gsub(".+: +","",v[j])),","))),ncol=2,byrow=T)
   rownames(gam.coeff)=gsub(" +","_",gsub(" +:.+","",v[j]))
   colnames(psi.coeff)=colnames(gam.coeff)=c('est','se')
   if (!multiseason) p.coeff=gam.coeff else {
     k=grep("^C.+ : ",v); eps.coeff=matrix(as.numeric(unlist(strsplit(gsub(" +",",",gsub(".+: +","",v[k])),","))),ncol=2,byrow=T)
     colnames(eps.coeff)=c('est','se')
     l=grep("^D.+ : ",v); p.coeff=matrix(as.numeric(unlist(strsplit(gsub(" +",",",gsub(".+: +","",v[l])),","))),ncol=2,byrow=T)
   }
   ii=grep("^[A-F].+ : ",v); names=gsub(" .+","",v[ii]); npar=length(names)



    ### get_VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  VC<-get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC=VC[1:length(i), 1:length(i)]; p.VC=VC[length(i)+1:length(j),length(i)+1:length(j)]

  gam.VC=eps.VC=NULL
  if (multiseason) {
    gam.VC<-p.VC; eps.VC<-VC[length(c(i,j))+1:length(k),length(c(i,j))+1:length(k)]
    p.VC=VC[length(c(i,j,k))+1:length(l),length(c(i,j,k))+1:length(l)]
  }

  ## Get real parameter estimates from output
  rnames=temp.pao$unitnames[1:samp.units]
  psiA.est<-get_real(real="<psiA>",v=v,con.offset=1,n1=1,nunits=samp.units,row.names=rnames,fixed=fixed)
  psiBa.est<-get_real(real="<psiBa>",v=v,con.offset=1,n1=1,nunits=samp.units,row.names=rnames,fixed=fixed)
  if(param=="nu") {
    nu.est<-get_real(real="<nu>",v=v,con.offset=1,n1=1,nunits=samp.units,row.names=rnames,fixed=fixed)
    ## derive psiBA
    psiBA.mat<-as.matrix(psi_mat[(samp.units+1):(2*samp.units),]+psi_mat[(2*samp.units+1):(3*samp.units),])
    logit.psiBA<-psiBA.mat%*%psi.coeff[,1]
    psiBA<-plogis(logit.psiBA)
    logit.var<-apply(psiBA.mat,1,function(xx) xx%*%psi.VC%*%xx)
    se<-sqrt(logit.var)*(psiBA*(1-psiBA))
    lower<-plogis(logit.psiBA-1.96*sqrt(logit.var))
    upper<-plogis(logit.psiBA+1.96*sqrt(logit.var))

    psiBA.est<-cbind(psiBA,se,lower,upper)
    colnames(psiBA.est)<-c("est","se","lower","upper")
    rownames(psiBA.est)<-temp.pao$unitnames[1:samp.units]
  } else {
    psiBA.est<-get_real(real="<psiBA>",v=v,con.offset=1,n1=1,nunits=samp.units,row.names=rnames,fixed=fixed)
    ## derive nu
    nu.mat<-as.matrix(psi_mat[(samp.units+1):(2*samp.units),]-psi_mat[(2*samp.units+1):(3*samp.units),])
    log.nu<-nu.mat %*% psi.coeff[,1]
    nu<-exp(log.nu)
    log.var<-apply(nu.mat,1,function(xx) xx%*%psi.VC%*%xx)
    se<-sqrt(log.var)*nu
    lower<-exp(log.nu-1.96*sqrt(log.var))
    upper<-exp(log.nu+1.96*sqrt(log.var))

    nu.est<-cbind(nu,se,lower,upper)
    colnames(nu.est)<-c("est","se","lower","upper")
    rownames(nu.est)<-temp.pao$unitnames[1:samp.units]
  }

  rnames=paste0(rep(temp.pao$unitnames[1:samp.units],temp.pao$nsurveys),"_",
            rep(temp.pao$surveynames,each=samp.units))
  pA.est<-get_real(real="<pA[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,rnames,fixed,temp.pao$nseasons+1)
  pB.est<-get_real(real="<pB[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,rnames,fixed,temp.pao$nseasons+1)
  rA.est<-get_real(real="<rA[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,rnames,fixed,temp.pao$nseasons+1)
  rBa.est<-get_real(real="<rBa[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,rnames,fixed,temp.pao$nseasons+1)
  if(param=="nu") {
    rho.est<-get_real(real="<rho[1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=samp.units,rnames,fixed,temp.pao$nseasons+1)
    ## derive rBA
    rBA.mat<-as.matrix(p_mat[(3*len+1):(4*len),]+p_mat[(4*len+1):(5*len),])
    logit.rBA<-rBA.mat%*%p.coeff[,1]
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
    log.rho<-rho.mat%*%p.coeff[,1]
    rho<-exp(log.rho)
    i=grep("^[bB]",rownames(VC)); if (multiseason) i=grep("^[dD]",rownames(VC))
    p.VC<-VC[i,i]
    log.var<-apply(rho.mat,1,function(xx) xx%*%p.VC%*%xx)
    se<-sqrt(log.var)*rho
    lower<-exp(log.rho-1.96*sqrt(log.var))
    upper<-exp(log.rho+1.96*sqrt(log.var))

    rho.est<-cbind(rho,se,lower,upper)
    colnames(rho.est)<-c("est","se","lower","upper")
    rownames(rho.est)<-rnames
}

  if (multiseason) {

    rnames=paste0(rep(temp.pao$unitnames[1:samp.units],temp.pao$nseasons-1),"_",
                  rep(1:(temp.pao$nseasons-1),each=samp.units))
    gamAB.est<-get_real(real="<gamAB",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,row.names=rnames,fixed=fixed)
    gamAb.est<-get_real(real="<gamAb",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    gamBAA.est<-get_real(real="<gamBAA",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    gamBAa.est<-get_real(real="<gamBAa",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    gamBaA.est<-get_real(real="<gamBaA",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    gamBaa.est<-get_real(real="<gamBaa",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)

    epsAB.est<-get_real(real="<epsAB",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    epsAb.est<-get_real(real="<epsAb",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    epsBAA.est<-get_real(real="<epsBAA",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    epsBAa.est<-get_real(real="<epsBAa",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    epsBaA.est<-get_real(real="<epsBaA",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
    epsBaa.est<-get_real(real="<epsBaa",v=v,con.offset=1,n1=data$nseasons-1,nunits=samp.units,rnames,fixed=fixed)
  }
  ##### check for warnings
  warn.conv<-check_conv_warn(v); warn.VC<-check_VC_warn(v)
   result<-list(modname=modname,
                model=list(psi=psi,p=p),dmat=list(psi=psi.dm,p=p.dm),
                data=temp.pao,outfile=outfile,
                neg2loglike=neg2loglike,npar=npar, aic=aic,
                beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,p=as.data.frame(p.coeff),p.VC=p.VC,VC=VC),
                real=list(psiA=as.data.frame(psiA.est),
                          psiBa=as.data.frame(psiBa.est), psiBA=as.data.frame(psiBA.est),
                          nu=as.data.frame(nu.est),
                          pA=as.data.frame(pA.est), pB=as.data.frame(pB.est),
                          rA=as.data.frame(rA.est), rBa=as.data.frame(rBa.est),
                          rBA=as.data.frame(rBA.est), rho=as.data.frame(rho.est)),
                warnings=list(conv=warn.conv,VC=warn.VC))
   if (multiseason) {
     result$model$gamma=gamma;                   result$model$epsilon=epsilon
     result$beta$gamma=as.data.frame(gam.coeff); result$beta$epsilon=as.data.frame(eps.coeff)
     result$beta$gamma.VC=gam.VC;                result$beta$epsilon.VC=eps.VC
     result$real$gammaAB=as.data.frame(gamAB.est);  result$real$epsilonAB=as.data.frame(epsAB.est)
     result$real$gammaAb=as.data.frame(gamAb.est);  result$real$epsilonAb=as.data.frame(epsAb.est)
     result$real$gammaBAA=as.data.frame(gamBAA.est); result$real$epsilonBAA=as.data.frame(epsBAA.est)
     result$real$gammaBAa=as.data.frame(gamBAa.est); result$real$epsilonBAa=as.data.frame(epsBAa.est)
     result$real$gammaBaA=as.data.frame(gamBaA.est); result$real$epsilonBaA=as.data.frame(epsBaA.est)
     result$real$gammaBaa=as.data.frame(gamBaa.est); result$real$epsilonBaa=as.data.frame(epsBaa.est)
     result$dmat$gamma=gamma.dm; result$dmat$epsilon=epsilon.dm
   }
   result.save<<-result
   if(param=="psiBA") class(result)<-c("occMod","so2sp1") else class(result)<-c("occMod","so2sp2")
   if(multiseason) class(result)<-c("occMod","do2sp1")
   return(result)
}
