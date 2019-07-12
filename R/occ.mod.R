
occMod<-function(model,cov.list=NULL,data,type=NULL,conf=0.95,
                  modname=NULL,paoname=NULL,outfile=NULL,
                  fixed=NULL,initvals=NULL, modfitboot=NULL,
                  VCoutopt="realbetavc",noDerived=F,randinit=0,maxfn='32000',quiet=1,limit.real=F,threads=1,chsq=FALSE,...){
  #'  Fit an occupancy model using Program PRESENCE
  #'
  #' \code{occMod} returns a list variable containing detection-history and covariate data.
  #'
  #' @param model list of formula for each paramter (eg., psi, p, gam, eps...)
  #' @param cov.list list of covariates used in the model
  #' @param data pao data object
  #' @param type a string indicating the type of occupancy model to fit. Currently implemented \code{type} options are:
  #'  \tabular{llll}{
  #'    \code{type} \tab Descriptionxxx \tab Function (see for details) \tab Reserved covariate names \cr
  #'    \code{"so"} \tab Static occupancy, single season \tab \code{\link{occMod_SO}}\tab \code{SURVEY} \cr
  #'    \code{"so.fp"} \tab Static occupancy, single season, false-positive detections\tab \code{\link{occMod_SO_fp}}\tab \code{SURVEY} \cr
  #'    \code{"so.mm"} \tab Static occupancy, single season, multi method (or multi-scale) \tab \code{\link{occMod_mm}}\tab \code{SURVEY}\cr
  #'    \code{"so.cd"} \tab Static occupancy, single season, correlated detections \tab \code{\link{occMod_SO_cd}}\tab \code{SURVEY}\cr
  #'    \code{"2sp"}\tab Static or Dynamic occupancy, 2 species, "psiBA" or "nu" param \tab \code{\link{occMod_2SP}}\tab \code{SURVEY,SP,INT,INT_o,INT_d,SEASON}\cr
  #'    \code{"so.2sp.fp"}\tab Static occupancy, single season, 2 species, false positive detections \tab \code{\link{occMod_SO_2SP_cd}}\tab \code{SURVEY,SP,INT,INT_o,INT_d}\cr
  #'    \code{"do.1"} \tab Dynamic occupancy, multi season, 1st param \tab \code{\link{occMod_DO}}\tab \code{SURVEY,SEASON}\cr
  #'    \code{"do.fp"} \tab Dynamic occupancy, multi season, with false-positive detections \tab \code{\link{occMod_DO_fp}}\tab \code{SURVEY,SEASON}\cr
  #'    \code{"do.4"} \tab Dynamic occupancy, multi season, 4th param \tab \code{\link{occMod_DO4}}\tab \code{SURVEY,SEASON}\cr
  #'    \code{"do.ms.1"}\tab Dynamic occupancy, multi season, multi state, 1st param\tab \code{\link{occMod_DO_ms1}}\tab \code{SURVEY,SEASON,DYN,PREV_STATE,STATE}\cr
  #'    \code{"do.ms.2"}\tab Dynamic occupancy, multi season, multi state, 2nd param\tab \code{\link{occMod_DO_ms2}}\tab \code{SURVEY,SEASON,DYN,PREV_STATE,STATE}\cr
  #'    \code{"do.2sp.1"}\tab Dynamic occupancy, multi season, 2 species, psiBA param \tab \code{\link{occMod_2SP}}\tab \code{SURVEY,SP,INT,INT_o,INT_d}\cr
  #'  }
  #' @param conf level for confidence interval (may be vector valued).
  #' @param modname (optional) a string with the user supplied name for a model.
  #'     If \code{NULL}, \code{modname} is created from the formulae supplied in \code{model}.
  #' @param paoname (optional) a string with the user supplied filename for PRESENCE data and
  #'    output files. If \code{NULL}, a generated name is used.
  #' @param outfile name for output file (use outfile='modname') for outfile named via model name
  #' @param fixed a data.frame with variables \code{param} and \code{value}. \code{param}
  #' should be a vector of string values containing the names of the real parameters (eg., "p(2)") to fix, and
  #' \code{value} be a numeic vector of the values the real parameter should be fixed at.
  #'  Ordering depends on model \code{type}.
  #' @param initvals a vector providing initial values for regression coefficients for
  #' the optimization procedure used by PRESENCE. Ordering depends on model \code{type}.
  #' @param modfitboot number of bootstraps for assessing model fit. Only works for \code{type = "so"}.
  #' @param VCoutopt option for PRESENCE computation/output of beta/real var-cov matrices
  #' @param noDerived if TRUE, doesn't print derived estimates from model
  #' @param maxfn max number of function calls before aborting optimization
  #' @param randinit number of random initial starting value vectors to try (defalut=0)
  #' @details choices for \code{VCoutopt} are:
  #' \itemize{
  #'  \item{"nose"}{only compute likelihood and beta estimates, no beta var-cov matrix or beta SE's (useful for faster model selection of big models or simulations)}
  #'  \item{"betavc"}{compute but don't print beta var-cov matrix so we get beta SE's, don't compute real params}
  #'  \item{"noreal"}{compute and print beta beta var-cov matrix, don't compute real params (I think this is what you want for RPresence)}
  #'  \item{"novcs"}{compute real params and real var-cov matrix, don't print either var-cov matrix.}
  #'  \item{"realbetavc"}{(default) compute real params and real var-cov matrix, print only beta var-cov matrix.}
  #'  \item{"bothvc"}{print both var-cov matrices.}
  #'}
  #' @param quiet 0=normal output from Presence,  1=no output (for simulations or repeated runs)
  ##@param threads number of threads to use for multi-core computers
  #' @param ... additional arguments to occMod.XX function.

  #'@export
  #' @return Returns a list of class "occMod" and class \code{type}.
  #'   The \code{occMod} object has the following objects:
  #'    \item{modname}{the name of the model}
  #'    \item{model}{a named list containing the right-hand formula for each real parameter type.}
  #'    \item{data}{a copy of the temporary \code{pao} object that is created from the \code{model} argument then used to create the temporary PRESENCE data file for fitting the model.}
  #'    \item{outfile}{name of PRESENCE output file.}
  #'    \item{neg2loglike}{-2*log(likelihood) value for the model evaluated at the maximum likelihood estimates.}
  #'    \item{npar}{number of parameters in the model, calculated as the number of regression coefficients.}
  #'    \item{aic}{value of Akaike's Information Criterion.}
  #'    \item{beta}{list of estimated regression coefficents (sometimes called beta parameters), standard errors and variance-covariance matrix for each parameter type. Exact content depends on the \code{type} of model being fit to the data. Values will be on either the logit or log scale depending on whether the parameter type is a probability or not, respectively.}
  #'    \item{real}{list of real parameter estimates on their natural scale (i.e., 0-1 for probabilities). A data frame for each parameter type the contains the estimate, standard error and limits of the requested confidence intervals. Exact content depends on the \code{type} of model being fit to the data.}
  #'    \item{derived}{list of parameter estimates that are derived from real parameters and not directly estimated. A data frame for each parameter type the contains the estimate, standard error and limits of the requested confidence intervals. Exact content depends on the \code{type} of model being fit to the data.}
  #'    \item{modfit}{list of statistics from parametric bootstrap-based model fit.}
  #'    \item{warnings}{a list containing any warnings from the PRESENCE output file; \code{conv}=number of significant digits convergence was achieved to (=\code{NULL} if >7); and \code{VC}=1 if there is a problem with the variance-covariance matrix (likely non-invertiable observed Hessian matrix), =\code{NULL} otherwise. If \code{conv}<3 then may want to attempt alternative initial values or rescaling covariates to see if different results are obtained. If \code{VC}=1 then all standard errors and confidence intervals should be ignored.}
  #'    \item{version}{a list containing the version numbers for PRESENCE and RPresence used when fitting the model}
  #'
  #' @author Darryl MacKenzie
  #'
  #' @seealso \code{\link{occMod_SO}},\code{\link{occMod_mm}},\code{\link{occMod_2SP}},
  #'          \code{\link{occMod_DO}},\code{\link{occMod_DO4}},\code{\link{occMod_DO_ms2}},\code{\link{occMod_2SPfp}}
  #'
  #' @examples
  #'# load a PRESENCE data file
  #'filename<-system.file("extdata/weta.pao",package="RPresence")
  #'weta.data<-readPao(filename)
  #'
  #'## convert indicator variables to categorical covariates
  #'weta.data$unitcov$Habitat[weta.data$unitcov$Browsed==1] <- "browsed"
  #'weta.data$unitcov$Habitat[weta.data$unitcov$Unbrowsed==1] <- "unbrowsed"
  #'weta.data$survcov$Obs<-as.factor(
  #'      1*weta.data$survcov$Obs1+
  #'      2*weta.data$survcov$Obs2+
  #'      3*weta.data$survcov$Obs3)
  #'
  #'## fit some models
  #'mod1<-occMod(model=list(psi~Habitat,p~SURVEY+Obs),data=weta.data,type="so")
  #'mod2<-occMod(model=list(psi~Habitat,p~Obs),data=weta.data,type="so")
  #'mod3<-occMod(model=list(psi~Habitat,p~SURVEY),data=weta.data,type="so")
  #'mod4<-occMod(model=list(psi~Habitat,p~1),data=weta.data,type="so")
  #'mod5<-occMod(model=list(psi~1,p~SURVEY+Obs),data=weta.data,type="so")
  #'mod6<-occMod(model=list(psi~1,p~Obs),data=weta.data,type="so")
  #'mod7<-occMod(model=list(psi~1,p~SURVEY),data=weta.data,type="so")
  #'mod8<-occMod(model=list(psi~1,p~1),data=weta.data,type="so")
  #'#
  #'## create AIC table
  #'models<-list(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
  #'results<-createAicTable(models)
  #'summary(results)
  #'
  #'## perform model averageing
  #'psi.ma<-modAvg(aic.tab=results,param="psi")
  #'p.ma<-modAvg(aic.tab=results,param="p")
  #'
  #' \dontrun{
  #' m1<-occMod(model=list(psi~habitat,gamma~habitat+SEASON,epsilon~habitat+SEASON,
  #'                        p~habitat+SEASON),data=data,type="do.1")
  #' }

  ## start with some checking
  #check data is pao object
  if(!("pao"%in%class(data))){ stop("data must be a pao object!!")  }
  # check implemented model type
  valid.mod.types<-data.frame(
    id=c("so","so.fp","do.1","do.fp","do.4","so.mm","so.cd","so.2sp.1","so.2sp.cd","so.2sp.2",
         "do.ms.1","do.ms.2","do.2sp.1","so.2sp.fp","do.mm"),
    descr=c("Static occupancy, single season",
            "Static occupancy, single season false-positive",
            "Dynamic occupancy, multi season, 1st param",
            "Dynamic occupancy, multi season, false positive",
            "Dynamic occupancy, multi season, 4th param",
            "Static occupancy, single season, multi method",
            "Static occupancy, single season, correlated detectinos",
            "Static occupancy, single season, 2 species, psiBA param",
            "Static occupancy, single season, 2 species, psiBA param, corr.detections",
            "Static occupancy, single season, 2 species, nu param",
            "Dynamic occupancy, multi season, multi state, 1st param",
            "Dynamic occupancy, multi season, multi state, 2nd param",
            "Dynamic occupancy, multi season, 2 species, psiBA param",
            "Static occupancy, single season, 2 species, psiBA param, false pos detections",
            "Dynamic occupancy, multi season, multi method"
            )
  )
  if (sum(type==valid.mod.types)==0) {
    cat("\nInvalid model type specified. Currently implemented model types are:\n")
    print(valid.mod.types)
    return(NULL)
  }
  #create model name if none supplied
  if(is.null(modname)) {
    modname<-paste0(lapply(model,function(xx){paste0(deparse(xx[[2]]),"(",Reduce(paste,deparse(xx[[3]])),')')}),collapse='')
    modname<-gsub('(1)','()',gsub(" ","",modname),fixed=T)
    modname<-gsub("+"," P ",modname,fixed=TRUE)
    modname<-gsub("*"," X ",modname,fixed=TRUE)
    modname<-gsub(":"," T ",modname,fixed=TRUE)
  }

  if (is.null(outfile)) outfile=paste0('del_',Sys.getpid(),'_',sample(1:99999,1),'.out') else
    if (outfile=="modname") outfile=paste0(modname,".out")
  if (length(grep('_modname',outfile))>0) outfile=gsub('_modname',paste0('_',modname,'.out'),outfile)

  miscopts=c(VCoutopt,paste0('maxfn=',maxfn),paste0('name=',modname))
  if (quiet) miscopts=c(miscopts,'quiet')
  if (limit.real) miscopts=c(miscopts,'lmt=1')
  if (noDerived) miscopts=c(miscopts,'noderived')
  if (!is.null(modfitboot)) miscopts=c(miscopts,paste0("boot2=",modfitboot))
  if (threads>1) miscopts=c(miscopts,paste0('nthreads=',threads))
  if (randinit>0) miscopts=c(miscopts,paste0('nrand=',randinit))
  if (chsq) miscopts=c(miscopts,'chsq')
  miscopts=c(miscopts,'')

  bothcovs=data$survcov; bnames=colnames(bothcovs)
  if(!is.null(data$unitcov))
    if (nrow(data$unitcov)>0)
      for (i in 1:ncol(data$unitcov)) {
        bothcovs=cbind(bothcovs,rep(data$unitcov[,i],data$nsurveys))
        bnames=c(bnames,colnames(data$unitcov)[i])
      }
  colnames(bothcovs)=bnames
  e=new.env(); res=NA
  #################################################################
  # run specified model
  #################################################################
  if(type=="so"){ ## static occupancy, single season
    params<-c("psi","p"); extract_model(model,params,e,type)
    #check psi and p both defined
    if(is.null(e$psi) | is.null(e$p)){
      stop("psi and p must both be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov<-data$survcov
    }
    res<-occMod_SO(psi=e$psi,psi.cov=e$psi.cov,         p=e$p,p.cov=e$p.cov,
                     modname=modname,paoname=paoname,outfile=outfile,
                     model=100,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################
  if(type=="so.fp") { # static occupancy, false positive
    params<-c("psi","p11","p10","b"); extract_model(model,params,e,type)

    #check psi, p, p11, p10 all defined
    if(is.null(e$psi) | is.null(e$b) | is.null(e$p11) | is.null(e$p10)){
      stop("psi, b, p11 and p10 must all be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$b.cov)) e$b.cov=cov.list$b.cov else
	    if(!is.null(data$unitcov)) e$b.cov=bothcovs else e$b.cov<-data$survcov
    if(!is.null(cov.list$p11.cov)) e$p11.cov=cov.list$p11.cov else
                if(!is.null(data$unitcov)) e$p11.cov=bothcovs else e$p11.cov=data$survcov
    if(!is.null(cov.list$p10.cov)) e$p10.cov=cov.list$p10.cov else
                if(!is.null(data$unitcov)) e$p10.cov=bothcovs else e$p10.cov<-data$survcov

    res<-occMod_SO_fp(psi=e$psi,psi.cov=e$psi.cov,
                       p11=e$p11,p11.cov=e$p11.cov,
                       p10=e$p10,p10.cov=e$p10.cov,
                       b=e$b,b.cov=e$b.cov,
                       modname=modname,paoname=paoname,outfile=outfile,
                       model=103,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################
  if(type=="so.mm" | type=="do.mm") { # static or dynamic occupancy, multi method
    params<-c("psi","theta","p","gamma","epsilon"); extract_model(model,params,e,type)

    #check psi, theta and p all defined
    if(is.null(e$psi) | is.null(e$theta) | is.null(e$p)){
      stop("psi, theta and p must all be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$theta.cov)) e$theta.cov=cov.list$theta.cov else {
      if(!is.null(data$unitcov)) e$theta.cov=bothcovs else e$theta.cov=data$survcov
    }
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    res<-occMod_mm(psi=e$psi,psi.cov=e$psi.cov,
                    theta=e$theta,theta.cov=e$theta.cov,
                    p=e$p,p.cov=e$p.cov,
                    gamma=e$gamma, gamma.cov=e$gamma.cov,
                    epsilon=e$epsilon, epsilon.cov=e$epsilon.cov,
                    modname=modname,paoname=paoname,outfile=outfile,
                    fixed=fixed,initvals=initvals,data=data,miscopts)
  }
  #################################################################
  if(type=="so.cd") { # static occupancy, correlated detections
    params<-c("psi","theta","p","th0pi"); extract_model(model,params,e,type)

    #check psi, theta and p all defined
    if(is.null(e$psi) | is.null(e$theta) | is.null(e$p)){
      stop("psi, theta and p must all be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$theta.cov)) e$theta.cov=cov.list$theta.cov else e$theta.cov=bothcovs
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    if(!is.null(cov.list$th0pi.cov)) e$th0pi.cov<-cov.list$th0pi.cov else e$th0pi.cov<-data$unitcov
    res<-occMod_SO_cd(psi=e$psi,psi.cov=e$psi.cov, theta=e$theta, theta.cov=e$theta.cov,
                       p=e$p,    p.cov=e$p.cov,     th0pi=e$th0pi, th0pi.cov=e$th0pi.cov,
                       modname=modname,paoname=paoname,outfile=outfile,
                       model=103,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################
  if(type=="so.2sp.1" | type=="so.2sp.2" | type=="do.2sp.1"){ # static or dynamic occupancy, 2 species, psiBA param
    param=ifelse(type != "so.2sp.2","psiBA","nu")
    params<-c("psi","p","gamma","epsilon"); extract_model(model,params,e,type)
    #check psi and p all defined
    if(is.null(e$psi) | is.null(e$p)) stop("psi and p must both be defined in the model list!")
    if(!is.null(cov.list$psi.cov)) e$psi.cov=cov.list$psi.cov else e$psi.cov=data$unitcov
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    if(!is.null(cov.list$gamma.cov)) e$gamma.cov=cov.list$gamma.cov else {
      gcov=NULL; for (i in 2:data$nseasons) gcov=rbind(gcov,data$unitcov); #  replicate site covariates for
      e$gamma.cov=gcov                                                     #   each season (except last)
    }
    if(!is.null(cov.list$epsilon.cov)) e$epsilon.cov=cov.list$epsilon.cov else {
      ecov=NULL; for (i in 2:data$nseasons) ecov=rbind(ecov,data$unitcov)   #  replicate site covariates for
      e$epsilon.cov=ecov                                                    #   each season (except last)
    }
    # if(sum(cov.list$psi.cov=="-")>0 || sum(cov.list(p.ocv)=="-")>0 ||) {
    #    cat("\n******* covariate cannot contain missing values ***\n\n")
    #  }

    res<-occMod_2SP(psi=e$psi,psi.cov=e$psi.cov,
                        p=e$p,p.cov=e$p.cov,
                        gamma=e$gamma, gamma.cov=e$gamma.cov,
                        epsilon=e$epsilon, epsilon.cov=e$epsilon.cov,
                        param=param,sp.contr=TRUE,
                        modname=modname,paoname=paoname,outfile=outfile,
                        model=3000,fixed=fixed,initvals=initvals,data=data,miscopts)
  }

  #################################################################
  if(type=="so.2sp.cd"){ # static or dynamic occupancy, 2 species, psiBA param, corr.detections
    param="psiBA"
    params<-c("psi","theta","p","th0pi"); extract_model(model,params,e,type)

    #check psi and p all defined
    if(is.null(e$psi) | is.null(e$p)) stop("psi and p must both be defined in the model list!")
    if(!is.null(cov.list$psi.cov)) e$psi.cov=cov.list$psi.cov else e$psi.cov=data$unitcov
    if(!is.null(cov.list$theta.cov)) e$theta.cov=cov.list$theta.cov else {
      if(!is.null(data$unitcov)) e$theta.cov=bothcovs else e$theta.cov=data$survcov
    }
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    if(!is.null(cov.list$th0pi.cov)) e$th0pi.cov=cov.list$th0pi.cov else e$th0pi.cov=data$unitcov

    res<-occMod_SO_2SP_cd(psi=e$psi,psi.cov=e$psi.cov, theta=e$theta, theta.cov=e$theta.cov,
                    p=e$p,p.cov=e$p.cov, th0pi=e$th0pi, th0pi.cov=e$th0pi.cov,
                    param="psiBA",sp.contr=TRUE,
                    modname=modname,paoname=paoname,outfile=outfile,
                    model=3000,fixed=fixed,initvals=initvals,data=data,miscopts)
  }
  #################################################################
  if(type=="so.2sp.fp"){ # static occupancy, 2 species, psiBA param
    params<-c("psi","p","omeg","conf"); extract_model(model,params,e,type)
    #check psi and p all defined
    if(is.null(e$psi) | is.null(e$p)) stop("psi and p must both be defined in the model list!")
    if(!is.null(cov.list$psi.cov)) e$psi.cov=cov.list$psi.cov else e$psi.cov=data$unitcov
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    if(!is.null(cov.list$omeg.cov)) e$omeg.cov=cov.list$omeg.cov else {
      if(!is.null(data$unitcov)) e$omeg.cov=bothcovs else e$omeg.cov=data$survcov
    }
    if(!is.null(cov.list$conf.cov)) e$conf.cov=cov.list$conf.cov else {
      if(!is.null(data$unitcov)) e$conf.cov=bothcovs else e$conf.cov=data$survcov
    }
    res<-occMod_2SPfp(psi=e$psi,psi.cov=e$psi.cov,p=e$p,p.cov=e$p.cov,
                    omeg=e$omeg, omeg.cov=e$omeg.cov,
                    conf=e$conf, conf.cov=e$conf.cov,
                    param=param,sp.contr=TRUE,
                    modname=modname,paoname=paoname,outfile=outfile,
                    model=3000,fixed=fixed,initvals=initvals,data=data,miscopts)
  }

  #################################################################

  if(type=="do.1"){ # dynamic occupancy, multi season initial param
    params<-c("psi","gamma","epsilon","p","theta","th0pi"); extract_model(model,params,e,type)

    #check psi,gamma,epsilon and p all defined
    if(is.null(e$psi) | is.null(e$gamma) | is.null(e$epsilon) | is.null(e$p)){
      stop("psi, gamma, epsilon and p must both be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov=cov.list$psi.cov else e$psi.cov=data$unitcov
    if(!is.null(cov.list$gamma.cov)) e$gamma.cov=cov.list$gamma.cov else e$gamma.cov=bothcovs
    if(!is.null(cov.list$epsilon.cov)) e$epsilon.cov=cov.list$epsilon.cov else e$epsilon.cov=bothcovs
    if(!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if(!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    res<-occMod_DO(psi=e$psi,psi.cov=e$psi.cov,
                     gamma=e$gamma,gamma.cov=e$gamma.cov,
                     epsilon=e$epsilon,epsilon.cov=e$epsilon.cov,
                     p=e$p,p.cov=e$p.cov,
                     theta=e$theta, theta.cov=e$theta.cov, th0pi=e$th0pi, th0pi.cov=e$th0pi.cov,
                     modname=modname,paoname=paoname,outfile=outfile,
                     model=200,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################
  if(type=="do.fp") { # dynamic occupancy, false positive
    params<-c("psi","epsilon","gamma","p11","p10","b"); extract_model(model,params,e,type)

    #check psi, p, p11, p10 all defined
    if(is.null(e$psi) | is.null(e$b) | is.null(e$p11) | is.null(e$p10) | is.null(e$gamma) | is.null(e$epsilon)){
      stop("psi, gamma, epsilon, b, p11 and p10 must all be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$gamma.cov)) e$gamma.cov=cov.list$gamma.cov else {
      e$gamma.cov=bothcovs[1:(data$nunits*(data$nseasons-1)),]
    }
    if(!is.null(cov.list$epsilon.cov)) e$epsilon.cov=cov.list$epsilon.cov else {
      e$epsilon.cov=bothcovs[1:(data$nunits*(data$nseasons-1)),]
    }
    if(!is.null(cov.list$b.cov)) e$b.cov=cov.list$b.cov else
      if(!is.null(data$unitcov)) e$b.cov=bothcovs else e$b.cov<-data$survcov
      if(!is.null(cov.list$p11.cov)) e$p11.cov=cov.list$p11.cov else
        if(!is.null(data$unitcov)) e$p11.cov=bothcovs else e$p11.cov=data$survcov
        if(!is.null(cov.list$p10.cov)) e$p10.cov=cov.list$p10.cov else
          if(!is.null(data$unitcov)) e$p10.cov=bothcovs else e$p10.cov<-data$survcov

          res<-occMod_DO_fp(psi=e$psi,psi.cov=e$psi.cov,
                            gamma=e$gamma,gamma.cov=e$gamma.cov,epsilon=e$epsilon,epsilon.cov=e$epsilon.cov,
                            p11=e$p11,p11.cov=e$p11.cov,
                            p10=e$p10,p10.cov=e$p10.cov,
                            b=e$b,b.cov=e$b.cov,
                            modname=modname,paoname=paoname,outfile=outfile,
                            model=203,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################
  if(type=="do.4"){ # dynamic occupancy, multi season 4th param
    params<-c("psi","p"); extract_model(model,params,e,type)

    #check psi and p all defined
    if(is.null(e$psi) | is.null(e$p))  stop("psi and p must both be defined in the model list!")
    if (!is.null(cov.list$psi.cov)) e$psi.cov=cov.list$psi.cov else {
      if (!is.null(data$survcov)) e$psi.cov=bothcovs else e$psi.cov=data$unitcov
    }
    if (!is.null(cov.list$p.cov)) e$p.cov=cov.list$p.cov else {
      if (!is.null(data$unitcov)) e$p.cov=bothcovs else e$p.cov=data$survcov
    }
    res<-occMod_DO4(psi=e$psi,psi.cov=e$psi.cov,p=e$p,p.cov=e$p.cov,
                     modname=modname,paoname=paoname,outfile=outfile,
                     model=240,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################

  if(type=="do.ms.1"){ # dynamic occupancy, multi season multi state 1st param
    params<-c("psi","phi","p"); extract_model(model,params,e,type)

    #check psi, p defined
    if(is.null(e$psi) | is.null(e$p) | is.null(e$phi)) stop("\npsi, phi and p must be defined in the model list!\n\n")

    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$phi.cov)) e$phi.cov<-cov.list$phi.cov else e$phi.cov<-data$unitcov
    if(!is.null(cov.list$p.cov)){ e$p.cov<-cov.list$p.cov } else {
      if(!is.null(data$unitcov)) e$p.cov<-bothcovs else e$p.cov<-data$survcov
    }

    res<-occMod_DO_ms1(psi=e$psi,psi.cov=e$psi.cov, phi=e$phi, phi.cov=e$phi.cov,
                     p=e$p,p.cov=e$p.cov,
                     modname=modname,paoname=paoname,outfile=outfile,
                     model=6000,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts=miscopts)
  }
  #################################################################

  if(type=="do.ms.2"){ # dynamic occupancy, multi season multi state 2nd param
    params<-c("psi","r","p","delta"); extract_model(model,params,e,type)

    #check psi,r, p and delta all defined
    if(is.null(e$psi) | is.null(e$r) | is.null(e$p) | is.null(e$delta)){
      stop("psi and p must both be defined in the model list!")
    }
    if(!is.null(cov.list$psi.cov)) e$psi.cov<-cov.list$psi.cov else e$psi.cov<-data$unitcov
    if(!is.null(cov.list$r.cov)) e$r.cov<-cov.list$r.cov else e$r.cov<-data$unitcov
    if(!is.null(cov.list$p.cov)){
      e$p.cov<-cov.list$p.cov
    } else {
      if(!is.null(data$unitcov)) e$p.cov<-bothcovs else e$p.cov<-data$survcov
    }
    if(!is.null(cov.list$delta.cov)) e$delta.cov<-cov.list$delta.cov else e$delta.cov<-data$survcov

    res<-occMod_DO_ms2(psi=e$psi,psi.cov=e$psi.cov,
                       r=e$r,r.cov=e$r.cov,
                       p=e$p,p.cov=e$p.cov,
                       delta=e$delta,delta.cov=e$delta.cov,
                       modname=modname,paoname=paoname,outfile=outfile,
                       model=6100,fixed=fixed,initvals=initvals,data=data,conf=conf,miscopts)
  }
  #################################################################
  if (substr(outfile,1,4)=="del_") file.remove(outfile)
  rm(e); return(res)  #tidy up environment
}
