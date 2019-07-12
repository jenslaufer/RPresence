
#   R-PRES.r : functions to use PRESENCE from R  (Ver 3 - Nov 2013)
### DM has modifed code originally supplied by JH

readPao <- function(f) {
  #'  Read a pao data object.
  #'
  #' \code{readPao} Creates a pao data object from a PRESENCE data file.
  #'
  #'@param f A string specifying the name of the PRESENCE data file.
  #'@export
  #' @import stats
  #' @importFrom graphics plot polygon points lines abline
  #' @importFrom grDevices grey
  #' @importFrom utils combn
  #'@return list of class "\code{pao}" containing the data required to fit an
  #'  occupancy model using \code{\link{occMod}}.
  #'  The list contains the following objects:
  #'
  #'   \item{nunits}{number of units or sampling units, i.e., \code{nrow(data)}.}
  #'   \item{nsurveys}{total number of surveys, i.e., \code{ncol(data)}.}
  #'   \item{nseasons}{number of seasons of data.}
  #'   \item{det.data}{detection/nondetection data, i.e., \code{data}.}
  #'   \item{unitcov}{see above. If no covariates are supplied a vector of 1's named TEMP is added.}
  #'   \item{survcov}{see above. A factor called SURVEY is automatically added to allow survey-specific detection probabilities (but equal across units).}
  #'   \item{nsurveyseason}{see above.}
  #'   \item{title}{see above.}
  #'   \item{unitnames}{see above.}
  #'   \item{surveynames}{a vector of strings providing labels for the survey occasions. Of little practical consequence.}
  #'   \item{paoname}{see above. Equals \code{f} for \code{createPao}.}
  #'
  #'@author Darryl MacKenzie and Jim Hines
  #'
  #'@seealso \code{\link{occMod}}
  #'
  #'@examples
  #' readPao(system.file('extdata/skinks.pao',package='RPresence'))
  #'

  str2mat <- function(s,nr,nc,sepchr) {  # convert vector of strings to matrix
    v=matrix(suppressWarnings(as.numeric(unlist(strsplit(s,sepchr)))),nr,nc,byrow=T)
    return(as.data.frame(v))
  }
  a=scan(f,sep='\n',what='c',quiet=TRUE);                        #  pao$nunits
  b=as.numeric(unlist(strsplit(a[1],'\t')))                      #  pao$nsurveys
  nunits=b[1]; nsurveys=b[2];      nc=length(unlist(strsplit(a[2],'\t')))
  d=str2mat(a[1+1:nunits],nunits,nc,'\t')                  #  pao$det.data   - detection history matrix
  unitcov=covname1=NULL;
  nunitcov=as.numeric(a[nunits+2]);  fp=nunits+2                 #  pao$nunitcov - num unit covariates
  if (nunitcov>0) {
    for (i in 1:nunitcov) {
      covname1=c(covname1,a[fp+1])  ;                        #  pao$covname1 - vector of unit cov names
      unitcov=cbind(unitcov,a[fp+1+1:nunits]);               #  pao$unitcov     - matrix of unit covs
      fp=fp+nunits+1
    }
    colnames(unitcov)=covname1;
    unitcov=as.data.frame(unitcov,stringsAsFactors=F)
    for (i in 1:nunitcov)
      if (length(grep("[.][0-9]",unitcov[,i]))>1) unitcov[,i]=as.numeric(unitcov[,i])
  }
  covname2="SURVEY"
  survcov<-data.frame(SURVEY=as.factor(rep(1:nsurveys,each=nunits)))
  nsurvcov=as.numeric(a[fp+1])+1                             #  pao$nsurvcov - num survey covs
  fp=fp+1;
  if (nsurvcov>1) {
    for (i in 2:nsurvcov) {
      covname2=c(covname2,a[fp+1])                          #  pao$covname2 -
      aa=str2mat(a[fp+1+1:nunits],nunits,nsurveys,'\t')
      survcov<-data.frame(cbind(survcov,as.numeric(c(as.matrix(aa)))))
      fp=fp+nunits+1
    }
  }
  names(survcov)=covname2;
  s=a[fp+1]; nmethods=1; s1=gsub(';.+','',s);
  if (s1!=s) nmethods=as.numeric(gsub('.+;','',s))
  nsurveyseason=as.numeric(unlist(strsplit(s1,','))); fp=fp+1
  if (length(nsurveyseason)==1) {
    nseasons=max(1,nsurveys/nsurveyseason);
    nsurveyseason=rep(nsurveyseason,nseasons)
  } else nseasons=length(nsurveyseason);
  title=a[fp+1]; fp=fp+1; unitnames=paste0('unit',1:nunits)
  if ((fp+nunits)<=length(a)) {
    unitnames=a[fp+1:nunits]; fp=fp+nunits;  #  pao$unitnames
  }
  surveynames=NULL;
  for (i in 1:length(nsurveyseason)) surveynames=c(surveynames,paste0(i,"-",1:nsurveyseason[i]))
  if ((fp+nsurveys)<=length(a)) {
    surveynames=a[fp+1:nsurveys]; fp=fp+1;   #  pao$surveynames
  }
  rownames(d)=unitnames
  if (ncol(d)>nsurveys) { frq=d[,ncol(d)]; d=d[,-ncol(d)]} else frq=rep(1,nrow(d))

  result<-list(nunits=nunits,nsurveys=nsurveys,nseasons=nseasons,det.data=d,
               nunitcov=nunitcov,unitcov=unitcov,nsurvcov=nsurvcov,survcov=survcov,
               nsurveyseason=nsurveyseason,nmethods=nmethods,title=title,unitnames=unitnames,
               surveynames=surveynames,paoname=f,frq=frq)
  class(result)<-c("pao")
  return(result)
}

writePao <- function(pao) {        #  writes pao file
  #'  Write a pao data object to disk file.
  #'
  #' \code{writePao} Writes a pao data object to disk.
  #'
  #'@param pao A pao object created with \code{readPao} or \code{createPao}.
  #'@export
  #'@return name of file written.
  #'
  #'@author Darryl MacKenzie and Jim Hines
  #'
  #'@seealso \code{\link{occMod}}
  #'
  #'@examples
  #'\dontrun{
  #' writePao(paodata)
  #'}
  #'
  #### edited for long-storage of p covariates and R missing values
  f=pao$paoname; if (is.null(f)) f="tmp.pao"

  nsurveys=ifelse(is.null(pao$nsurveys),dim(pao$det.data)[2],pao$nsurveys)
  nunits=dim(pao$det.data)[1];
  cat(sprintf('%d\t%d',nunits,nsurveys),file=f,sep='\n')    #  write no. units, surveys
  dd=pao$det.data; if (!is.null(pao$frq)) dd=cbind(pao$det.data,pao$frq)
  write.table(dd,file=f,na='-',sep='\t',quote=F,row.names=F,col.names=F,append=TRUE)   #  write detection data matrix
  nunitcov=nsurvcov=0;

  if (is.vector(pao$unitcov)) pao$unitcov=matrix(pao$unitcov,ncol=1)
  if (is.matrix(pao$unitcov) | is.data.frame(pao$unitcov)) {
    nunitcov=ncol(pao$unitcov)
    if (nrow(pao$unitcov)<2) nunitcov=0
    cat(nunitcov,sep='\n',file=f,append=TRUE)                                       #  write no. unit covariates
    cnames=colnames(pao$unitcov)
    if (!is.null(pao$covname1))cnames=pao$covname1
    if (nunitcov>0) {
      for (i in 1:nunitcov) {
        cat(cnames[i],sep='\n',file=f,append=TRUE)                                  #  write unit cov. name
        write.table(pao$unitcov[,i],na="-",file=f,quote=F,row.names=F,col.names=F,append=TRUE)      #  write sit cov. vector
      }
    }
  } else cat('0',sep='\n',file=f,append=T)

  if (is.vector(pao$survcov)) pao$survcov=matrix(pao$survcov,ncol=1)
  if (is.matrix(pao$survcov)| is.data.frame(pao$survcov)) {
    nsurvcov=ncol(pao$survcov)
    if (nrow(pao$survcov)<2) nsurvcov=0
    cat(nsurvcov,sep='\n',file=f,append=TRUE)                                       #  write no. survey covariates
    if (nsurvcov>0) {
      for (i in 1:nsurvcov) {
        cat(colnames(pao$survcov)[i],sep='\n',file=f,append=TRUE)                   #  write survey cov. name
        write.table(array(pao$survcov[,i],dim=dim(pao$det.data)),na="-",file=f,sep='\t',quote=F,row.names=F,col.names=F,append=TRUE) #  write survey cov. matrix
      }
    }#   nunit rows, one on top of another
  } else cat('0',sep='\n',file=f,append=T)
    #  write no. surveys per season
  s1=paste(nsurveys); s2="1"
  if (length(grep("nsurveyseason",names(pao)))>0) s1=paste0(pao$nsurveyseason,collapse=',');
  if (length(grep("nmethods",names(pao)))>0) s2=pao$nmethods
  s=paste(s1,s2,sep=';')
  cat(s,'\n',file=f,append=TRUE)
  s=paste('no title file=',f)
     #  write title
  if (length(grep("title",names(pao)))>0) s=pao$title
  cat(s,'\n',file=f,append=TRUE)
   #  write unitnames
  if (is.null(pao$unitnames)) unitnames=dimnames(pao$det.data)[1] else unitnames=pao$unitnames
  write.table(unitnames,file=f,quote=F,row.names=F,col.names=F,append=TRUE)
  return(f)
}

runPresence <- function(a,dm,model,modname='m1',fixed=NULL,initvals=NULL,outfile,miscopts='') {
  #'  call PRESENCE to produce model estimates (replaces "write_dm_and_run")
  #'
  #' \code{runPresence} calls PRESENCE to produce model estimates.
  #'
  #' Input is passed
  #' to PRESENCE via arguments (no pao or design matrix input files needed).  There are
  #' three choices for output:
  #'
  #' \itemize{
  #' \item outfile=NULL - no output file saved, results in R variables only,
  #' \item outfile="modname" - output file saved with name of model and ".out" appended,
  #' \item outfile="xyz" - output file saved with user-specified name, (eg., "xyz.out")
  #' }
  #' The first option is best for running simulations of a model using "Snowfall"
  #' for multi-tasking.  The 2nd option is useful when you would like to find
  #' the full text output of a model from Windows explorer.  The 3rd option allows
  #' you to specify any output file name.
  #'
  #'@param a pao object
  #'@param dm list of design matrices (should be 6 in list)
  #'@param model model number (not used anymore)
  #'@param modname name of model (eg., "psi(.)p(.)")
#  #'@param fixed vector of fixed values (eg., fixed(1)=0, fixed(4)=.1 -> fixed=c(0,.1); names(fixed)=c(1,4))
  #'@param fixed data.frame of fixed values (e.g., fixed$params<-c("p[1]","p[3]"),fixed$value<-0)
  #'@param initvals vector of beta initial values
  #'@param miscopts options (VCoutput, quiet, chsq,...)
  #'@param outfile name for output file (use outfile='modname') for outfile named via model name
  #'@export
  #'@return list(LL=log-likelihood, beta=beta estimates, betaVC=var-cov matrix of betas).
  #'
  #'@author Darryl MacKenzie and Jim Hines
  #'
  #'@seealso \code{\link{occMod}}
  #'
  #'@examples
  #' runPresence("test.pao",list(dm1,dm2,dm3,dm4,dm5,dm6),100,"psi(.)p(.)")
  #'
  fnname='_Z9rpresencePiS_S_S_S_PdS_S0_PPcS_S0_S2_S_S2_S2_S2_S_S2_S2_S0_S0_S_S0_S2_S2_S0_'

  dmdims=unlist(lapply(list(dm[[1]],dm[[2]],dm[[3]],dm[[4]],dm[[5]],dm[[6]]),
                       function(x)  if (length(dim(x))==2) dim(x) else rep(0,2)))
  for (i in 2:6) if (is.null(dm[[i]])) dm[[i]]=matrix(0,1,1)  #  cant pass null pointer to C function
  dms=c(t(dm[[1]]),t(dm[[2]]),t(dm[[3]]),t(dm[[4]]),t(dm[[5]]),t(dm[[6]]));
  nreal=sum(dmdims[c(1,3,5,7,9,11)]); nbeta=sum(dmdims[c(2,4,6,8,10,12)])
  fv=rep(-999,nreal); if (!is.null(fixed)) fv[fixed$idx]=fixed$value
  if (is.null(initvals)) initvals=rep(0,nbeta)
  unitcovnames=survcovnames='asdf'
  if (is.null(a$unitcov)) a$unitcov=0 else unitcovnames=colnames(a$unitcov)
  if (is.null(a$survcov)) a$survcov=0 else survcovnames=colnames(a$survcov)
  neighbors=matrix(0,a$nunits,a$nunits); nbr_wgts=rep(0,a$nunits)

  rvals=rep(0,1+nbeta+nbeta^2);
  if (is.null(a$frq)) a$frq=rep(1,a$nunits)
  a$det.data[is.na(a$det.data)]=-1; a$unitcov[is.na(a$unitcov)]=-9999; a$survcov[is.na(a$survcov)]=-9999
  dd=unlist(a$det.data[,1:a$nsurveys])
  if (a$nsurvcov>0)
    for (j in 1:ncol(a$survcov)) {
      if (sum(( dd!= -1) & (a$survcov[,j]==-9999))>0)
        warning('\nMissing value for survey covariate (',survcovnames[j],') for non-missing detection data\n')
    }
  if (outfile=='savepao.out') {
    writePao(a); fnm=gsub('.pao','.dm',a$paoname); file.remove(fnm)
    for (i in 1:6) {
      if (nrow(dm[[i]])==1 && ncol(dm[[i]])==1 && dm[[i]][1,1]==0) cat(i-1,0,0,'\n',file=fnm,append=T)
      else {
        cat(i-1,1+dim(dm[[i]]),'\n',file=fnm,append=T)
        cat('-',colnames(dm[[i]]),sep=',','\n',file=fnm,append=T)
        write.table(dm[[i]],file=fnm,append=T,quote=F,col.names=F,sep=',')
      }
    }
  }
  i=.C(fnname,
   as.integer(a$nunits),as.integer(a$nsurveys),as.integer(a$nseasons),as.integer(a$nmethods),  #  N,T,nseasons,nmeth
   as.integer(t(a$det.data)), as.numeric(a$frq),                                               # det.dat, freq
   as.integer(a$nunitcov),as.numeric(t(a$unitcov)),as.character(unitcovnames),             # nsitecovs, sitecovs, sitecovnames
   as.integer(a$nsurvcov),as.numeric(t(a$survcov)),as.character(survcovnames),             # nsampcovs, sampcovs, sampcovnames
   as.integer(a$nsurveyseason),                                                                # nsurveys per season
   as.character(a$title),                                                                      # title
   as.character(a$unitnames),                                                                  # sitenames
   as.character(a$surveynames),                                                                # surveynames
   as.integer(dmdims),                                                                         # design matrices dimensions
   as.character(dms),                                                                          # design matrices
   as.character(c(rownames(dm[[1]]),rownames(dm[[2]]),rownames(dm[[3]]),rownames(dm[[4]]),rownames(dm[[5]]),rownames(dm[[6]]))), # real parm names
   as.numeric(initvals), as.numeric(fv),                                                     # init vals, fixed vals
   as.integer(neighbors), as.numeric(nbr_wgts),                                                # neighborhood matrix, weights
   as.character(outfile),                                                                        # output filename,
   as.character(miscopts),                                                                     # run options (eg., lmt=,maxfn=...)
   rvals)                                       #  return values: LL, betas, VC

  l=length(i); rvals=i[[l]]
  return(rvals)
}

normalize <- function(x) {      #   function to normalize a vector of covariates
  a=mean(x); b=x-a; bb=sqrt(sum(b*b)/(length(b)-1))
  return((x-a)/bb)
}

get_line<-function(v,line) unlist(strsplit(gsub("[ ()]+","@",v[line]),"@"))

check_conv_warn<-function(v,quiet=1){
  warn.conv=NULL; instns<-grep("Numerical convergence",v);   l<-length(instns)
  if(l>0) {
    warn.conv<-as.numeric(get_line(v,instns[1]+2)[2])
    if (quiet==0) cat("\nWarning: Numerical convergence may not have been reached. ",
                "Parameter esimates converged to approximately ", warn.conv," signifcant digits.\n\n")
  }
  return(warn.conv)
}

check_VC_warn<-function(v){
  warn.VC=NULL; instns<-grep("neg. std.err",v);  l<-length(instns)
  if(l>0) warn.VC<-TRUE
  return(warn.VC)
}

getbetas <- function(vv,x.dm) {
  if (length(vv)<1) return(NULL) else {
    vals=as.numeric(unlist(strsplit(gsub(" +",",",gsub(".+: +","",vv)),",")))
    vals=matrix(vals,ncol=2,byrow=T); rn=paste(1:ncol(x.dm))
    for (i in 1:ncol(x.dm)) {
      j=which.max(x.dm[,i]!="0")
      rn[i]=x.dm[j,i]; if (rn[i]=="1") rn[i]=paste0(rownames(x.dm)[j],".Int")
    }
    rownames(vals)=rn[1:nrow(vals)]
    return(data.frame(est=vals[,1],se=vals[,2]))
  }
}
get_coeff<-function(v,index,n.coeff,cov.names,b.string,b.count){
  out<-unlist(lapply(1:n.coeff,function(ii) {
    z<-get_line(v[index],ii)
    lz<-length(z)
    return(as.numeric(c(z[lz-1],z[lz],index[ii])))
  }))
  out<-t(array(out,dim=c(3,n.coeff)))
  colnames(out)<-c("est","se","index")
  rownames(out)<-paste(b.string,(b.count+1:n.coeff),"_",cov.names,sep="")
  return(out)
}

get_VC<-function(v,offset,n.coeff,cov.names){
  VC<-unlist(lapply(1:n.coeff,function(ii) {
    z<-get_line(v,offset+ii);  lz<-length(z)
    return(as.numeric(z[3:lz]))
  } ))
  VC<-t(array(VC,dim=c(n.coeff,n.coeff))); colnames(VC)<-cov.names; rownames(VC)<-cov.names
  return(VC)
}

get_real<-function(real,v,con.offset,n1,nunits,row.names,fixed,real.offset=0){
  instns<-grep(real,v,fixed=TRUE);
  l<-length(instns)
  if (l>0) {
    if(v[instns[1]+ n1+con.offset+1] == "") {
      out<-unlist(lapply(1:n1,function(jj){
        z<-get_line(v,instns[1]+jj+1); lz<-length(z)
        i=0+(z[lz]=="fixed")
        return(as.numeric(rep(c(z[lz-4-i],z[lz-3-i],z[lz-2-i],z[lz-i]),nunits)))
      }))
#    out<-t(array(out,dim=c(4,n1*nunits)))
    } else {
      out<-unlist(lapply(1:n1,function(jj){
        # index<-match((jj+real.offset),rownames(fixed))
        # fixedpar<-sum((1:jj+real.offset)%in%rownames(fixed))
        index<-match((jj+real.offset),fixed$idx)
        fixedpar<-sum((1:jj+real.offset)%in%fixed$idx)
        if(!is.na(index)){
          return(rep(c(fixed[index,1],NA,NA,NA),nunits))
        } else {
          return(
            unlist(lapply(1:nunits,function(ii){
              z<-get_line(v,instns[1]+con.offset+(jj-1-fixedpar)*nunits+ii+fixedpar); lz<-length(z)
              if(lz<5 || z[lz-4]==".") {
                return(c(NA,NA,NA,NA))
              } else {
                return(as.numeric(c(z[lz-4],z[lz-3],z[lz-2],z[lz])))
              }
            }))
          )
        }
      }))
    }
    out<-t(array(out,dim=c(4,n1*nunits))); colnames(out)<-c("est","se","lower","upper"); rownames(out)<-row.names
  }
  else out=NULL
  return(out)
}

### generic function
summary.pao<-function(object){
  #' Summary of occupancy data object
  #'
  #' Provides a brief summary of an object of class \code{pao}, which is the input data for occupancy analysis.
  #'
  #'@param object an object of class \code{pao} created by \code{readPao} or \code{createPao}.
  #'@export
  #'@return Provides a brief summary of the input object, including number of sites, surveys.
  #'
  #'@author Jim Hines
  #'
  #'@keywords pao

  x=matrix(as.numeric(unlist(object$det.data)),nrow(object$det.data))
  naiveocc=sum(0+(rowSums(x,na.rm=T)>0)*object$frq)/sum(object$frq)
  naiveR=sum(0+(rowSums(x==2,na.rm=T)>0)*object$frq)/sum(object$frq)
  v=c(object$nunits,object$nsurveys,object$nseasons,paste(object$nsurveyseason,collapse=','),
      object$nmethods,object$nunitcov,object$nsurvcov)
  cat("paoname=",object$paoname,'\ntitle=',object$title,'\nNaive occ=',naiveocc,'\n',sep="")
  if (naiveR>0) cat('naiveR   =',naiveR,'\n',sep='')
  names(v)=c("nunits","nsurveys","nseasons","nsurveyseason","nmethods","nunitcov","nsurvcov")
  print(v)
  cat('unit covariates  :',colnames(object$unitcov),'\n')
  cat('survey covariates:',colnames(object$survcov),'\n\n')
}

### generic function
is.pao <- function(object) {
  #' check if object is pao object
  #'@param object object to check
  #'@export
  #'@return T or F
  inherits(object,"pao")
}

### generic function
summary.occMod<-function(object,...){
  #' Summary of occupancy model
  #'
  #' Provides a brief summary of an object of class \code{occMod}, which is the result from fitting amodel using the \code{\link{occMod}} function.
  #'
  #'@param object an object of class \code{occMod} which is the result from fitting a model using the \code{\link{occMod}} function.
  #'@param ... possible additional arguments
  #'@export
  #'@return Provides a brief summary of the fitted model, including model name, AIC value and number of parameters.
  #'
  #'@author Darryl MacKenzie
  #'
  #'@keywords htest

  occmod<-object
  cat(paste("Model name=",occmod$modname,
            "\nAIC=",occmod$aic,
            "\n-2*log-likelihood=",occmod$neg2loglike,
            "\nnum. par=",occmod$npar,"\n",sep=""))
  if(length(occmod$warnings$conv)>0) cat(paste("Warning: Numerical convergence may not have been reached. ",
                                               "Parameter esimates converged to approximately ", occmod$warnings$conv,
                                               " signifcant digits.\n",sep=""))
  if(length(occmod$warnings$VC)>0) cat("Warning: Problem with estimation of VC covariance. Ignore all SE's and values in VC matrix.\n")
}

### generic function
summary.aicTable<-function(object,...){
  #' Summary of  \code{aic.table} object
  #'
  #' Provides a brief summary of an object of class \code{aicTable}, which is the result of
  #' the function \code{\link{createAicTable}}.
  #'
  #'@param object an object of class \code{aicTable} which is the result of the function
  #' \code{\link{createAicTable}}
  #'@param ... possible additional arguments
  #'@export
  #'@return Produces a summary of the model comparison using AIC as a data frame that
  #' includes the model name, relative difference in AIC, model weight, number of parameters,
  #' twice the negative log-likelihood and any warnings.
  #'
  #'@author Darryl MacKenzie
  #'
  #'@keywords htest
  aictable<-object
  summ<-data.frame(Model=aictable$table$Model,
                   DAIC=as.numeric(format(aictable$table$DAIC,digits=2,nsmall=2)),
                   wgt=as.numeric(format(aictable$table$wgt,digits=2,nsmall=2)),
                   npar=as.numeric(format(aictable$table$npar,digits=2,nsmall=0)),
                   neg2ll=as.numeric(format(aictable$table$neg2ll,digits=2,nsmall=2)),
                   warn.conv=as.numeric(format(aictable$table$warn.conv,digits=2,nsmall=2)),
                   warn.VC=as.numeric(format(aictable$table$warn.VC,digits=2,nsmall=0))
  )
  return(summ)
}

createAicTable<-function(mod.list, use.aicc=FALSE){
  #' Create AIC table
  #'
  #' Create an AIC table of results from a list of previously run models
  #'
  #'@param mod.list a list of objects of class \code{occMod}, i.e., the output from
  #' the function \code{\link{occMod}}.
  #'@param use.aicc logical argument for whether small-sample adjustment should be used
  #' with AIC. Effective sample set to number of units.
  #' ... OR... user-defined value for effective sample size.  So, possibilities are:
  #' \itemize{
  #' \item{\code{0 or FALSE}}{: Use AIC}
  #' \item{\code{1 or TRUE}} {: Use number of sites for eff. sample size to calculate AICc}
  #' \item{\code{n}} {: Use \code{n} for eff. sample size to calculate AICc}
  #' }
  #'@export
  #'@return Returns a list with the following components:
  #'  \itemize{
  #'   \item{\code{table}}{ - A data frame that includes:
  #'   \itemize{
  #'     \item{the name}
  #'     \item{AIC value}
  #'     \item{relative difference in AIC}
  #'     \item{AIC model weight}
  #'     \item{-2*log-likelihood}
  #'     \item{number of parameters}
  #'     \item{any warnings for each model}
  #'   }
  #'   Table entries are ranked according to AIC.
  #'   }
  #'   \item{\code{models}}{ - List of models supplied in \code{mod.list}, ordered by AIC}
  #'   \item{\code{ess}}{ - Effective sample size for AICc calculation (0=use AIC, >0=eff. samp. size)}
  #'  }
  #'@author Darryl MacKenzie
  #'
  #'@seealso \code{\link{occMod}}

  table<-unlist(lapply(mod.list,function(xx){
    c(xx$aic,xx$neg2loglike,xx$npar,ifelse(length(xx$warnings$conv)>0,xx$warnings$conv,0),!is.null(xx$warnings$VC))
  }))
  table<-as.data.frame(t(array(table,dim=c(5,length(mod.list)))))
  colnames(table)<-c("AIC","neg2ll","npar","warn.conv","warn.VC")
  table<-data.frame(Model=unlist(lapply(mod.list,function(xx) xx$modname)),table)
  k=table$npar; n=0; if (use.aicc==1) n=mod.list[[1]]$data$nunits;
  if (use.aicc>1) n=use.aicc
  aicc=table$AIC+2*k*(k+1)/(n-k-1)
  if (use.aicc>0) table$AIC=aicc
  table$DAIC<-table$AIC-min(table$AIC)
  table$modlike<-exp(-table$DAIC/2)
  table$wgt<-round(table$modlike/sum(table$modlike),4)
  table$modlike=round(table$modlike,4)
  table$DAIC=round(table$DAIC,4)
  order<-order(table$DAIC)
  if (use.aicc) {colnames(table)[2]="AICc"; colnames(table)[7]="DAICc"}
  #  table<-cbind(as.character(mod.list),table)
  result<-list(table=table[order,],models=mod.list[order], ess=n)
  class(result)<-"aicTable"
  return(result)
}

# predict.p<-function(occmod,newdata){
#   pred.p.mat<-occmod.matrix(occmod$model$p,newdata)
#   p<-plogis(pred.p.mat%*%occmod$beta$p.coeff[,"est"])
#   logit.var<-apply(pred.p.mat,1,function(xx) xx%*%occmod$beta$p.VC%*%xx)
#   se<-sqrt(logit.var)*(p*(1-p))
#   return(data.frame(p=p,se=se))
# }
#
# predict.psi<-function(occmod,newdata){
#   pred.psi.mat<-occmod.matrix(occmod$model$psi,newdata)
#   psi<-plogis(pred.psi.mat%*%occmod$real$psi.coeff[,"est"])
#   logit.var<-apply(pred.psi.mat,1,function(xx) xx%*%occmod$beta$psi.VC%*%xx)
#   se<-sqrt(logit.var)*(psi*(1-psi))
#   return(data.frame(psi=psi,se=se))
# }

### generic function
predict.occMod<-function(object,newdata,param=NULL,conf=0.95){
  #' Predict real parameter values.
  #'
  #' Predict real parameter values from a model that has been fitted to the data, using the
  #' variable or covariates values that are supplied by the user. \bold{Currently has only been coded
  #' to predict probabilities and checked for the static, or single-season, occupancy model.}
  #'
  #' Note that can not make predictions directly from models that have been imported from a PRESENCE
  #' project using the \code{\link{importProject}} function.
  #'
  #'@param object and object of class \code{occMod} which is the model that has been fitted to
  #'  the data from which predictions are required. The result of a call to \code{\link{occMod}}.
  #'@param newdata a data frame containing the values of the variables for which predictions are
  #' required. Variable names must exactly match those used in the fitted model. May contain
  #' additional variables not used in the fitted model.
  #'@param param a string containing the name of the real parameter for which predictions
  #' are required, e.g., \code{type="psi"}.
  #'@param conf level for confidence interval (may be vector valued).
  #'
  #'@export
  #'
  #'@return returns a data frame containing the predicted value, standard error and limits of a
  #' 95\% confidence interval.
  #'
  #'@examples
  #' psi.real=predict(mod1,salmdr.data$survcov,param="psi",conf=0.95); print(head(psi.real))
  #'@author Darryl MacKenzie
  #'@seealso \code{\link{occMod}}
  #'

#  if(!(param%in%c("psi","p"))) stop("Only predictions for psi or p are supported at the moment!")
  occmod<-object
  if(is.null(occmod$model)) stop("Predictions not supported for imported projects.")
  coeff<-get(param,occmod$beta)[,"est"]
  VC<-get(paste0(param,".VC"),occmod$beta)
  pred.mat<-model.matrix(get(param,occmod$model),newdata)
  logit.prob<-pred.mat%*%coeff;
  logit.var<-apply(pred.mat,1,function(xx) xx%*%VC%*%xx)

  alpha<-(1-conf)/2; z<--qnorm(alpha)
  logit.lower<-sapply(z,function(zz) logit.prob-zz*sqrt(logit.var))
  if (length(logit.lower)==1) names(logit.lower)<-paste("lower",conf,sep="_") else colnames(logit.lower)<-paste("lower",conf,sep="_")

  logit.upper<-sapply(z,function(zz) logit.prob+zz*sqrt(logit.var))
  if (length(logit.upper)==1) names(logit.upper)<-paste("upper",conf,sep="_") else colnames(logit.upper)<-paste("upper",conf,sep="_")

  prob<-plogis(logit.prob);  se<-sqrt(logit.var)*(prob*(1-prob))
  lower<-plogis(logit.lower);  upper<-plogis(logit.upper)

  v=data.frame(est=prob,se=se,lower,upper);
  return(v)
}

modAvg<-function(aic.tab,param=NULL,index=1:nrow(aic.tab$table),conf=0.95,predict=FALSE,newdata){
  #' Calculate model averaged estimates
  #'
  #' Calculates model averaged estimates of real parameters stored in a \code{aicTable} object;
  #' the result of \code{\link{createAicTable}}. Does not perform model averaging of the
  #' regression coefficients, or beta parameters. Can be used to model average predicted
  #' probabilities too.
  #'
  #'@param aic.tab an \code{aicTable} object; the result of \code{\link{createAicTable}}.
  #'@param param a string with the name of the real parameter to model average, e.g.,
  #'\code{param = "psi"}.  Presently, only works for parameters "psi" and "p".
  #'@param index a vector of integer values indexing which ranked models are to be model
  #'averaged across. Can be used if only model averaging across a subset of models.
  #'Note that AIC weights are recalcuated such that their sum = 1 for the subset.
  #'@param conf level for confidence interval (may be vector valued).
  #'@param predict if \code{FALSE} (default) then model averages observed real estimates.
  #' Otherwise model averages predicted values for \code{param} from the supplied \code{newdata}.
  #'@param newdata a data frame containing the values of the covariates for which
  #'predictions are required. Variable names must exactly match those used in the data.
  #'May contain varaibles that are not used in some models.
  #'@export
  #'@return returns a data frame containing the model averaged estimates, standard errors
  #'and limts of a 95\% confidence interval.
  #'
  #'@author Darryl MacKenzie
  #'@seealso \code{\link{createAicTable}}, \code{\link{predict.occMod}}
  #'
  #'@keywords htest

  #make an internal copy, mainly ofr only using a subset of full models
  aic.table<-list(table=aic.tab$table[index,],models=aic.tab$models[index])
  class(aic.table)<-"aic.table"
  aic.table$table$DAIC<-aic.table$table$AIC-min(aic.table$table$AIC)
  aic.table$table$modlike<-exp(-aic.table$table$DAIC/2)
  aic.table$table$wgt<-aic.table$table$modlike/sum(aic.table$table$modlike)

  if(!predict){ ## model averaging estimates from the observed data
    est<-array(unlist(lapply(aic.table$models, function(xx) get(param,xx$real)[,"est"])),
               dim=c(nrow(get(param,aic.table$models[[1]]$real)),length(aic.table$models)))
    se<-array(unlist(lapply(aic.table$models, function(xx) get(param,xx$real)[,"se"])),
              dim=c(nrow(get(param,aic.table$models[[1]]$real)),length(aic.table$models)))
  } else {
    pred<-lapply(aic.table$models, function(xx) predict(xx,newdata,param=param))
    est<-array(unlist(lapply(pred, function(xx) xx[,"est"])),
               dim=c(nrow(newdata),length(aic.table$models)))
    se<-array(unlist(lapply(pred, function(xx) xx[,"se"])),
              dim=c(nrow(newdata),length(aic.table$models)))

#     if(predict=="p") {
#       pred<-lapply(aic.table$models, function(xx) predict.p(xx,newdata))
#       est<-array(unlist(lapply(pred, function(xx) xx[,"p"])),
#                  dim=c(nrow(newdata),length(aic.table$models)))
#       se<-array(unlist(lapply(pred, function(xx) xx[,"se"])),
#                 dim=c(nrow(newdata),length(aic.table$models)))
#     }
#     if(predict=="psi") {
#       pred<-lapply(aic.table$models, function(xx) predict.psi(xx,newdata))
#       est<-array(unlist(lapply(pred, function(xx) xx[,"psi"])),
#                  dim=c(nrow(newdata),length(aic.table$models)))
#       se<-array(unlist(lapply(pred, function(xx) xx[,"se"])),
#                 dim=c(nrow(newdata),length(aic.table$models)))
#     }
  }

  ma<-est%*%aic.table$table$wgt
  ma.se<-sqrt((se^2+(est-as.vector(ma))^2)%*%aic.table$table$wgt)

  ## calculate confidence intervals. Presuming dealing with probabilities here
  logit.est<-log(ma/(1-ma))

  logit.se<-ma.se # not on logit scale yet, just taking a copy and used below
  on.bounds<-ma%in%1|ma%in%0
  #  on.bounds<-ma==1|ma==0
  logit.se[on.bounds]<-0
  logit.se[!on.bounds]<-ma.se[!on.bounds]/(ma[!on.bounds]*(1-ma[!on.bounds]))

  alpha<-(1-conf)/2;  z<--qnorm(alpha)
  lower<-sapply(z,function(zz) logit.est-zz*logit.se)
  colnames(lower)<-paste("lower",conf,sep="_")

  upper<-sapply(z,function(zz) logit.est+zz*logit.se)
  colnames(upper)<-paste("upper",conf,sep="_")

  lower<-plogis(lower)
  upper<-plogis(upper)

  result<-data.frame(est=ma,se=ma.se,lower,upper)
  if(!predict) {
    rownames(result)<-rownames(get(param,aic.table$models[[1]]$real))
  } else {
    rownames(result)<-rownames(newdata)
  }
  return(result)
}

createPao<-function(data,unitcov=NULL,survcov=NULL,nsurveyseason=ncol(data),nmethods=1,frq=rep(1,nrow(data)),
                     title="PRESENCE Analysis",unitnames=paste0("unit",1:nrow(data)),paoname="pres.pao") {
#'  Create a pao data object.
#'
#' \code{createPao} Creates a pao data object from R variables.
#'
#'@param data data frame containing the detection (=1) and nondetection (=0) data
#'   with observations from a sampling unit along a row, and surveys in the columns.
#'   Missing values are allowed (=\code{NA}).
#'@param frq vector of frequencies of each detection history (optional)
#'@param unitcov data frame for unit-specific covariates (e.g., habitat,
#'   elevation, etc.) with 1 row per unit (so \code{nrow(unitcov) == nrow(data)}).
#'   Each column of the data frame should contain a single covariate, with the
#'   column name used as the covariate name. Covariatecolumn names shoud not contain
#'   spaces, and avoid using all uppercase as  some terms are used by the package
#'   to define common effects used in models, e.g., SURVEY. If there are no
#'   unit-specific covariates then \code{unitcov=NULL}. Note that additional
#'   unit-specific covariates can be supplied when fitting a model with
#'   \code{\link{occMod}}.
#'@param survcov data frame for survey-specific covariates (e.g., survey time,
#'   weather conditions, observer, etc.) with 1 row per survey (so
#'   \code{nrow(survcov) == nrow(data)*ncol(data)}). Each column of the data
#'   frame should contain a single covariate, with the column name used as the
#'   covariate name. Covariate/column names shoud not contain spaces, and avoid
#'   using all uppercase as  some terms are used by the package to define common
#'   effects used in models, e.g., SURVEY. If there are no survey-specific
#'   covariates then \code{survcov=NULL}. Note that additional survey-specific
#'   covariates can be supplied when fitting a model with \code{\link{occMod}}.
#'   Note that the format of survey-specific covariates in \code{survcov} is
#'   different from that used in PRESENCE (i.e., long vs flat format respectively.
#'@param nsurveyseason vector of the number of surveys per season (e.g., if
#'   three seasons, with 4 surveys per season, \code{nsurveyseason=c(4,4,4)}).
#'   The default assumes that \code{data} is for a single-season or static
#'   occupancy model.
#'@param nmethods integer number of methods per survey
#'@param title (optional) a label included in the PRESENCE output files,
#'   of little practical consequence.
#'@param unitnames a vector of string containing labels for each sampling
#'   unit. Included in the results as row names for the real parameter estimates
#'   of the \code{\link{occMod}} object.
#'@param paoname string containing a name for the temporary pao
#'@export
#'@return list of class "\code{pao}" containing the data required to fit an
#'  occupancy model using \code{\link{occMod}}. The list contains the following
#'  objects:
#'   \item{nunits}{number of units or sampling units, i.e., \code{nrow(data)}.}
#'   \item{nsurveys}{total number of surveys, i.e., \code{ncol(data)}.}
#'   \item{nseasons}{number of seasons of data.}
#'   \item{det.data}{detection/nondetection data, i.e., \code{data}.}
#'   \item{unitcov}{see above. If no covariates are supplied a vector of 1's
#'   named TEMP is added.}
#'   \item{survcov}{see above. A factor called SURVEY is automatically added to
#'   allow survey-specific detection probabilities (but equal across units).}
#'   \item{nsurveyseason}{see above.}
#'   \item{title}{see above.}
#'   \item{unitnames}{see above.}
#'   \item{surveynames}{a vector of strings providing labels for the survey occasions. Of little practical consequence.}
#'   \item{paoname}{see above. Equals \code{f} for \code{createPao}.}
#'
  nunits<-nrow(data); nsurveys<-ncol(data); nseasons<-length(nsurveyseason)
  if(is.null(unitcov)) unitcov<-data.frame(TEMP=rep(1:nunits)) # adding temp variable
  nunitcov<-ncol(unitcov); srvcov=data.frame(SURVEY=as.factor(rep(1:nsurveys,each=nunits)))
  if (is.null(survcov)) survcov=srvcov else survcov=data.frame(survcov,srvcov)
  nsurvcov<-ncol(survcov)
  surveynames<-paste0(
    unlist(lapply(1:nseasons,function(xx) {rep(xx,each=nsurveyseason[xx])})),  "-",
    unlist(lapply(1:nseasons,function(xx) {1:nsurveyseason[xx]}))  )
  if (is.null(rownames(data))) rownames(data)=unitnames

  ## nunits, nsurveys, nseasons, d,
  ## nunitcov, unitcov,nsurvcov,survcov
  ## nsurveyseason, title, unitnames, surveynames, paoname
  pao<-list(nunits=nunits,nsurveys=nsurveys,nseasons=nseasons,nmethods=nmethods,det.data=data,
            nunitcov=nunitcov,unitcov=unitcov,nsurvcov=nsurvcov,survcov=survcov,
            nsurveyseason=nsurveyseason,title=title,unitnames=unitnames,
            surveynames=surveynames,paoname=paoname,frq=frq)
  class(pao)<-c("pao")
  return(pao)

}

extract_model<-function(model.list,params,e,model.type){
  lapply(1:length(model.list),function(xx){
    form<-model.list[[xx]]
    #check lh-side is valid
    if(!(deparse(form[[2]]) %in% params)){
      stop(paste0(deparse(form[[2]])," is not a valid parameter for model.type=",model.type,"!"))
    }
    assign(deparse(form[[2]]),nlme::getCovariateFormula(form),envir=e)
  })
}

read_output<-function(outfile, type, aux=NULL){
  # aux is a list containing the design matrices for the model, and other info
  # aux = NULL if reading importing PRESENCE results run from outside of R

  ### extract results from output file
  v=readLines(outfile); l=length(v);

  ## find AIC value for each model
  instns=grep("-2log(likelihood)",v,fixed=TRUE)
  z<-get_line(v,instns[1])
  neg2loglike<-as.numeric(z[4])
  z<-get_line(v,instns[1]+1)
  aic<-as.numeric(z[3])

  # find number of parameters
  instns=grep("Number of parameters",v,fixed=TRUE)
  z<-get_line(v,instns[1])
  npar<-as.numeric(z[length(z)])

  ## need to split off here for the different model types
  if(type=="so"){
    estimates<-get_so(v,aux,npar)
  }

  if(type=="do.1"){
    estimates<-get_do_1(v,aux,npar)
  }

  if(type=="do.4"){
    estimates<-get_do_4(v,aux,npar)
  }

  ##### check for warnings
  warn.conv<-check_conv_warn(v)
  warn.VC<-check_VC_warn(v)

  return(list(
    neg2loglike=neg2loglike,
    npar=npar, aic=aic,
    beta=estimates$beta,
    real=estimates$real,
    warnings=list(conv=warn.conv,VC=warn.VC)
  ))
}

importProject<-function(pa3){
  #' Import a PRESENCE project into R.
  #'
  #' Imports the results stored in a PRESENCE project folder into R as an \code{aic.table}
  #' object. Note that in this package version the code has only been developed to import
  #' static occupancy/single-season models, and the first and fourth parameterisation of the
  #' dynamic occupancy/multi-season models. Because design matrices in PRESENCE are constructed
  #' manually, and can be done in multiple ways, no attempt has been made to second-guess the
  #' intended model stucture with respect to covariates. Therefore predictions can not be made
  #' from imported models (unless manually reconfigured), although model averaged estimates of
  #' real parameters of the surveyed units can be obtained.
  #'
  #'@param pa3 a string containing the path (if outside of the current working directory)
  #' and file name of the PRESENCE project summary file with the extension .pa3. e.g.,
  #' \code{import.project("bird_project/bird_project.pa3")}
  #'@export
  #'@return  Returns an object of class \code{aic.table}, the same as the result of
  #'the \code{\link{createAicTable}} function. The only difference being that the
  #'\code{model} component of the stored \code{\link{occMod}} objects = NULL. i.e.,
  #'the model structure is not imported and left undefined.
  #'
  #'@author Darryl MacKenzie.
  #' @importFrom utils read.table
  #'
  #'@seealso \code{\link{createAicTable}},\code{\link{occMod}},
  #'\code{\link{predict.occMod}},\code{\link{modAvg}}
  #'
  #'@keywords models
  #'
  olddir<-getwd()
  setwd(dirname(pa3))

  v=readLines(pa3)
  ind<-grep("DATAFILE",v,fixed=TRUE)
  paofile<-unlist(strsplit(v[ind],":"))[2]
  #  paofile<-file.path(dirname(pa3),paofile)
  pao<-readPao(paofile)

  i=grep("^Model.AIC",v)
  table<-read.table(pa3,sep="\t",skip=i-1,header=TRUE)

  models<-lapply(table$Model,function(modname){
    outfile<-gsub(",","_",modname,fixed=TRUE)
    outfile<-gsub(" ","_",outfile,fixed=TRUE)
    outfile<-gsub("+","P",outfile,fixed=TRUE)
    outfile<-gsub("*","X",outfile,fixed=TRUE)
    outfile<-paste0("pres_",outfile,".out")

    # peek at outfile to determine type
    v<-readLines(outfile); cat('reading',outfile,'\n')
    ind<-grep("==>model",v,fixed=TRUE)
    z<-unlist(strsplit(v[ind],"="))
    modtype<-as.numeric(z[length(z)])
    type=""

    if (length(grep('Custom Model:',v))>0) type="so" # if(modtype==100) {type="so"}
    if (length(grep('^Multi-season model - ',v))>0) {
        type="do.1" # if(modtype==200) {type="do.1"}
        if (length(grep(", eps=1-gam param",v)>0)) type="do.4" # if(modtype==240) {type="do.4"}
    }

    if(type=="") stop(paste0(modname, " is of an unsupported type in the current package. Stopping import."))

    import<-read_output(outfile,type=type,aux=NULL)
    result<-list(modname=modname,
                 model=NULL,
                 data=pao,outfile=outfile,
                 neg2loglike=import$neg2loglike,
                 npar=import$npar, aic=import$aic,
                 beta=import$beta,
                 real=import$real,
                 warnings=import$warnings)
    class(result)<-c("occMod",type)
    return(result)
  })

  setwd(olddir)


  return(createAicTable(models))
}


get_so<-function(v,aux,npar){
  ### get coefficients
  instns<-grep("Untransformed",v,fixed=TRUE)
  extract<-v[(instns[1]+3):(instns[1]+2+npar)]
  if(is.null(aux)){ # so importing results from a non R-based run
    # get number of parameters for psi and p
    instns<-grep("Matrix 1",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.psi<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1

    instns<-grep("Matrix 2",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.p<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1

    # make indices for psi and p
    index.psi<-1:npar.psi
    index.p<-(npar.psi+1):npar

    # get covariate names
    covnames<-unlist(lapply(1:length(extract),function(zz){
      get_line(extract,zz)[2]
    }))
    covnames.psi<-covnames[index.psi]
    covnames.p<-covnames[index.p]

    # get pao info
    instns<-grep("==>i",v,fixed=TRUE)
    z<-unlist(unlist(strsplit(v[instns[1]],"=")))
    file<-z[length(z)]
    temp.pao<-readPao(file)

    # get fixed parameter info
    instns<-grep("fixed",v,fixed=TRUE)
    if(length(instns)>0){
      fxd<-data.frame()
      for(ii in 1:length(instns)){
        rbind(fxd,get_line(v[instns][ii]))
      }
      fixed<-matrix(fxd[,4],ncol=1)
      rownames(fixed)<-fxd[,2]
    } else{
      fixed=NULL
    }

  } else { # get this info from aux
    # get number of parameters for psi and p
    npar.psi<-ncol(aux$psi_mat)
    npar.p<-ncol(aux$p_mat)

    # make indices for psi and p
    index.psi<-grep(".psi.",extract,fixed=TRUE)
    index.p<-grep(".p.",extract,fixed=TRUE)

    # get covariate names
    covnames.psi<-colnames(aux$psi_mat)
    covnames.p<-colnames(aux$p_mat)

    # get pao info
    temp.pao<-aux$temp.pao

    # get fixed parameter info
    fixed=aux$fixed
  }

  psi.coeff<-get_coeff(v=extract,index=index.psi,n.coeff=npar.psi,
                       cov.names=covnames.psi,b.string="a",b.count=0)
  p.coeff<-get_coeff(v=extract,index=index.p,n.coeff=npar.p,
                     cov.names=covnames.p,b.string="b",b.count=0)
  ### get VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  names<-c(rownames(psi.coeff),rownames(p.coeff))
  VC<-get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC<-VC[psi.coeff[,"index"],psi.coeff[,"index"]]

  p.VC<-VC[p.coeff[,"index"],p.coeff[,"index"]]

  ## Get real parameter estimates from output
  i=grep("<[Pp]si>",v); tmps=gsub(">.+",">",gsub(".+<","<",v[i]))
  psi.est<-get_real(real=tmps,v=v,con.offset=1,
                    n1=1,nunits=temp.pao$nunits,
                    row.names=temp.pao$unitnames,
                    fixed=fixed,real.offset=0)
  i=grep("<p1>",v); if (length(i)<1) i=grep("<p[1",v,fixed=T);
  tmps=gsub(">.+",">",gsub(".+<","<",v[i]))
  p.est<-get_real(real=tmps,v=v,con.offset=1,
                  n1=temp.pao$nsurveys,nunits=temp.pao$nunits,
                  row.names=paste(rep(temp.pao$unitnames,temp.pao$nsurveys),"_",
                                  rep(temp.pao$surveynames,each=temp.pao$nunits),
                                  sep=""),
                  fixed=fixed,real.offset=temp.pao$nseasons+1)
  psi_c.est<-get_real(real="Psi-conditional",v=v,con.offset=2,
                      n1=1,nunits=temp.pao$nunits,
                      row.names=temp.pao$unitnames,
                      fixed=fixed,real.offset=0)

  return(list(
    beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
              p=as.data.frame(p.coeff),p.VC=p.VC, VC=VC),
    real=list(psi=as.data.frame(psi.est),
              p=as.data.frame(p.est),
              psi_c=as.data.frame(psi_c.est))
  ))

}

get_do_1<-function(v,aux,npar){
  ### get coefficients
  instns<-grep("Untransformed",v,fixed=TRUE)
  extract<-v[(instns[1]+3):(instns[1]+2+npar)]
  if(is.null(aux)){ # so importing results from a non R-based run
    # get number of parameters for each parameter type
    instns<-grep("Matrix 1",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.psi<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1
    nrow.psi<-as.numeric(unlist(strsplit(unlist(strsplit(z[length(z)-1],"="))[2],",")))-1

    instns<-grep("Matrix 2",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.gamma<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1
    nrow.gamma<-as.numeric(unlist(strsplit(unlist(strsplit(z[length(z)-1],"="))[2],",")))-1

    instns<-grep("Matrix 3",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.epsilon<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1
    nrow.epsilon<-as.numeric(unlist(strsplit(unlist(strsplit(z[length(z)-1],"="))[2],",")))-1

    instns<-grep("Matrix 4",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.p<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1

    # make indices for psi and p
    index.psi<-1:npar.psi
    index.gamma<-(npar.psi+1):(npar.psi+npar.gamma)
    index.epsilon<-(npar.psi+npar.gamma+1):(npar.psi+npar.gamma+npar.epsilon)
    index.p<-(npar.psi+npar.gamma+npar.epsilon+1):(npar)

    # get covariate names
    covnames<-unlist(lapply(1:length(extract),function(zz){
      get_line(extract,zz)[2]
    }))
    covnames.psi<-covnames[index.psi]
    covnames.gamma<-covnames[index.gamma]
    covnames.epsilon<-covnames[index.epsilon]
    covnames.p<-covnames[index.p]

    # get pao info
    instns<-grep("==>i",v,fixed=TRUE)
    z<-unlist(unlist(strsplit(v[instns[1]],"=")))
    file<-z[length(z)]
    temp.pao<-readPao(file)

    # get fixed parameter info
    instns<-grep("fixed",v,fixed=TRUE)
    if(length(instns)>0){
      fxd<-data.frame()
      for(ii in 1:length(instns)){
        rbind(fxd,get_line(v[instns][ii]))
      }
      fixed<-matrix(fxd[,4],ncol=1)
      rownames(fixed)<-fxd[,2]
    } else{
      fixed=NULL
    }

    # labelling quirks
    gamma.string<-"<gam1>"
    epsilon.string<-"<eps1>"

    # create psi, gamma, and epsilon matrices for derived estimates
    #     instns<-grep("Matrix 1",v,fixed=TRUE)
    #     psi_mat<-data.frame()
    #     for (ii in 1:nrow.psi){
    #       z<-get_line(v,instns[1]+1+ii)
    #       temp<-unlist(lapply(z[-1],function(xx){
    #         if(suppressWarnings(is.na(as.numeric(xx)))){
    #           # must be a covariate name
    #           temp.pao$unitcov[,xx]
    #         } else{
    #           rep(as.numeric(xx),temp.pao$nunits)
    #         }
    #       }))
    #       temp<-array(temp,dim=c(temp.pao$nunits,npar.psi))
    #       psi_mat<-rbind(psi_mat,temp)
    #     }
    #     psi_mat<-as.matrix(psi_mat)
    #
    #     instns<-grep("Matrix 2",v,fixed=TRUE)
    #     gamma_mat<-data.frame()
    #     for (ii in 1:nrow.gamma){
    #       z<-get_line(v,instns[1]+1+ii)
    #       temp<-unlist(lapply(z[-1],function(xx){
    #         if(suppressWarnings(is.na(as.numeric(xx)))){
    #           # must be a covariate name
    #           temp.pao$unitcov[,xx]
    #         } else{
    #           rep(as.numeric(xx),temp.pao$nunits)
    #         }
    #       }))
    #       temp<-array(temp,dim=c(temp.pao$nunits,npar.gamma))
    #       gamma_mat<-rbind(gamma_mat,temp)
    #     }
    #     gamma_mat<-as.matrix(gamma_mat)
    #
    #     instns<-grep("Matrix 3",v,fixed=TRUE)
    #     epsilon_mat<-data.frame()
    #     for (ii in 1:nrow.epsilon){
    #       z<-get_line(v,instns[1]+1+ii)
    #       temp<-unlist(lapply(z[-1],function(xx){
    #         if(suppressWarnings(is.na(as.numeric(xx)))){
    #           # must be a covariate name
    #           temp.pao$unitcov[,xx]
    #         } else{
    #           rep(as.numeric(xx),temp.pao$nunits)
    #         }
    #       }))
    #       temp<-array(temp,dim=c(temp.pao$nunits,npar.epsilon))
    #       epsilon_mat<-rbind(epsilon_mat,temp)
    #     }
    #     epsilon_mat<-as.matrix(epsilon_mat)

  } else { # get this info from aux
    # get number of parameters for psi and p
    npar.psi<-ncol(aux$psi_mat)
    npar.gamma<-ncol(aux$gamma_mat)
    npar.epsilon<-ncol(aux$epsilon_mat)
    npar.p<-ncol(aux$p_mat)

    # make indices for psi and p
    index.psi<-grep(".psi.",extract,fixed=TRUE)
    index.gamma<-grep(".gamma.",extract,fixed=TRUE)
    index.epsilon<-grep(".epsilon.",extract,fixed=TRUE)
    index.p<-grep(".p.",extract,fixed=TRUE)

    # get covariate names
    covnames.psi<-colnames(aux$psi_mat)
    covnames.gamma<-colnames(aux$gamma_mat)
    covnames.epsilon<-colnames(aux$epsilon_mat)
    covnames.p<-colnames(aux$p_mat)

    # get pao info
    temp.pao<-aux$temp.pao

    # get fixed parameter info
    fixed=aux$fixed

    # labelling quirks
    gamma.string<-"<gamma1>"
    epsilon.string<-"<epsilon1>"

    ## get design matrices for prediction
    psi_mat<-aux$psi_mat
    gamma_mat<-aux$gamma_mat
    epsilon_mat<-aux$epsilon_mat
  }
  cat('\n****************',v[8],' ************\n')
  psi.coeff<-get_coeff(v=extract,index=index.psi,n.coeff=npar.psi,
                       cov.names=covnames.psi,b.string="a",b.count=0)
  gamma.coeff<-get_coeff(v=extract,index=index.gamma,n.coeff=npar.gamma,
                         cov.names=covnames.gamma,b.string="b",b.count=0)
  epsilon.coeff<-get_coeff(v=extract,index=index.epsilon,n.coeff=npar.epsilon,
                           cov.names=covnames.epsilon,b.string="c",b.count=0)
  p.coeff<-get_coeff(v=extract,index=index.p,n.coeff=npar.p,
                     cov.names=covnames.p,b.string="d",b.count=0)
  ### get VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  names<-c(rownames(psi.coeff),rownames(gamma.coeff),rownames(epsilon.coeff),rownames(p.coeff))
  VC<-get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC<-VC[psi.coeff[,"index"],psi.coeff[,"index"]]
  gamma.VC<-VC[gamma.coeff[,"index"],gamma.coeff[,"index"]]
  epsilon.VC<-VC[epsilon.coeff[,"index"],epsilon.coeff[,"index"]]
  p.VC<-VC[p.coeff[,"index"],p.coeff[,"index"]]

  ## Get real parameter estimates from output
  i=grep('std.error',v)[1]; vv=v[i:length(v)]   #  get lines after design matrix stuff...

  i=grep('^psi1 .+:',vv); v1=gsub('.+: +','',vv[i]); v2=gsub('[ -]+',',',v1); v3=unlist(strsplit(v2,','))
  psi.est=matrix(as.numeric(v3),ncol=4,byrow=T);
  if (nrow(psi.est)==1) psi.est=psi.est[rep(1,temp.pao$nunits),]
  rownames(psi.est)=temp.pao$unitnames
  colnames(psi.est)=c('est','se','lower','upper')
  #psi.est<-get_real(real="<psi1>",v=v,con.offset=1,n1=1,nunits=temp.pao$nunits,row.names=temp.pao$unitnames,fixed=fixed,real.offset=0)

  i=grep('^gam\\(.+\\)',vv); v1=gsub('.+: +','',vv[i]); v2=gsub('[ -]+',',',v1); v3=unlist(strsplit(v2,','))
  gamma.est=matrix(as.numeric(v3),ncol=4,byrow=T); colnames(gamma.est)=colnames(psi.est)
  if (nrow(gamma.est)==(temp.pao$nseasons-1)) gamma.est=gamma.est[rep(1:nrow(gamma.est),each=temp.pao$nunits),]
  rownames(gamma.est)=paste0(temp.pao$unitnames,"_",rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits))
  #gamma.est<-get_real(real=gamma.string,v=v,con.offset=1,n1=temp.pao$nseasons-1,nunits=temp.pao$nunits,
  #                    row.names=paste0(temp.pao$unitnames,"_",rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits)),fixed=fixed,real.offset=1)

  i=grep('^eps\\(.+\\)',vv); v1=gsub('.+: +','',vv[i]); v2=gsub('[ -]+',',',v1); v3=unlist(strsplit(v2,','))
  epsilon.est=matrix(as.numeric(v3),ncol=4,byrow=T); colnames(epsilon.est)=colnames(psi.est)
  if (nrow(epsilon.est)==(temp.pao$nseasons-1)) epsilon.est=epsilon.est[rep(1:nrow(epsilon.est),each=temp.pao$nunits),]
  rownames(epsilon.est)=paste0(temp.pao$unitnames,"_",rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits))
  #epsilon.est<-get_real(real=epsilon.string,v=v,con.offset=1,n1=temp.pao$nseasons-1,nunits=temp.pao$nunits,
  #                      row.names=paste0(temp.pao$unitnames,"_",rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits)),
  #                      fixed=fixed,real.offset=temp.pao$nseasons)

  i=grep('^P\\[',vv); v1=gsub('.+: +','',vv[i]); v2=gsub('[ -]+',',',v1); v3=unlist(strsplit(v2,','))
  p.est=matrix(as.numeric(v3),ncol=4,byrow=T); colnames(p.est)=colnames(psi.est)
  rownames(p.est)=paste0(rep(temp.pao$unitnames,temp.pao$nsurveys),"_",rep(temp.pao$surveynames,each=temp.pao$nunits))
  #p.est<-get_real(real="<P[1-1]>",v=v,con.offset=1,n1=temp.pao$nsurveys,nunits=temp.pao$nunits,
  #                row.names=paste0(rep(temp.pao$unitnames,temp.pao$nsurveys),"_",rep(temp.pao$surveynames,each=temp.pao$nunits)),
  #                fixed=fixed,real.offset=2*(temp.pao$nseasons-1)+1)

  instns<-grep("psi2,psi3,",v,fixed=TRUE); ii <- instns[1]+3;
  temp=data.frame()
  repeat{
    z<-get_line(v,ii); ll<-length(z); if(ll==0) break
    temp<-rbind(temp,as.numeric(c(z[2],z[ll-7],z[ll-4:2],z[ll])))
    ii<-ii+1
  }
  colnames(temp)=c('season','unit','est','se','lower','upper')
  for(tt in 2:temp.pao$nseasons){
    d.temp<-temp[temp$season==tt,]
    nunique<-nrow(d.temp)
    if(nunique>1){
      for(ii in 1:(nunique-1)){
        psi.est<-rbind(psi.est,d.temp[rep(ii,d.temp$unit[ii+1]-d.temp$unit[ii]),3:6])
      }
    }
    psi.est<-rbind(psi.est,d.temp[rep(nunique,temp.pao$nunits-d.temp$unit[nunique]+1),3:6])
  }


  #   temp.psi<-temp.SE<-rep(NA,temp.pao$nunits*(temp.pao$nseasons-1))
  #
  #   for(ii in 1:temp.pao$nunits){
  #     index<-seq(ii,by=temp.pao$nunits,length.out=(temp.pao$nseasons-1))
  #     deriv<-rbind(c(psi_mat[ii,]*psi.est[ii,"est"]*(1-psi.est[ii,"est"]),
  #                    rep(0,npar.gamma),
  #                    rep(0,npar.epsilon)),
  #                  as.matrix(data.frame(array(0,dim=c(temp.pao$nseasons-1,npar.psi)),
  #                                       gamma_mat[index,]*gamma.est[index,"est"]*(1-gamma.est[index,"est"]),
  #                                       array(0,dim=c(temp.pao$nseasons-1,npar.epsilon)))),
  #                  as.matrix(data.frame(array(0,dim=c(temp.pao$nseasons-1,npar.psi)),
  #                                       array(0,dim=c(temp.pao$nseasons-1,npar.gamma)),
  #                                       epsilon_mat[index,]*epsilon.est[index,"est"]*(1-epsilon.est[index,"est"])))
  #     )
  #     temp.VC<-deriv%*%VC[1:(npar.psi+npar.gamma+npar.epsilon),1:(npar.psi+npar.gamma+npar.epsilon)]%*%t(deriv)
  #     psi.temp<-psi.est[ii,"est"]
  #     psi.deriv<-array(0,dim=c(temp.pao$nseasons,2*(temp.pao$nseasons-1)+1))
  #     psi.deriv[1,1]<-1
  #
  #     for(jj in 1:(temp.pao$nseasons-1)){
  #       temp.psi[index[jj]]<-psi.temp*(1-epsilon.est[index[jj],"est"])+(1-psi.temp)*gamma.est[index[jj],"est"]
  #       temp.ind<-c(1,2:(2+jj-1),(temp.pao$nseasons+1):(temp.pao$nseasons+1+jj-1))
  #       psi.deriv[jj+1,temp.ind]<-psi.deriv[jj,temp.ind]*((1-epsilon.est[index[jj],"est"])-gamma.est[index[jj],"est"])
  #       psi.deriv[jj+1,1+jj]<-(1-psi.temp)
  #       psi.deriv[jj+1,temp.pao$nseasons+jj]<- -psi.temp
  #
  #       psi.temp<-temp.psi[index[jj]]
  #     }
  #
  #     temp.psi.VC<-psi.deriv%*%temp.VC%*%t(psi.deriv)
  #     temp.SE[index]<-sqrt(diag(temp.psi.VC)[-1])
  #
  #   }
  #   lower<-plogis(qlogis(temp.psi)-1.96*temp.SE/(temp.psi*(1-temp.psi)))
  #   upper<-plogis(qlogis(temp.psi)+1.96*temp.SE/(temp.psi*(1-temp.psi)))


  # get derived parameters
  #   psi.est<-rbind(psi.est,as.matrix(data.frame(temp.psi,temp.SE,lower,upper)))
  rownames(psi.est)<-paste(rep(temp.pao$unitnames,temp.pao$nseasons),
                           rep(1:temp.pao$nseasons,each=temp.pao$nunits),
                           sep="_")

  return(list(
    beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
              gamma=as.data.frame(gamma.coeff),gamma.VC=gamma.VC,
              epsilon=as.data.frame(epsilon.coeff),epsilon.VC=epsilon.VC,
              p=as.data.frame(p.coeff),p.VC=p.VC,
              VC=VC),
    real=list(psi=as.data.frame(psi.est),
              gamma=as.data.frame(gamma.est),
              epsilon=as.data.frame(epsilon.est),
              p=as.data.frame(p.est))
  ))

}

get_do_4<-function(v,aux,npar){
  ### get coefficients
  instns<-grep("Untransformed",v,fixed=TRUE)
  extract<-v[(instns[1]+3):(instns[1]+2+npar)]
  if(is.null(aux)){ # so importing results from a non R-based run
    # get number of parameters for each parameter type
    instns<-grep("Matrix 1",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.psi<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1

    instns<-grep("Matrix 4",v,fixed=TRUE)
    z<-get_line(v,instns[1])
    npar.p<-as.numeric(unlist(strsplit(z[length(z)],"="))[2])-1

    # make indices for psi and p
    index.psi<-1:npar.psi
    index.p<-(npar.psi+1):(npar)

    # get covariate names
    covnames<-unlist(lapply(1:length(extract),function(zz){
      get_line(extract,zz)[2]
    }))
    covnames.psi<-covnames[index.psi]
    covnames.p<-covnames[index.p]

    # get pao info
    instns<-grep("==>i",v,fixed=TRUE)
    z<-unlist(unlist(strsplit(v[instns[1]],"=")))
    file<-z[length(z)]
    temp.pao<-readPao(file)

    # get fixed parameter info
    instns<-grep("fixed",v,fixed=TRUE)
    if(length(instns)>0){
      fxd<-data.frame()
      for(ii in 1:length(instns)){
        rbind(fxd,get_line(v[instns][ii]))
      }
      fixed<-matrix(fxd[,4],ncol=1)
      rownames(fixed)<-fxd[,2]
    } else{
      fixed=NULL
    }

    # labelling quirks
  } else { # get this info from aux
    # get number of parameters for psi and p
    npar.psi<-ncol(aux$psi_mat)
    npar.p<-ncol(aux$p_mat)

    # make indices for psi and p
    index.psi<-grep(".psi.",extract,fixed=TRUE)
    index.p<-grep(".p.",extract,fixed=TRUE)

    # get covariate names
    covnames.psi<-colnames(aux$psi_mat)
    covnames.p<-colnames(aux$p_mat)

    # get pao info
    temp.pao<-aux$temp.pao

    # get fixed parameter info
    fixed=aux$fixed

  }

  psi.coeff<-get_coeff(v=extract,index=index.psi,n.coeff=npar.psi,
                       cov.names=covnames.psi,b.string="a",b.count=0)
  p.coeff<-get_coeff(v=extract,index=index.psi,n.coeff=npar.p,
                     cov.names=covnames.p,b.string="d",b.count=0)

  ### get VC matrix
  instns<-grep("Variance-Covariance Matrix of",v,fixed=TRUE)
  names<-c(rownames(psi.coeff),rownames(p.coeff))
  VC<-get_VC(v=v,offset=instns+1,n.coeff=npar,cov.names=names)

  psi.VC<-VC[psi.coeff[,"index"],psi.coeff[,"index"]]
  p.VC<-VC[p.coeff[,"index"],p.coeff[,"index"]]

  ## Get real parameter estimates from output
  psi.est<-get_real(real="<psi1>",v=v,con.offset=1,
                    n1=1,nunits=temp.pao$nunits,
                    row.names=temp.pao$unitnames,
                    fixed=fixed,real.offset=0)
  ## check whether <gam1> or <psi2> string used
  string<-ifelse(length(grep("<gam1>",v,fixed=TRUE))>0,"<gam1>","<psi2>")
  gamma.est<-get_real(real=string,v=v,con.offset=1,
                      n1=temp.pao$nseasons-1,nunits=temp.pao$nunits,
                      row.names=paste(temp.pao$unitnames,"_",
                                      rep(1:(temp.pao$nseasons-1),each=temp.pao$nunits),
                                      sep=""),
                      fixed=fixed,real.offset=1)

  epsilon.est<-gamma.est
  epsilon.est[,"est"]<- 1-gamma.est[,"est"]
  epsilon.est[,"lower"]<- 1-gamma.est[,"upper"]
  epsilon.est[,"upper"]<- 1-gamma.est[,"lower"]


  psi.est<-rbind(psi.est,gamma.est)
  row.names(psi.est)<-paste(temp.pao$unitnames,"_",
                            rep(1:temp.pao$nseasons,each=temp.pao$nunits),
                            sep="")


  p.est<-get_real(real="<P[1-1]>",v=v,con.offset=1,
                  n1=temp.pao$nsurveys,nunits=temp.pao$nunits,
                  row.names=paste(rep(temp.pao$unitnames,temp.pao$nsurveys),"_",
                                  rep(temp.pao$surveynames,each=temp.pao$nunits),
                                  sep=""),
                  fixed=fixed,real.offset=temp.pao$nseasons+1)

  return(list(
    beta=list(psi=as.data.frame(psi.coeff),psi.VC=psi.VC,
              p=as.data.frame(p.coeff),p.VC=p.VC,
              VC=VC),
    real=list(psi=as.data.frame(psi.est),
              gamma=as.data.frame(gamma.est),
              epsilon=as.data.frame(epsilon.est),
              p=as.data.frame(p.est))
  ))

}

summedWgt<-function(covnames,param,aic.tab){
  #' Calculate summed model weights for a set of variables for a specific parameter.
  #'
  #' Calculates the summed model weights (and evidence ratios) for a set of variables
  #'   \code{covnames}, for a parameter \code{param} and a set of models stored in an
  #'   \code{aic.tab} object.
  #'
  #'@param covnames a vector of string values containing the names of the variables of interest.
  #'@param param a string with the name of the parameter of interest, e.g., \code{param="psi"}.
  #'@param aic.tab an \code{aic.tab} object that includes the set of models from which the
  #'   summed model weights are to be calculated.
  #'@export
  #'@return returns a data frame with the summed weights and evidence ratios.
  #'
  #'@author Darryl MacKenzie.
  #'
  #'@seealso \code{\link{createAicTable}}
  #'
  #'@keywords htest

  temp<-unlist(lapply(covnames,function(cvn){
    ind<-unlist(lapply(aic.tab$models,function(model){
      grepl(cvn,deparse(get(param,model$model)))
    }))
    return(sum(ind*aic.tab$table$wgt))
  }))
  return(data.frame(covnames=covnames,sum.wgt=temp,ER=temp/(1-temp)))
}

modCombos<-function(param,covs){
  #' Create a list of formulae for all possible first-order combinations for a set of variables.
  #'
  #' Creates a list of formulae for the parameter \code{param} which contains all possible
  #' first-order (i.e., no interactions) combinations for the set of variables \code{cov}.
  #'
  #'@param param a string with the name of the parameter of interest, e.g., \code{param="psi"}.
  #'@param covs a vector of string values containing the names of the variables of interest.
  #'@export
  #'@return  returns a list of formulae of each combination of variables.
  #'
  #'@author Darryl MacKenzie.
  #'
  #'@seealso \code{\link{createAicTable}}, \code{\link{predict.occMod}}
  #'
  #'@keywords models
  #'
  combs <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), length(covs)))
  vars <- apply(combs, 1, function(i) covs[i])
  vars[[1]]<-1
  form <- paste(param,"~", lapply(vars, paste, collapse=" + "), sep = "")
  form <- lapply(form, as.formula)
  return(form)
}

### extract fitted real values from occMod object
### generic function
fitted.occMod<-function(object,param=NULL){
  #' Extract fitted values from occupancy model, i.e., real parameter estimates
  #'
  #'@param object an object of class \code{occMod} which is the result from fitting a model using the \code{\link{occMod}} function.
  #'@param param name of real parameter to extract
  #'@export
  #'@return a data.frame of the estimates, standard errors and confidence intervals limits
  #'
  #'@author Darryl MacKenzie

  if(is.null(param)){
    stop(cat("Error. Need to specify parameter type. Available parameters are:\n",names(object$real)))
  }
  if(param%in%names(object$real)){
    temp<-get(param,object$real)
    return(temp)
  } else {
    stop(cat("Error. Need to specify parameter type. Available parameters are:\n",names(object$real)))
  }
}

### generic function
coef.occMod<-function(object,param=NULL){
  #' Extract regression coefficients from occupancy model, i.e., real parameter estimates
  #'
  #'@param object an object of class \code{occMod} which is the result from fitting a model using the \code{\link{occMod}} function.
  #'@param param name of model component to extract regression coefficients for
  #'@export
  #'@return a data.frame of the estimates and standard errors
  #'
  #'@author Darryl MacKenzie

  if(is.null(param)){
    stop(cat("Error. Need to specify parameter type. Available parameters are:\n",names(object$beta)))
  }
  if(param%in%names(object$beta)){
    temp<-get(param,object$beta)[,-3]
    return(temp)
  } else {
    stop(cat("Error. Need to specify parameter type. Available parameters are:\n",names(object$beta)))
  }
}

calc_real<-function(cov,coeff,VC,conf=0.95,link="logit",rownames,fixed=NULL){

  est<-cov%*%coeff
  var<-apply(cov,1,function(xx) xx%*%VC%*%xx)
  alpha<-(1-conf)/2;  z=qnorm(alpha); lower=est+z*sqrt(var); upper=est-z*sqrt(var)
  colnames(lower)<-paste("lower",conf,sep="_"); colnames(upper)<-paste("upper",conf,sep="_")

  real<-data.frame(est=est,se=sqrt(var),lower,upper)

  if(link=="logit") {
    real$est<-plogis(real$est)
    real$se<-real$se*(real$est*(1-real$est))
    real[,-c(1,2)]<-plogis(as.matrix(real[,-c(1,2)]))
  }
  rownames(real)<-rownames

  if (!is.null(fixed)) {
    s=gsub('_.+','',rownames)
    for (i in 1:nrow(fixed)) {
      j=which(fixed$param[i]==s)
      if (length(j)>0) { real[j,]=NA; real[j,1]=fixed$value[i] }
    }
  }
  return(real)
}

calc_real_ms<-function(cov,coeff,VC,conf=0.95,rownames,fixed=NULL){

  est<-cov%*%coeff
  var<-apply(cov,1,function(xx) xx%*%VC%*%xx)
  alpha<-(1-conf)/2;  z=qnorm(alpha); lower=est+z*sqrt(var); upper=est-z*sqrt(var)
  colnames(lower)<-paste("lower",conf,sep="_"); colnames(upper)<-paste("upper",conf,sep="_")

  real<-data.frame(est=est,se=sqrt(var),lower,upper)

  if(link=="logit") {
    real$est<-plogis(real$est)
    real$se<-real$se*(real$est*(1-real$est))
    real[,-c(1,2)]<-plogis(as.matrix(real[,-c(1,2)]))
  }
  rownames(real)<-rownames

  if (!is.null(fixed)) {
    s=gsub('_.+','',rownames)
    for (i in 1:nrow(fixed)) {
      j=which(fixed$param[i]==s)
      if (length(j)>0) { real[j,]=NA; real[j,1]=fixed$value[i] }
    }
  }
  return(real)
}
calc_psi_c<-function(psi.cov,psi.coeff,p.cov,p.coeff,VC,det.data,rownames,conf=0.95,link="logit"){
  k<-ncol(det.data);   s<-nrow(det.data)
  psi_cond <- function(ii) {    #  function to compute psi-c for 1 site...
    if(sum(det.data[ii,],na.rm=TRUE)>0){     #  if at least 1 detection, retrun psi-c = 1.0
      return(c(1,0,rep(1,2*length(conf))))
    } else{                                         #  otherwise, compute psi-c via delta method...
      ## create matrix
      temp.psi<-c(psi.cov[ii,],rep(0,length(p.coeff)))
      temp.p<-cbind(array(0,dim=c(k,length(psi.coeff))),p.cov[ii+((0:(k-1))*s),])
      temp<-rbind(temp.psi,temp.p)
      probs<-plogis(temp%*%c(psi.coeff,p.coeff))
      is.na(probs)<-c(FALSE,is.na(det.data[ii,])); i=(!is.na(probs)) # i=surveys w/o missing data
      probs.VC<-(temp%*%VC%*%t(temp))[i,i]                           #  only include p's with data
      grad<-as.numeric((probs*(1-probs)))[i]
      probs.VC<-diag(grad)%*%probs.VC%*%diag(grad)
      p_star<-1-prod(1-probs[-1],na.rm=TRUE)
      est<-probs[1]*(1-p_star)/(1-probs[1]*p_star)
      grad[1]<-(1-p_star)/(1-probs[1]*p_star)^2;  i[1]=FALSE  #  i[1]=F so probs[i] only has p's w/ data
      grad[-1]<--(probs[1]*(1-probs[1])*(1-p_star))/((1-probs[1]*p_star)^2*(1-probs[i]))
      var<-grad%*%probs.VC%*%grad
      if (sum(dim(var))==0) var=0

      ## backtransform to logit to calc CIs
      logit.est<-qlogis(est);  logit.se<-sqrt(var)/(est*(1-est))
      alpha<-(1-conf)/2; z<--qnorm(alpha)
      lower<-plogis(logit.est-z*logit.se);  upper<-plogis(logit.est+z*logit.se)

      return(c(est,sqrt(var),lower,upper))
    }
  }
  psi_c=NULL; for (i in 1:s) psi_c=cbind(psi_c,psi_cond(i))

  psi_c<-t(psi_c);   rownames(psi_c)<-rownames
  colnames(psi_c)<-c("est","se",paste("lower",conf,sep="_"),paste("upper",conf,sep="_"))
  return(psi_c)
}
diagbind=function(a,b) {  #  diagonal bind of 2 matrices
  if (!is.matrix(a)) a=matrix(a,nrow=1)
  if (!is.matrix(b)) b=matrix(b,nrow=1)
  x=matrix(0,nrow(a)+nrow(b),ncol(a)+ncol(b))
  x[1:nrow(a),1:ncol(a)]=a
  x[nrow(a)+1:nrow(b),ncol(a)+1:ncol(b)]=b
  rownames(x)=c(rownames(a),rownames(b))
  return(x)
}

print_one_site_estimates <- function(mod,site=1) {
  #'  Print real param estimates for i-th site only
  #' \code{print_one_site_estimates} Prints real param estimates for ith site only
  #'@param mod model output from \code{\link{occMod}} function.
  #'@param site site number(s) to print (default=1)
  #'@export
  #'@return nothing
  w=NULL; s=names(mod$real); cat(mod$modname,'\n')
  for (i in 1:length(mod$real)) {
    x=as.matrix(mod$real[[i]][site,],ncol=4) #  for ith site only
    rownames(x)=paste(s[i],rownames(x),sep="_"); t=gsub('.est','',colnames(x))
    if (ncol(x)==4) w=rbind(w,x) else {
      x1=matrix(unlist(x[,1:4]),ncol=4); rownames(x1)=paste(t[1],site,sep="_")
      x2=matrix(unlist(x[,5:8]),ncol=4); rownames(x2)=paste(t[5],site,sep="_")
      w=rbind(w,x1,x2)
    }
  }
  print(w,quote=F)
}

simplify_dm <- function(dm,pao) {
  #'  simplify design matrix by replacing all individual covariates which are all one or zero
  #'  with '1' or '0' in the design matrix.  This makes the design matrix much easier to read
  #'  and faster for the model.
  #' \code{simplify_dm} Replaces ind. covariates which are all 0 or 1, with 0 or 1.
  #'@param dm model design matrix.
  #'@param pao pao file created by occMod.???
#  #'@export
  #'@return simplified design matrix
  for (i in 1:nrow(dm))
    for (j in 1:ncol(dm)) {
      s=dm[i,j];
      k=which(colnames(pao$unitcov)==s);
      if (length(k)>0) {
        t=pao$unitcov[,k]
        if (sum(t!=1,na.rm=TRUE)==0) dm[i,j]="1"
        if (sum(t!=0,na.rm=TRUE)==0) dm[i,j]="0"
      }
      k2=which(names(pao$survcov)==s)
      if (length(k2)>0) {
        t=pao$survcov[,k2]
        if (sum(t!=1,na.rm=TRUE)==0) dm[i,j]="1"
        if (sum(t!=0,na.rm=TRUE)==0) dm[i,j]="0"
      }
    }
  k=which(colSums(dm!="0")==0);  # del cols which are all zero
  if (length(k)>0) dm=dm[,-k]
  return(dm)
}

new_simplify_dm <- function(dm,pao) {
  #'  simplify design matrix by replacing all individual covariates which are all one or zero
  #'  with '1' or '0' in the design matrix.  This makes the design matrix much easier to read
  #'  and faster for the model.
  #' \code{simplify_dm} Replaces ind. covariates which are all 0 or 1, with 0 or 1.
  #'@param dm model design matrix.
  #'@param pao pao file created by occMod.???
  #  #'@export
  #'@return simplified design matrix
  for (i in 1:nrow(dm))
    for (j in 1:ncol(dm)) {
      s=dm[i,j];
      k=which(colnames(pao$unitcov)==s);
      if (length(k)>0) {
        t=pao$unitcov[,k]
        if (sum(t!=1,na.rm=TRUE)==0) dm[i,j]="1"
        if (sum(t!=0,na.rm=TRUE)==0) dm[i,j]="0"
      }
      k2=which(names(pao$survcov)==s)
      if (length(k2)>0) {
        t=pao$survcov[,k2]
        if (sum(t!=1,na.rm=TRUE)==0) dm[i,j]="1"
        if (sum(t!=0,na.rm=TRUE)==0) dm[i,j]="0"
      }
    }
  k=which(colSums(dm!="0")>0);  # del cols which are all zero
  if (length(k)>0) dm=dm[,k]
  return(list(dm=dm,nzcol=k))
}

get_not01_dm <- function(dm, covlst) {
  #'  get site/survey covarite names which are not all 0 or 1 in the design matrix
  #'  so the 0/1 covariate names (eg., Intercept) can be omitted from the temp.pao file.
  #'
  #' \code{get_not01_dm} Gets site/survey covariate names which are not all 0 or 1 in design matrix.
  #'@param dm model design matrix.
  #'@param covlst character vector of site/survey covariate names
#  #'@export
  #'@return character vector of site/survey covariate names
  lst=rep(0,length(covlst))
  for (dmi in dm) { k=which(dmi == covlst); if (length(k)>0) lst[k]=1 }
  return(lst)
}
psi_cond <- function(m0,paodata)  {
  #'  get conditional occupancy estimates from single-season ("so") model
  #'
  #' \code{psicond} Computes conditional (on detection history) occupancy estimates
  #'   from single-season ("so") model.  Uses "delta" method, with numerical gradient.
  #'@param m0 model object (result from running "occMod").
  #'@param paodata pao data object
#  #'@export
  #'@return numeric matrix with column of estimates and column of std. errors
  #'
  if (! "so" %in% class(m0)) cat('\n psicond only works for type, "so"\n') else {
    v=calc_psi_c(m0$dmat$psi, m0$beta$psi$est, m0$dmat$p, m0$beta$p$est, m0$beta$VC,
                 det.data=m0$data$det.data,rownames=m0$data$unitnames)
  }
  return(v)
}
genpresEV <- function(N=100,K=4,psi=.75,p=.4,fp=.0,b=.0,eps=.0,gam=.0,sps=4) {
  #' generate expected value detection histories
  #'
  #' \code{genpresEV} generates expected value detection history frequencies for
  #' single or multi-season models (false-positives optional)
  #'@param N number of sites
  #'@param K number of surveys
  #'@param psi probability of occupancy
  #'@param p probability of detection (not false positive)
  #'@param fp probability of false detection
  #'@param b probability of "sure" detection
  #'@param eps probability of local extinction between seasons
  #'@param gam probability of local colonization between seasons
  #'@param sps number of surveys per season
  #'@export
  #'@return list(hst=det. history matrix, frq=vector of frequencies)
  #'@note This function generates at least 2^K detection-histories, so if K>8, it may take a
  #'       very long time to run.
  #'
  #'@examples
  #'    # generate expected value data for standard single-season model (no false positives)
  #' hstfrq=genpresEV(N=1000, K=4, psi=.8, p=.3)  # generate det. hists and freqs
  #' EVpao=createPao(hstfrq$hst,frq=hstfrq$frq)                  #  create pao file
  #'
  #'    # generate expected value data for single-season model with false positive detections
  #' hstfrq=genpresEV(N=1000, K=4, psi=.8, p=.3, fp=0.1, b=0.5)  # generate det. hists and freqs
  #' EVpao=createPao(hstfrq$hst,frq=hstfrq$frq)                  #  create pao file
  #'
  #'    # generate expected value data for multi-season model (no false positive detections)
  #' hstfrq=genpresEV(N=1000, K=4, psi=.8, p=.3, eps=.2, gam=.1, sps=2)  # generate det. hists and freqs
  #' EVpao=createPao(hstfrq$hst,frq=hstfrq$frq,nsurveyseason=2)   #  create pao file
  #'
  p=c(fp,p); b=c(0,b); phi=c(gam,1-eps); newseasn=seq(sps,K,sps)
  srv <- function(k,n,hst,isocc) {
    if (n>0) {
      occ=n*phi[isocc+1]
      cap(k,occ,  hst,1) # generate det. histories for occupied sites
      cap(k,n-occ,hst,0) # generate det. histories for unoccupied sites
    }
  }
  cap <- function(j,n,hst,isocc) {  #  determine how many sites with detections or
    if (n>0) {                      #  non-detections (recursively) each survey
      ndet=n*p[isocc+1]; nsure=ndet*b[isocc+1]; unsure=ndet-nsure; notdet=n-ndet;
      hsure=hunsure=hnotdet=hst;  hunsure[j]=1; hsure[j]=2
      if (j<K) {
        if (j %in% newseasn) {
          srv(j+1,nsure,hsure,isocc); srv(j+1,unsure,hunsure,isocc); srv(j+1,notdet,hnotdet,isocc)
        } else {
          cap(j+1,nsure,hsure,isocc); cap(j+1,unsure,hunsure,isocc); cap(j+1,notdet,hnotdet,isocc)
        }
      } else {
        s1=paste(hunsure,collapse=''); s2=paste(hsure,collapse=''); s0=paste(hnotdet,collapse='')
        xhst<<-c(xhst,s0,s1,s2); frq<<-c(frq,notdet,unsure,nsure)
      }
    }
  }
  occ=N*psi; notocc=N-occ; frq=xhst=NULL;
  cap(1,occ,rep(0,K),1);   # generate det. histories for occupied sites
  cap(1,notocc,rep(0,K),0) # generate det. histories for unoccupied sites
  s=unique(xhst);
  sumfrq=sapply(1:length(s),function(i) sum(frq[xhst==s[i]])) # summarize freqs by det. hist.
  i=which(sumfrq>0);   #  omit det. hists with frq=0
  hst=matrix(as.numeric(unlist(strsplit(s[i],''))),ncol=K,byrow=T);
  colnames(hst)=paste(1:K)
  return(list(hst=hst, frq=sumfrq[i]))
}

genpresEVms <- function(N=100,K=4,psi=.5,R=.8,p1=.4,p2=.3,dlta=.6,nstates=3,
                         psi1=0,Cpsi=0,CR=0,phi=0,p=0,sps=4) {
  #' generate expected value detection histories for multi-state models
  #'
  #' \code{genpresEVms} generates expected value detection history frequencies for
  #' single or multi-season multi-state models.
  #'@param N number of sites
  #'@param K number of surveys
  #'@param psi probability of occupancy
  #'@param R probability of breeding (state=2), given occupancy
  #'@param p1 probability of detection, given occupied (state=1)
  #'@param p2 probability of detection, given occupied and breeding (state=2)
  #'@param dlta probability of detecting breeding, given occupancy state=2
  #'@param nstates number of states (including "not occupied state")
  #'@param psi1 probability vector of initial states (not occ, occ, occ_breeding)
  #'@param Cpsi conditional probability of occupancy, given previous state
  #'@param CR conditional probability of breeding, given occupancy and previous state
  #'@param phi transition probability matrix (an alternative to specifying psi,R,Cpsi,CR)
  #'@param p detection probability matrix (an alternative to specifying p1,p2,dlta)
  #'@param sps number of surveys per season
  #'@export
  #'@return list(hst=det. history matrix, frq=vector of frequencies)
  #'@note This function generates at least nstates^(K-1) detection-histories, so if K>6, it may take a
  #'       very long time to run.
  #'
  #'@examples
  #'    # generate expected value data for standard single-season model (psi,R parameterization)
  #' x1=genpresEVms(N=100,K=4,psi=.75,R=.8,p1=.4,p2=.3,dlta=.6,nstates=3,psi1=0,Cpsi=0,CR=0,phi=0,p=0,sps=4)
  #'
  #'    # generate expected value data for multi-season model (psi,R parameterization)
  #' x2=genpresEVms(N=100,K=4,psi=.75,R=.8,p1=.4,p2=.3,dlta=.6,nstates=3,psi1=0,
  #'                Cpsi=c(.66,.1,.15),CR=c(.72,.77,0),phi=0,p=0,sps=2)
  #'
  #'    # generate expected value data for multi-season model (phi parameterization)
  #' nstates=3; p=matrix(0,nstates,nstates); p[2,2]=.4; p[3,2]=.3*.4; p[3,3]=.3*.6
  #' psi1=c(.25,.75*(1-.8),.75*.8)
  #' phi=rbind(c(.66*(1-.72),.66*.72),c(.1*(1-.77),.1*.77),c(.15,.0)); phi=cbind(1-rowSums(phi),phi)
  #' x3=genpresEVms(N=100,K=4,nstates=3,psi1=psi1,phi=phi,p=p,sps=2)
  #'
  #'    # generate expected value data for multi-season model (phi parameterization)
  #' nstates=3; p=matrix(0,nstates,nstates); p[2,2]=p[3,2]=p[3,3]=.25
  #' psi1=c(.167,.5,.333)
  #' phi=rbind(c(.5,.333),c(.5,.333),c(.5,.333)); phi=cbind(1-rowSums(phi),phi)
  #' x4=genpresEVms(N=100,K=4,nstates=3,psi1=psi1,phi=phi,p=p,sps=2)
  #'
  #'    # generate expected value data for multi-season model (phi parameterization)
  #' nstates=4; p=matrix(0,nstates,nstates); p[2:4,2]=p[3:4,3]=p[4,4]=.2
  #' psi1=c(0,.4,.3,.2); psi1[1]=1-sum(psi1)
  #' phi=matrix(0,nstates,nstates); phi[,2]=.4; phi[,3]=.3; phi[,4]=.2; phi[,1]=1-rowSums(phi)
  #' x4=genpresEVms(N=100,K=4,nstates=4,psi1=psi1,phi=phi,p=p,sps=2)
  #' EVpao=createPao(x4$hst,frq=x4$frq,nsurveyseason=2,nstates=4)   #  create pao file
  #'
  if (length(psi1)==1) {
    psi1=c(1-psi, psi*(1-R), psi*R)
    p=matrix(c(1,0,0, 1-p1,p1,0, 1-p2,p2*(1-dlta),p2*dlta),3,3,byrow=T)
    phi=cbind(Cpsi*(1-CR),Cpsi*CR); phi=cbind(1-rowSums(phi),phi)
  }
  seasn_end=seq(sps,K,sps);  p[,1]=1-rowSums(p[,-1]); p[1,1]=1

  srv <- function(j,n,hst,istate) {
    if (n>0) for (i in 1:nstates) cap(j+1,n*phi[istate,i],hst,i)
  }
  cap <- function(j,n,hst,istate) {  #  determine how many sites with detections or
    if (n>0) {                      #  non-detections (recursively) each survey
      ndet=n*p[istate,]; h=matrix(rep(hst,nstates),nstates,byrow=T); i=2:nstates; h[i,j]=i-1
      if (j<K) {
        for (i in 1:nstates)
          if (j %in% seasn_end) srv(j,ndet[i],h[i,],istate) else cap(j+1,ndet[i],h[i,],istate)
      } else {
        for (i in 1:nstates)
          if (ndet[i]>0) { xhst<<-c(xhst,paste(h[i,],collapse='')); frq<<-c(frq,ndet[i])}
      }
    }
  }
  occ=N*psi1; frq=xhst=NULL;
  for (i in 1:nstates) cap(1,occ[i],rep(0,K),i);   # generate det. histories for occupied sites
  s=unique(xhst);
  sumfrq=sapply(1:length(s),function(i) sum(frq[xhst==s[i]])) # summarize freqs by det. hist.
  i=which(sumfrq>0);   #  omit det. hists with frq=0
  hst=matrix(as.numeric(unlist(strsplit(s[i],''))),ncol=K,byrow=T);
  colnames(hst)=paste(1:K)
  return(list(hst=hst, frq=sumfrq[i]))
}

clean_covs_from_pao <- function(psi.dm, temp.pao) {
  #  get rid of site covariates not in design matrices
  s1=as.character(colnames(temp.pao$unitcov))
  s2=unique(as.character(psi.dm))
  s3=NULL; for (i in 1:length(s1)) if (s1[i] %in% s2) s3=c(s3,s1[i])
  if (!is.null(s3)) {
    temp.pao$unitcov=temp.pao$unitcov[,s3,drop=F]; temp.pao$nunitcov=ncol(temp.pao$unitcov)
  } else temp.pao$unitcov=temp.pao$nunitcov=0

  #  get rid of survey covariates not in design matrix
  s1=as.character(colnames(temp.pao$survcov))
  s2=unique(as.character(psi.dm))
  s3=NULL; for (i in 1:length(s1)) if (s1[i] %in% s2) s3=c(s3,s1[i])
  if (!is.null(s3)) { temp.pao$survcov=temp.pao$survcov[,s3,drop=F]; temp.pao$nsurvcov=ncol(temp.pao$survcov)
  } else temp.pao$nsurvcov=temp.pao$survcov=0
  return(temp.pao)
}
