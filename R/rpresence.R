#'
#' R Interface for Program PRESENCE
#'
#' Provide an R interface for running some of the occupancy models available in
#' Program PRESENCE (\url{http://www.mbr-pwrc.usgs.gov/software/presence.html}),
#' plus some additional helpful routines.
#'
#' details:
#' \tabular{ll}{
#' Package: \tab RPresence\cr
#' Type: \tab Package\cr
#' Date: \tab 2017-11-27\cr
#' License: \tab GPL-2\cr
#' }
#'
#'@author Darryl MacKenzie and Jim Hines
#   importFrom Rcpp evalCpp
#'@useDynLib RPresence
#'
#' @examples
#' \dontrun{
#' # load a csv file with detection-histories
#' filename<-system.file("extdata/Blue_Ridge_pg99.csv",package="RPresence")
#' salmdr.csv<-read.csv(filename,header=FALSE)
#'
#' # Create PRESENCE input file object from csv
#' salmdr.data<-createPao(salmdr.csv,paoname="salmdr.pao")
#'
#' ## fit some models
#' mod1<-occMod(model=list(psi~1,p~SURVEY),data=salmdr.data,type="so")
#' mod2<-occMod(model=list(psi~1,p~1)     ,data=salmdr.data,type="so")
#'
#' ## create AIC table
#' models<-list(mod1,mod2)
#' results<-createAicTable(models)
#' summary(results)
#'
#' ## print real estimates (for 1st site only)
#' cat('========== Model: psi(.)p(SURVEY) -\n')
#' cat('   psi:\n'); print(fitted(mod1,'psi')[1,]) # print 1st site psi estimates
#' i=which(!duplicated(fitted(mod1,'p')$est))      #  get survey-specific p estimates
#' cat('   p  :\n'); print(fitted(mod1,'p')[i,])   # print 1st site p estimates
#' cat('========== Model: psi(.)p(.)      -\n')
#' cat('   psi:\n'); print(fitted(mod2,'psi')[1,]) # print 1st site psi estimates
#' cat('   p  :\n'); print(fitted(mod2,'p')[1,])   # print 1st site p estimates
#' }
# Maintainer: Darryl MacKenzie <darryl@proteus.co.nz>
"_PACKAGE"
