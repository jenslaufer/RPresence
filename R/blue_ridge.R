#' Blue-Ridge salamander data (single-season example with covariate)
#'
#' The Blue-Ridge two-lined salamander (Eurycea wilderae) is one of over 30
#' species of salamanders that occur within the Great Smoky Mountains National
#' Park (Dodd, 2003).  These data are from 2001 for one species
#' (Blue-Ridge two-lined salamander.  See MacKenzie et. al (2017) for more info.
#'
#' \itemize{
#' \item 39 sites
#' \item 5 surveys
#' \item 1 season
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Blue_Ridge_pg135.csv
## @usage read.csv(system.file("extdata/Blue_Ridge_pg135.csv",package="RPresence"),header=FALSE)
#' @format csv file with 39 rows and 5 columns
#' @examples
#' #   Exanple 1 - Single-season model with Salamander data (pg 135 in new occ. book)
#' # load a csv file with detection-histories
#' filename<-system.file("extdata/Blue_Ridge_pg135.csv",package="RPresence")
#' salmdr.csv<-read.csv(system.file("extdata/Blue_Ridge_pg135.csv",package="RPresence"),header=FALSE)
#'
#' # Create PRESENCE input file object from csv
#' salmdr.data<-createPao(salmdr.csv)
#'
#' ## fit some models
#' salmod1<-occMod(model=list(psi~1,p~SURVEY),data=salmdr.data,type="so")
#' salmod2<-occMod(model=list(psi~1,p~1)     ,data=salmdr.data,type="so")
#'
#' ## create AIC table
#' models<-list(salmod1,salmod2)
#' results<-createAicTable(models)
#' cat('Blue-Ridge example (table 4.2 in book)\n'); print(summary(results))
NULL