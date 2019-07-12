#' Occupancy detection data for California Spotted Owls (single-season, multi-state analysis)
#'
#' One of the motivating examples for Nichols et al. (2007a) was estimation of the
#' reproductive rate of California spotted owls (Strix occidentalis) in the central
#' Sierra Nevada (California, USA), while allowing for state classification uncertainty
#' (e.g., whether territories are unoccupied, occupied without breeding, or
#'   occupied with breeding) and variable sampling protocols.
#'
#' Survey data collected during April to mid-August 2004 from 54 potential
#' territories in the Eldorado study area were analyzed using a multi-state occupancy
#' modeling framework.
#'
#' \itemize{
#' \item 54 sites
#' \item 5 surveys
#' \item 1 season
#' \item 3 occupancy states ("0"=unoccupied, "1"=occupied, "2"=occupied with breeding)
#' \item contains missing data (denoted by ".")
#' }
#' @docType data
#' @keywords datasets
#' @name cal_owl_multistate_data.csv
## @usage read.csv(system.file("extdata/cal_owl_multistate_data.csv",package="RPresence"),as.is=TRUE)
#' @format csv file with 54 rows and 5 columns
#' @examples
#' #  Single-season multi-state example (2 occupancy states, psi-R parameteization) (pg 231 in book)
#'
#' #    read detection history data from csv file...
#' csv<-read.csv(system.file("extdata/cal_owl_multistate_data.csv",package="RPresence"),as.is=T)
#' csv[csv=="."]=NA
#'
#' sitenames=csv[,1]  #  sitenames in 1st column
#' dethist=csv[,-1];  #  get rid of 1st column (site name)
#' nsites=nrow(dethist) #  set number of sites,surveys from det. history data
#' nsrvys=ncol(dethist)
#'
#' #  create survey covariate to categorize surveys 1-2 as 1 period
#' #                                    and surveys 3-5 as another period
#' #    Since it's a survey covariate, it will be a NxT matrix (N=nsites, T=nsurveys)
#' #    The covariate matrix would be:
#' #         1 2 3 4 5
#' #  site1 [1 1 2 2 2]
#' #  site2 [1 1 2 2 2]
#' #   :    [: : : : :]
#' #  siteN [1 1 2 2 2]
#' #    Filling in by cols means we repeat "1" N times, then "1" N times,
#' #                  then "2" N times, then "2" N times, then "2" N times.
#' #    We save it as a data frame (without the matrix dimensions).
#' cov2=data.frame(PER=as.factor(c(rep(1,2*nsites),rep(2,3*nsites))))
#'
#' #          create input "pao" object, for use with occMod function
#' data=createPao(dethist,survcov=cov2,title="Cal Owl example")
#'
#' xmods=list(); i=1       #  run each model and save in list variable, "xmods"
#' xmods[[i]]=occMod(model=list(psi~1,r~1,p~1,     delta~1),data=data,type="do.ms.2");i=i+1
#' xmods[[i]]=occMod(model=list(psi~1,r~1,p~SURVEY,delta~1),data=data,type="do.ms.2");i=i+1
#' xmods[[i]]=occMod(model=list(psi~1,r~1,p~STATE, delta~1),data=data,type="do.ms.2");i=i+1
#' xmods[[i]]=occMod(model=list(psi~1,r~1,p~1,     delta~PER),data=data,type="do.ms.2");i=i+1
#' xmods[[i]]=occMod(model=list(psi~1,r~1,p~SURVEY,delta~PER),data=data,type="do.ms.2");i=i+1
#' xmods[[i]]=occMod(model=list(psi~1,r~1,p~STATE, delta~PER),data=data,type="do.ms.2");i=i+1
#'
#' #     create AIC table of model results and print
#' results2<-createAicTable(xmods); cat('Cal Owl example\n'); print(results2$table)
#'
#' #     print table 5.1 from book...
#' cat('CA spotted owl reproduction Table 5.1 (pg 234 in book)\n')
#' estmts=xmods[[6]]$real
#' estimate_table=rbind(estmts$psi[1,],estmts$r[1,],estmts$p1[1,],estmts$p2[1,],estmts$delta[1,],estmts$delta[109,])
#' rownames(estimate_table)=c('psi','R','p1','p2','delta1','delta2'); print(estimate_table)
NULL