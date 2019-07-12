#' HouseFinch data csv file (multi-season example with covariates)
#'
#' To investigate the westward expansion during the period 1976–2001, data
#' from 694 BBS routes within 2600 km from the Long Island point of release were
#' considered at 5-year intervals (i.e., 1976, 1981, 1986, . . . ).
#'
#' \itemize{
#' \item 694 sites
#' \item 300 surveys
#' \item 6 seasons (50 surveys per season)
#' \item covariate: d = distance from release point (Long Island) in 100k increments.
#' \item covariate: f = observed frequency of occurrence indicating whether house finches were detected on >10 stops
#' in the route in any previous year or not (i.e., were locally highly abundant).
#' }
#' @docType data
#' @keywords datasets
#' @name HouseFinch.csv
### @usage read.csv(system.file('extdata/housefinch.csv',package="RPresence"),stringsAsFactors=FALSE,header=FALSE)
#' @format csv file with 694 rows and 300 columns
#' @examples
#' #  test multi-season model using "House Finch" data (pg 363 in new occupancy book)
#' fname=system.file('extdata/housefinch.csv',package="RPresence")  #  get filename
#' #  read detection histories and replace "-" with NA
#' det_hist=read.csv(fname,stringsAsFactors=FALSE,header=FALSE); det_hist[det_hist=='-']=NA
#' #  read site covariate - "d"
#' d=read.csv(gsub('.csv','_sitecov.csv',fname),stringsAsFactors=FALSE,header=FALSE)
#' #  read survey covariate - "f"
#' f=unlist(read.csv(gsub('.csv','_survcov.csv',fname),stringsAsFactors=FALSE,header=FALSE))
#' f=data.frame(f,stringsAsFactors=FALSE); colnames(d)='d';
#' #  create pao object
#' pao=createPao(det_hist,unitcov=d,survcov=f,nsurveyseason=rep(50,6))
#' #  run model with psi(distance band),gam(year+dist),eps(dist),p(year*dist+frq-of-occurrence) (run-time=15.6 min)
#' hf1<-occMod(model=list(psi~d,gamma~SEASON+d,epsilon~d,p~SEASON*d+f),data=pao,type="do.1",outfile="modname",VCoutput="nose");
#' #  run model with psi(distance band),gam(year*dist),eps(dist),p(year*dist+frq-of-occurrence)
#' hf2<-occMod(model=list(psi~d,gamma~SEASON*d,epsilon~d,p~SEASON*d+f),data=pao,type="do.1",outfile="modname",VCoutput="nose");
#' ## create AIC table
#' results<-createAicTable(list(hf1,hf2))
#' cat('House Finch example\n'); print(summary(results))
