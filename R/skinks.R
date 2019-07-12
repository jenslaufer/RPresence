#' Grand Skinks data -multi-season example with a covariate
#'
#' Data have been collected near Macraes Flat, Otago, NZ by the Department of Conservation.
#' Rock outcrops were surveyed by DOC up to 3 times per year for 5 years.
#' See MacKenzie et al. (2006) for more details.
#'
#' \itemize{
#' \item 352 sites
#' \item 15 surveys
#' \item 5 seasons
#' }
#'
#' @docType data
#' @keywords datasets
#' @name skinks.pao
#' @usage mypao=readPao(system.file("extdata/skinks.pao",package="RPresence"))
#' @format pao object with detection history data, and a covariate
#' @examples
#' #  test multi-season model using "Skinks" data
#' #   read pao data object (assuming already created using PRESENCE)
#' data<-readPao(system.file('extdata/skinks.pao',package="RPresence"))
#' #   run a model...
#' m0<-occMod(model=list(psi~Pasture,gamma~1,epsilon~1,p~1),data=data,type="do.1")
#' cat("skinks AIC:",m0$aic,"\n")
#' print(summary(m0))
