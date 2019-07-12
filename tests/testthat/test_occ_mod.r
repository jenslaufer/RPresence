rm(list=ls()); library(RPresence); setwd('h:/a/aatmp'); steptime=NULL
#===========================================================
t=Sys.time()
s='"c:/hines/y/presence2/genpres/genpres8.exe" N=100 T=4 psi=.75 p=0.60 P0MISS nrand=0 QUIET'
i=system(s)
a=readPao('genpres.pao')
a$surveynames=c(paste0("1-",a$surveynames),paste0("2-",a$surveynames))
m1=occMod(data=a, type='so', model=list(psi~1,p~1),outfile='savepao.out',quiet=FALSE)
fixed=data.frame(param=paste0('p',1:4),value=rep(1,4))
m2=occMod(data=a, type='so', model=list(psi~1,p~1),fixed=fixed,outfile='savepao.out',quiet=FALSE)
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
s='"c:/hines/y/presence2/genpres/genpres8.exe" N=100 T=4 psi=.75 theta=.5 p1=0.50 p2=0.60 eps=.2 gam=.1 P0MISS SRVYPERSEASN=2,2 METH=2 nrand=0 QUIET'
i=system(s)
a=readPao('genpres.pao')
a$surveynames=c(paste0("1-",a$surveynames),paste0("2-",a$surveynames))
m1=occMod(data=a, type='do.mm', model=list(psi~1,theta~1,epsilon~1,gamma~1,p~-1+DEVICE),outfile='savepao.out',quiet=FALSE,chsq=TRUE)
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
s='"c:/hines/y/presence2/genpres/genpres8.exe" N=100 T=4 psi=.75 theta=.5 p1=0.50 p2=0.60 P0MISS METH=2 nrand=0 QUIET'
i=system(s)
a=readPao('genpres.pao')
a$surveynames=c(paste0("1-",a$surveynames),paste0("2-",a$surveynames))
m1=occMod(data=a, type='so.mm', model=list(psi~1,theta~1,p~-1+DEVICE),outfile='savepao.out',quiet=FALSE)
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#  test corr-detections (single-season) model using expected-value data from genpres...
setwd('h:/a/aatmp'); library(RPresence)
N=100; TT=6; psi=.75; p=.4; th0=.3; th1=.6; th0pi=0
args=paste0("N=",N," T=",TT," psi=",psi," p=",p," th0=",th0," th1=",th1," th0pi=",th0pi)
cmd=paste('"c:/hines/y/presence2/genpres/genpres8.exe"',args,'QUIET P0MISS nrand=0')
cat(cmd,'\n'); system(cmd)
CDpao=readPao('genpres.pao'); CDpao$nmethods=2
CDmod1<-occMod(model=list(psi~1,theta~PRIME,p~1,th0pi~1),data=CDpao,type="so.cd",outfile='modname')  # generate results
CDmod2<-occMod(model=list(psi~1,theta~SURVEY,p~1,th0pi~1),data=CDpao,type="so.cd",outfile='modname')  # generate results
CDmod3<-occMod(model=list(psi~1,theta~PRIME+SURVEY,p~1,th0pi~1),data=CDpao,type="so.cd",outfile='modname')  # generate results
CDmod4<-occMod(model=list(psi~1,theta~PRIME*SURVEY,p~1,th0pi~1),data=CDpao,type="so.cd",outfile='modname')  # generate results
CDresults=createAicTable(list(CDmod1,CDmod2,CDmod3,CDmod4)); print(CDresults$table)
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2
#===========================================================
#  test multi-season  false-positive model with expected-value data

# generate expected value data for multi-season model with false positive detections
hstfrq=genpresEV(N=100, K=6, psi=.75, p=.4, fp=0.3, b=0.6, eps=.2, gam=.1, sps=2)  # generate det. hists and freqs
EVpao=createPao(hstfrq$hst,frq=hstfrq$frq,nsurveyseason=rep(2,(ncol(hstfrq$hst)+1)/2))                  #  create pao file

EVmod2<-occMod(model=list(psi~1,p11~1,p10~1,b~1,gamma~1,epsilon~1),data=EVpao,type="do.fp",outfile="savepao.out")  # generate results
#  compare specified input parameters with estimated parameters...
reslt=cbind(c(.75,.4,.3,.6,.2,.1,.1),
  c(EVmod2$real$psi$est[1],EVmod2$real$p11$est[1],EVmod2$real$p10$est[1],EVmod2$real$b$est[1],
    EVmod2$real$epsilon$est[1],EVmod2$real$gamma$est[c(1,EVpao$nunits+1)]))
rownames(reslt)=c('psi','p11','p10','b','epsilon','gamma1','gamma2'); colnames(reslt)=c('input_parameter','estimated_parameter')
cat('testing multi-season model (using expected values)...\n'); print(round(reslt,4))
sumerr=round(sum(reslt[,1]-reslt[,2]),2)
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#   Exanple 1 - Single-season model with Salamander data (pg 135 in new occ. book)
# load a csv file with detection-histories
salmdr.csv<-read.csv(system.file("extdata/Blue_Ridge_pg135.csv",package="RPresence"),header=FALSE)

# Create PRESENCE input file object from csv
salmdr.data<-createPao(salmdr.csv,paoname="salmdr.pao")

## fit some models
salmod1<-occMod(model=list(psi~1,p~SURVEY),data=salmdr.data,type="so",VCoutopt="noreal")
salmod2<-occMod(model=list(psi~1,p~1)     ,data=salmdr.data,type="so")

## create AIC table
models<-list(salmod1,salmod2)
results<-createAicTable(models)
cat('Blue-Ridge example (table 4.2 in book)\n'); print(summary(results))
salmod1_psi.real=predict(salmod1,salmdr.data$survcov,param="psi",conf=0.95);
print(head(salmod1_psi.real))
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2
#===========================================================
#  test false-positive model with expected-value data

hstfrq=genpresEV(N=1000, K=4, psi=.8, p=.3, fp=0.1, b=0.5)  #  generate expected-value data...
EVpao=createPao(hstfrq$hst,frq=hstfrq$frq)                  #  create pao file
EVmod1<-occMod(model=list(psi~1,p11~1,p10~1,b~1),data=EVpao,type="so.fp",outfile='modname')  # generate results
#  compare specified input parameters with estimated parameters...
reslt=cbind(c(.8,.3,.1,.5),c(EVmod1$real$psi$est[1],EVmod1$real$p11$est[1],EVmod1$real$p10$est[1],EVmod1$real$b$est[1]))
rownames(reslt)=c('psi','p11','p10','b'); colnames(reslt)=c('input_parameter','estimated_parameter')
cat('testing false-positive model (using expected values)...\n'); print(round(reslt,4))
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#  test multi-season model with expected-value data

hstfrq=genpresEV(N=1000, K=4, psi=.8, p=.3, eps=.2, gam=.1, sps=2)  #  generate expected-value data...
EVpao=createPao(hstfrq$hst,frq=hstfrq$frq,nsurveyseason = rep(2,2))                  #  create pao file
iv=qlogis(c(.8,.1,.2,.3));
EVmod2<-occMod(model=list(psi~1,p~1,gamma~1,epsilon~1),data=EVpao,type="do.1",initvals=iv,maxfn=1,outfile="modname",threads=2)  # generate results
#  compare specified input parameters with estimated parameters...
reslt=cbind(c(.8,.3,.2,.1),c(EVmod2$real$psi$est[1],EVmod2$real$p$est[1],EVmod2$real$epsilon$est[1],EVmod2$real$gamma$est[1]))
rownames(reslt)=c('psi','p','epsilon','gamma'); colnames(reslt)=c('input_parameter','estimated_parameter')
cat('testing multi-season model (using expected values)...\n'); print(round(reslt,4))
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#  test multi-season model using "Skinks" data

data<-readPao(system.file('extdata/skinks.pao',package="RPresence")) # read pao object (created w/ PRESENCE)
data$unitcov$habitat[data$unitcov$Tussock==1]<-"tussock"             # create "habitat" covariate
data$unitcov$habitat[data$unitcov$Pasture==1]<-"pasture"

m0<-occMod(model=list(psi~1,gamma~1,epsilon~1,p~1),data=data,type="do.1");  #  run model (psi()gam()eps()p())
m1<-occMod(model=list(psi~habitat,gamma~1,epsilon~1,p~1),data=data,type="do.1"); # run model psi(hab)...
m2<-occMod(model=list(psi~habitat,gamma~SEASON,epsilon~1,p~1),data=data,type="do.1",outfile='modname',maxfn="32123",chsq=TRUE); # run model psi(hab)...

cat('skinks AIC:',m0$aic,'\n') #print(summary(m0))
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2


#===========================================================
#  test single-season, two-species model (data unknown)
# load a csv data file
dethist<-read.csv(system.file("extdata/twosp_exmpl.csv",package="RPresence"),as.is=TRUE,header=FALSE)
##
nsites=nrow(dethist); nsrvys=ncol(dethist)         #  set number of sites,surveys from det. history data
dethist=matrix(as.integer(unlist(dethist)),nrow=nsites) # replace missing values (-) with NA
cov1=cov2=NULL

##          create input "pao" object, for use with occMod function
data=createPao(dethist,unitcov=cov1,survcov=cov2,title="twosp example")

## fit some models
mod1<-occMod(model=list(psi~SP,p~SP),data=data,type="so.2sp.1",param="PsiBA",outfile='crap1.out')
fixed=data.frame(param=c('nu',paste0('rho[',1:5,']')),value=rep(1,6))
mod2<-occMod(model=list(psi~SP,p~SP),data=data,type="so.2sp.2",param="nu",fixed=fixed,outfile='crap.out')
#mod2<-occMod(model=list(psi~SP,p~SP),data=data,type="so.2sp.2",param="nu")
cat('two species AIC=',mod1$aic,'\n')
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#  bobcat_example.r

#    read detection history data from csv file,
csv<-read.csv("h:/x/workshops/Smithsonian2016CameraTrap/exercises/Bobcat-3day.csv",as.is=TRUE)

dethist=csv[,-1];  #  get rid of 1st column (site name)
sitenames=csv[,1]  #  sitenames in 1st column
nsites=nrow(dethist); nsrvys=ncol(dethist)  #  set number of sites,surveys from det. history data

#  read covariate table, exclding 1st column (sitenames)
cov1=read.csv("h:/x/workshops/Smithsonian2016CameraTrap/exercises/Bobcat_site_covar_3Day_class.csv",as.is=TRUE,nrows=nsites)[,-1]
cov1=as.data.frame(cov1)
cov2=data.frame(SRVY=as.factor(rep(1:nsrvys,each=nsites)))

#          create input "pao" object, for use with occMod function
data=createPao(dethist,unitcov=cov1,survcov=cov2,title="Bobcat example",unitnames=sitenames)

mods=list(); i=1
mods[[i]]=occMod(model=list(psi~1,p~1),            data=data,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1,p~Trail),data=data,paoname=NULL,type="so");i=i+1
mods[[i]]=occMod(model=list(psi~1,p~Trail+Det_dist),data=data,paoname=NULL,type="so");i=i+1
#     create AIC table of model results and print
results1<-createAicTable(mods); cat('bobcat example\n'); print(results1$table);
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#  exp. value data - multi-state, psi/phi parameterization...
t=Sys.time()
# generate expected value data for multi-season model (phi parameterization)
K=8; nstates=3; p=matrix(0,nstates,nstates); p[2,2]=p[3,2]=p[3,3]=.25; psi1=c(.167,.5,.333)
phi=rbind(c(.5,.333),c(.5,.333),c(.5,.333)); phi=cbind(1-rowSums(phi),phi)

x4=genpresEVms(N=100,K=K,nstates=3,psi1=psi1,phi=phi,p=p,sps=2)  # generate exp. value data using RPresence function

data=createPao(x4$hst,frq=x4$frq,nsurveyseason=rep(2,K/2))
fixed=data.frame(param=paste0('p21(',1:K,')'),value=rep(0,K))
t=Sys.time()
ymods=occMod(model=list(psi~STATE,phi~FROM_TO,p~OTSTATE),data=data,fixed=fixed,type="do.ms.1",outfile="modname",threads=2)
#print(round(ymods$real$psi[seq(1,nrow(ymods$real$psi),data$nunits),],6))
#print(round(ymods$real$phi[seq(1,nrow(ymods$real$phi),data$nunits),],6))
#print(round(ymods$real$p[seq(1,nrow(ymods$real$p),data$nunits),],6))
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
#  Single-season multi-state example (3 occupancy states, psi-R parameteization) (pg 231 in book)
t=Sys.time()
#    read detection history data from csv file...
csv<-read.csv(system.file("extdata/cal_owl_multistate_data.csv",package="RPresence"),as.is=TRUE)
csv[csv=="."]=NA

sitenames=csv[,1]  #  sitenames in 1st column
dethist=csv[,-1];  #  get rid of 1st column (site name)
nsites=nrow(dethist) #  set number of sites,surveys from det. history data
nsrvys=ncol(dethist)

#  create survey covariate to categorize surveys 1-2 as 1 period
#                                    and surveys 3-5 as another period
#    Since it's a survey covariate, it will be a NxT matrix (N=nsites, T=nsurveys)
#    The covariate matrix would be:
#         1 2 3 4 5
#  site1 [1 1 2 2 2]
#  site2 [1 1 2 2 2]
#   :    [: : : : :]
#  siteN [1 1 2 2 2]
#    Filling in by cols means we repeat "1" N times, then "1" N times,
#                  then "2" N times, then "2" N times, then "2" N times.
#    We save it as a data frame (without the matrix dimensions).
cov2=data.frame(PER=as.factor(c(rep(1,2*nsites),rep(2,3*nsites))))

#          create input "pao" object, for use with occMod function
data=createPao(dethist,survcov=cov2,title="Cal Owl example")

xmods=list(); i=1       #  run each model and save in list variable, "xmods"
xmods[[i]]=occMod(model=list(psi~1,r~1,p~1,     delta~1),data=data,type="do.ms.2",outfile='tmp.out');i=i+1
xmods[[i]]=occMod(model=list(psi~1,r~1,p~SURVEY,delta~1),data=data,type="do.ms.2",outfile='tmp.out');i=i+1
xmods[[i]]=occMod(model=list(psi~1,r~1,p~STATE, delta~1),data=data,type="do.ms.2",outfile='tmp.out');i=i+1
xmods[[i]]=occMod(model=list(psi~1,r~1,p~1,     delta~PER),data=data,type="do.ms.2",outfile='tmp.out');i=i+1
xmods[[i]]=occMod(model=list(psi~1,r~1,p~SURVEY,delta~PER),data=data,type="do.ms.2",outfile='tmp.out');i=i+1
xmods[[i]]=occMod(model=list(psi~1,r~1,p~STATE, delta~PER),data=data,type="do.ms.2");i=i+1

#     create AIC table of model results and print
results2<-createAicTable(xmods); cat('Cal Owl example\n'); print(results2$table)

#     print table 5.1 from book...
cat('CA spotted owl reproduction Table 5.1 (pg 234 in book)\n')
estmts=xmods[[6]]$real
estimate_table=rbind(estmts$psi[1,],estmts$r[1,],estmts$p1[1,],estmts$p2[1,],estmts$delta[1,],estmts$delta[109,])
rownames(estimate_table)=c('psi','R','p1','p2','delta1','delta2'); print(estimate_table)
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2

#===========================================================
# single-season 2-species corr. det. model
setwd('h:/a/aatmp')
s="c:/hines/y/presence2/genpres/genpres8 N=100 T=4 psiA=.8, psiBA=.3 psiBa=.7 thA=.4 thA'=.9"
s=paste(s,"thBA=.3 thBA'=.9 thBa=.3 thBa'=.9 pA=.6")
s=paste(s,"pB=.66 rA=.7 rBA=.4 rBa=.5 QUIET P0MISS 2SP nrand=0")
print(s)
i=system(s,intern=TRUE); print(i)
library(RPresence)
pao=readPao('genpres.pao')
amod=occMod(model=list(psi~-1+SP+INT,theta~-1+ALLDIFF,p~SP+INT_o+INT_d+INT_so,th0pi~1),data=pao,type="so.2sp.cd",outfile='tmp.out')
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2
#===========================================================

# multi-season-multi-state example using 1st parameterization (psi,phi,p)
t=Sys.time()
genpres='c:/progra~1/presence/genpres8.exe';        #  use Genpres to generate expected-value data...
if (!file.exists(genpres)) genpres='c:/progra~2/presence/genpres8.exe'
if (!file.exists(genpres)) genpres='h:/y/presence2/genpres/genpres8'
nsrvys=6;
s=paste(c(
  paste0('N=100 T=',nsrvys,' Psi=0.500,0.333'),       # 100 sites, 6 surveys, init occupancy=(.5, .333)
  'phi(0-1)=0,0.5,0,0.5,0,',                          # trans. prob (not occ to occ-state-1) = .5 between seasons
  'phi(0-2)=0,0.3,0,0.3,0,',                          # trans. prob (not occ to occ-state-2) = .3 between seasons
  'phi(1-1)=1,0.6,1,0.6,1,',                          # trans. prob (occ-st-1 to occ-state-1) = .6 between seasons
  'phi(1-2)=0,0.311,0,0.311,0,',                      # trans. prob (occ-st-1 to occ-state-2) = .311 between seasons
  'phi(2-1)=0,0.55,0,0.55,0,',                        # trans. prob (occ-st-2 to occ-state-1) = .55 between seasons
  'phi(2-2)=1,0.377,1,0.377,1,',                      # trans. prob (occ-st-2 to occ-state-2) = .377 between seasons
  paste0('p(1,1)=',paste(rep('0.25',nsrvys),collapse=',')),     #  Pr(detect state 1, given true state=1) = .25
  paste0('p(1,2)=',paste(rep('0.25',nsrvys),collapse=',')),     #  Pr(detect state 1, given true state=2) = .25
  paste0('p(2,1)=',paste(rep('0.00',nsrvys),collapse=',')),     #  Pr(detect state 2, given true state=1) = .00
  paste0('p(2,2)=',paste(rep('0.25',nsrvys),collapse=',')),     #  Pr(detect state 2, given true state=2) = .25
  'QUIET P0MISS  PAOF NSTATES=3 SRVYPERSEASN=2,2,2 nrand=0'),collapse=' ')
cat(genpres,s,'\n'); ii=system(paste(genpres,s),intern=TRUE)

library(RPresence)

data=readPao('genpres.pao')                    #  Create Presence input object from Genpres output file (genpres.pao)
fixed=data.frame(param=paste0('p21(',1:nsrvys,')'),value=rep(0,6))  #  fix Pr(detect state=2, true state=1)=0

#   run a model w/ psi,phi different for each transition, p different for each obs and true state...
zmods=occMod(model=list(psi~-1+STATE,phi~-1+FROM_TO,p~-1+OTSTATE),data=data,type="do.ms.1",fixed=fixed,outfile="modname")
#print(zmods$real)    # print output parameters, psi,phi,p
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2
#===========================================================
#   2 species single season false positive model...
y=read.table('h:/y/presence2/twospfp_detdata.csv',header=FALSE,sep=',')
w=read.table('h:/y/presence2/twospfp_confdata.csv',header=FALSE,sep=',')
pao=createPao(data=data.frame(y),survcov=data.frame(conf=as.numeric(unlist(w))),paoname='twospfp.pao')
m2fp<-occMod(list(psi~SP+INT,p~SP+INT_o+INT_so,omeg~SP,conf~SP),type='so.2sp.fp',data=pao,outfile='modname')
t2=Sys.time(); cat('time:',t2-t,'\n'); steptime=c(steptime,t2-t); t=t2


#========================================================

cat('\nTesting results match expectations...\n=====================================\n')
test_that("model results", {
  expect_equal(round(m1$aic,3),1751.535)
  expect_equal(round(m1$real$psi$est[1],4),0.2211)
  expect_equal(round(m1$real$psi$se[1],4),0.0400)
  expect_equal(round(EVmod1$aic,3),6017.259)
  expect_equal(round(EVmod1$real$psi$est[1],4),0.8)
  expect_equal(round(EVmod1$real$p11$est[1],4),0.3)
  expect_equal(round(EVmod1$real$p10$est[1],4),0.1)
  expect_equal(round(EVmod1$real$b$est[1],4),0.5)
  expect_equal(round(EVmod2$aic,3),4166.165)
  expect_equal(round(EVmod2$real$psi$est[1],4),0.8)
  expect_equal(round(EVmod2$real$p$est[1],4),0.3)
  expect_equal(round(EVmod2$real$eps$est[1],4),0.2)
  expect_equal(round(EVmod2$real$gam$est[1],4),0.1)
  expect_equal(round(EVmod2$derived$psi$est[1],4),0.66)
  expect_equal(round(EVmod2$derived$psi$se[1],4),0.0469)
  expect_equal(round(mod1$aic,3), 1172.482)
  expect_equal(round(mod1$real$psiA$est[1],4), 0.7171)
  expect_equal(round(mod1$real$psiA$se[1],4), 0.0668)
  expect_equal(round(results$table$DAIC[2],2),1.96)
  expect_equal(round(salmod1$aic,4),167.7144)
  expect_equal(round(salmod2$aic,4),165.7586)
  expect_equal(round(salmod1_psi.real$est[1],4),0.5799)
  expect_equal(round(salmod1_psi.real$se[1],4),0.1176)
  expect_equal(round(mods[[1]]$aic,3),2996.497)
  expect_equal(round(results1$table$DAIC[3],4),146.7203)
  expect_equal(round(xmods[[4]]$real$psi$est[1],2),0.98)
  expect_equal(round(xmods[[4]]$real$r$est[1],2),0.45)
  expect_equal(round(xmods[[4]]$real$delta$est[1],2),0.0)
  expect_equal(round(xmods[[4]]$real$delta$est[109],2),0.87)
  expect_equal(round(xmods[[4]]$real$p1$est[1],2),0.74)
  expect_equal(round(CDmod1$aic,3),438.354)
  expect_equal(round(amod$aic,3),706.729)
  expect_equal(round(m2fp$aic,3),1405.587)
  expect_equal(sumerr,0)
})

