# setwd("/home/darryl/Documents/CUSUM paper/weta")
# library(RPresence)
# 
# 
# # load a PRESENCE data file
# filename<-system.file("extdata/weta.pao",package="RPresence")
# weta.data<-read.pao(filename)
# 
# ## convert indicator variables to categorical covariates
# weta.data$unitcov$Habitat[weta.data$unitcov$Browsed==1] <- "browsed"
# weta.data$unitcov$Habitat[weta.data$unitcov$Unbrowsed==1] <- "unbrowsed"
# weta.data$survcov$Obs<-as.factor(
#   1*weta.data$survcov$Obs1+
#     2*weta.data$survcov$Obs2+
#     3*weta.data$survcov$Obs3)
# 
# ## fit some models
# mod1<-occ.mod(model=list(psi~Habitat,p~SURVEY+Obs),data=weta.data,type="so")
# mod2<-occ.mod(model=list(psi~Habitat,p~Obs),data=weta.data,type="so")
# mod3<-occ.mod(model=list(psi~Habitat,p~SURVEY),data=weta.data,type="so")
# mod4<-occ.mod(model=list(psi~Habitat,p~1),data=weta.data,type="so")
# mod5<-occ.mod(model=list(psi~1,p~SURVEY+Obs),data=weta.data,type="so")
# mod6<-occ.mod(model=list(psi~1,p~Obs),data=weta.data,type="so")
# mod7<-occ.mod(model=list(psi~1,p~SURVEY),data=weta.data,type="so")
# mod8<-occ.mod(model=list(psi~1,p~1),data=weta.data,type="so")
# 
# r.8<-resid.test(mod8)
# boot.8<-bootstrap.residuals(mod8,100,weta.data)
# 
# jpeg("CUSUM_constant_%d.jpg",width=800,height=400,res=144)
# par(mar=c(4,4,0,1)+0.3)
# cusum.plot(weta.data$unitcov$Habitat,r.8$occ,boot.8$occ,ylab="Occupancy Residual",xlab="Habitat",xaxt="n",
#            ylim=c(-10,10))
# axis(1,at=1:2,labels=c("Browsed","Unbrowsed"))
# 
# cusum.plot(weta.data$survcov$Obs,r.8$det,boot.8$det,ylab="Detection Residual",xlab="Observer",xaxt="n",
#            ylim=c(-20,20))
# axis(1,at=1:3,labels=LETTERS[1:3])
# 
# cusum.plot(weta.data$survcov$SURVEY,r.8$det,boot.8$det,ylab="Detection Residual",xlab="Survey",xaxt="n",
#            ylim=c(-20,20))
# axis(1,at=1:5,labels=1:5)
# 
# cusum.plot(rep(weta.data$unitcov$Habitat,5),r.8$det,boot.8$det,ylab="Detection Residual",xlab="Habitat",xaxt="n",
#            ylim=c(-25,25))
# axis(1,at=1:2,labels=c("Browsed","Unbrowsed"))
# dev.off()
# 
# r.1<-resid.test(mod1)
# boot.1<-bootstrap.residuals(mod1,100,weta.data)
# 
# jpeg("CUSUM_full_%d.jpg",width=800,height=400,res=144)
# par(mar=c(4,4,0,1)+0.3)
# 
# cusum.plot(weta.data$unitcov$Habitat,r.1$occ,boot.1$occ,ylab="Occupancy Residual",xlab="Habitat",xaxt="n",
#            ylim=c(-10,10))
# axis(1,at=1:2,labels=c("Browsed","Unbrowsed"))
# 
# cusum.plot(weta.data$survcov$Obs,r.1$det,boot.1$det,ylab="Detection Residual",xlab="Observer",xaxt="n",
#            ylim=c(-20,20))
# axis(1,at=1:3,labels=LETTERS[1:3])
# 
# cusum.plot(weta.data$survcov$SURVEY,r.1$det,boot.1$det,ylab="Detection Residual",xlab="Survey",xaxt="n",
#            ylim=c(-20,20))
# axis(1,at=1:5,labels=1:5)
# 
# cusum.plot(rep(weta.data$unitcov$Habitat,5),r.1$det,boot.1$det,ylab="Detection Residual",xlab="Habitat",xaxt="n",
#            ylim=c(-25,25))
# axis(1,at=1:2,labels=c("Browsed","Unbrowsed"))
# dev.off()
