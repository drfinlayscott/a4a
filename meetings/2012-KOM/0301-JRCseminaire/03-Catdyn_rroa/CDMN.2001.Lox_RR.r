DDir <- "C:/ruben/temp/tese2011/data"
WDir <- "C:/ruben/temp/tese2011/modeling"
setwd(DDir)
load("C:\\Ruben\\Temp\\TESE2011\\Data\\Lox.Database.1996.2011.RHRU.AZ.CM.JH.2.RData")
setwd(WDir)
################################################################################
library(CatDyn)
################################################################################
##################                   2001              #########################
################################################################################

SeasonData.2001.tot <- data.tot[data.tot$Year==2001,]
SeasonData.2001.sam <- data.sam[data.sam$Year==2001,]
plot(SeasonData.2001.sam$Eff.h,SeasonData.2001.tot$Eff.h) 

SeasonData.2001.tot <- SeasonData.2001.tot[,c("Week","Catch.m","Eff.h","Eff.h","MBm.g")]
names(SeasonData.2001.tot) <- c("period","obscat","obseff1","obseff2","obsmbm")
class(SeasonData.2001.tot) <- "CatDynData"
plot(x=SeasonData.2001.tot,
     tstep="Week",
     mult="Millions",
     unit1="Hours Diving",
     unit2="Hours Diving",
     bmunit="g",
     span=1,
     top.text="Estimated Effort 2001",
     hem='S')

SeasonData.2001.tot <- data.tot[data.tot$Year==2001 & data.tot$Week >13,]
SeasonData.2001.tot <- SeasonData.2001.tot[,c("Week","Catch.m","Eff.h","Eff.h","MBm.g")]
names(SeasonData.2001.tot) <- c("period","obscat","obseff1","obseff2","obsmbm")
class(SeasonData.2001.tot) <- "CatDynData"
plot(x=SeasonData.2001.tot,
     tstep="Week",
     mult="Millions",
     unit1="Hours Diving",
     unit2="Hours Diving",
     bmunit="g",
     span=1,
     top.text="Observed total 2001",
     hem='S')
###0 Perturbation Model
M             <- 0.00001
N0.ini        <- 455
k.ini         <- 1.05e-5
alpha.ini     <- 0.85
beta.ini      <- 0.95
pars.ini.0P   <- c(log(M),
                   log(N0.ini),
                   log(k.ini),
                   log(alpha.ini),
                   log(beta.ini))
#Dates
dates.2001.tot.0P <- c(head(SeasonData.2001.tot$period,1),
                       tail(SeasonData.2001.tot$period,1))
sealen.2001.tot   <- tail(dates.2001.tot.0P,1)-head(dates.2001.tot.0P,1)+1
##Catch Dynamics Matrix
tot.2001.0P.ini   <- CDMN0P(par=pars.ini.0P,
                            dates=dates.2001.tot.0P,
                            obscat=SeasonData.2001.tot$obscat,
                            obseff=SeasonData.2001.tot$obseff1,
                            obsmbm=SeasonData.2001.tot$obsmbm,
                            M.fixed=FALSE,
                            distr='normal')
#AIC
AIC.2001.tot.0P.ini <- 2*length(pars.ini.0P) -
                       2*(-((sealen.2001.tot-2)/2)*
                       log(sum(tot.2001.0P.ini$resids^2)))
#Plot ini
plot(x=tot.2001.0P.ini,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.0P.ini$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(AIC.2001.tot.0P.ini,1),
     top.text="Observed Catch-Effort 2001 - 0P Model - Normal",
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=0)
#Fit 1 - Normal
t.start             <- Sys.time()
tot.2001.0P.Normal  <- catdyn(p=0,
                              par=pars.ini.0P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.0P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="normal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 25.91164 secs
#Selected method
norm.meths.2001.0P <- 2
#Predict, plot
tot.2001.0P.Normal.fit   <- CDMN0P(par=log(tot.2001.0P.Normal[[norm.meths.2001.0P]]$bt.par),
                                   dates=dates.2001.tot.0P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='normal')

plot(x=tot.2001.0P.Normal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.0P.Normal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.0P.Normal[[norm.meths.2001.0P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.0P.Normal[[norm.meths.2001.0P]]$Model,tot.2001.0P.Normal[[norm.meths.2001.0P]]$Distr,tot.2001.0P.Normal[[norm.meths.2001.0P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=0)
#Fit 2 - LogNormal
t.start             <- Sys.time()
tot.2001.0P.LogNormal  <- catdyn(p=0,
                              par=pars.ini.0P,
                              itnmax=5000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.0P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="lognormal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 49.79529 secs
#Selected method
lnorm.meths.2001.0P <- 2
#Predict, plot
tot.2001.0P.LogNormal.fit   <- CDMN0P(par=log(tot.2001.0P.LogNormal[[lnorm.meths.2001.0P]]$bt.par),
                                   dates=dates.2001.tot.0P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='lognormal')
plot(x=tot.2001.0P.LogNormal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.0P.LogNormal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.0P.LogNormal[[lnorm.meths.2001.0P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.0P.LogNormal[[lnorm.meths.2001.0P]]$Model,tot.2001.0P.LogNormal[[lnorm.meths.2001.0P]]$Distr,tot.2001.0P.LogNormal[[lnorm.meths.2001.0P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=0)
######################
#1 Perturbation Model
M             <- 0.00001
N0.ini        <- 325.5
P1.ini        <- 145.6
k.ini         <- 0.4e-5
alpha.ini     <- 0.95
beta.ini      <- 0.95
pars.ini.1P   <- c(log(M),
                   log(N0.ini),
                   log(P1.ini),
                   log(k.ini),
                   log(alpha.ini),
                   log(beta.ini))
#Dates
P1.1P             <- 19
dates.2001.tot.1P <- c(head(SeasonData.2001.tot$period,1),
                       P1.1P,
                       tail(SeasonData.2001.tot$period,1))
sealen.2001.tot   <- tail(dates.2001.tot.1P,1)-head(dates.2001.tot.1P,1)+1
##Catch Dynamics Matrix
tot.2001.1P.ini   <- CDMN1P(par=pars.ini.1P,
                            dates=dates.2001.tot.1P,
                            obscat=SeasonData.2001.tot$obscat,
                            obseff=SeasonData.2001.tot$obseff1,
                            obsmbm=SeasonData.2001.tot$obsmbm,
                            M.fixed=FALSE,
                            distr='normal')
#AIC
AIC.2001.tot.1P.ini <- 2*length(pars.ini.1P) -
                       2*(-((sealen.2001.tot-2)/2)*
                       log(sum(tot.2001.1P.ini$resids^2)))
#Plot ini
plot(x=tot.2001.1P.ini,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.1P.ini$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(AIC.2001.tot.1P.ini,1),
     top.text="Raised Catch-Effort 2001 - 1P Model - Normal",
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=P1.1P)
#Fit 1 - Normal
t.start             <- Sys.time()
tot.2001.1P.Normal  <- catdyn(p=1,
                              par=pars.ini.1P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.1P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="normal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 31.93326 secs
#Selected Method
norm.meths.2001.1P <- 2
#Predict, plot
tot.2001.1P.Normal.fit   <- CDMN1P(par=log(tot.2001.1P.Normal[[norm.meths.2001.1P]]$bt.par),
                                   dates=dates.2001.tot.1P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='normal')

plot(x=tot.2001.1P.Normal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.1P.Normal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.1P.Normal[[norm.meths.2001.1P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.1P.Normal[[norm.meths.2001.1P]]$Model,tot.2001.1P.Normal[[norm.meths.2001.1P]]$Distr,tot.2001.1P.Normal[[norm.meths.2001.1P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=P1.1P)
#Fit 2 - LogNormal
P1.1P               <- 21
t.start             <- Sys.time()
tot.2001.1P.LogNormal  <- catdyn(p=1,
                              par=pars.ini.1P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.1P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="lognormal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 8.539975 mins
#Selected Method
lnorm.meths.2001.1P <- 1
#Predict, plot
tot.2001.1P.LogNormal.fit   <- CDMN1P(par=log(tot.2001.1P.LogNormal[[lnorm.meths.2001.1P]]$bt.par),
                                   dates=dates.2001.tot.1P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='lognormal')
plot(x=tot.2001.1P.LogNormal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.1P.LogNormal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.1P.LogNormal[[lnorm.meths.2001.1P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.1P.LogNormal[[lnorm.meths.2001.1P]]$Model,tot.2001.1P.LogNormal[[lnorm.meths.2001.1P]]$Distr,tot.2001.1P.LogNormal[[lnorm.meths.2001.1P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=P1.1P)
#####################
#2 Perturbation Model
M             <- 0.00001
N0.ini        <- 355.5
P1.ini        <- 215.6
P2.ini        <- 125.5
k.ini         <- 0.85e-5
alpha.ini     <- 0.95
beta.ini      <- 0.80
pars.ini.2P   <- c(log(M),
                   log(N0.ini),
                   log(P1.ini),
                   log(P2.ini),
                   log(k.ini),
                   log(alpha.ini),
                   log(beta.ini))
#Dates
P1.2P             <- 19
P2.2P             <- 29
dates.2001.tot.2P <- c(head(SeasonData.2001.tot$period,1),
                       P1.2P,
                       P2.2P,
                       tail(SeasonData.2001.tot$period,1))
sealen.2001.tot   <- tail(dates.2001.tot.2P,1)-head(dates.2001.tot.2P,1)+1
##Catch Dynamics Matrix
tot.2001.2P.ini   <- CDMN2P(par=pars.ini.2P,
                            dates=dates.2001.tot.2P,
                            obscat=SeasonData.2001.tot$obscat,
                            obseff=SeasonData.2001.tot$obseff1,
                            obsmbm=SeasonData.2001.tot$obsmbm,
                            M.fixed=FALSE,
                            distr='normal')
#AIC
AIC.2001.tot.2P.ini <- 2*length(pars.ini.2P) -
                       2*(-((sealen.2001.tot-2)/2)*
                       log(sum(tot.2001.2P.ini$resids^2)))
#Plot ini
plot(x=tot.2001.2P.ini,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.2P.ini$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(AIC.2001.tot.2P.ini,1),
     top.text="Raised Catch-Effort 2001 - 1P Model - Normal",
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(P1.2P,P2.2P))
#Fit 1 - Normal
t.start             <- Sys.time()
tot.2001.2P.Normal  <- catdyn(p=2,
                              par=pars.ini.2P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.2P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="normal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 49.17129 secs
#Selected Method
norm.meths.2001.2P <- 2
#Predict, plot
tot.2001.2P.Normal.fit   <- CDMN2P(par=log(tot.2001.2P.Normal[[norm.meths.2001.2P]]$bt.par),
                                   dates=dates.2001.tot.2P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='normal')

plot(x=tot.2001.2P.Normal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.2P.Normal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.2P.Normal[[norm.meths.2001.2P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.2P.Normal[[norm.meths.2001.2P]]$Model,tot.2001.2P.Normal[[norm.meths.2001.2P]]$Distr,tot.2001.2P.Normal[[norm.meths.2001.2P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(P1.2P,P2.2P))
#Fit 2 - LogNormal
P1.2P             <- 21
P2.2P             <- 35
t.start             <- Sys.time()
tot.2001.2P.LogNormal  <- catdyn(p=2,
                              par=pars.ini.2P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS", "Nelder-Mead"),
                              hessian=TRUE,
                              dates=dates.2001.tot.2P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="lognormal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 1.445343 mins
#Selected  method
lnorm.meths.2001.2P <- 1
#Predict, plot
tot.2001.2P.LogNormal.fit   <- CDMN2P(par=log(tot.2001.2P.LogNormal[[lnorm.meths.2001.2P]]$bt.par),
                                   dates=dates.2001.tot.2P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='lognormal')
plot(x=tot.2001.2P.LogNormal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.2P.LogNormal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.2P.LogNormal[[lnorm.meths.2001.2P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.2P.LogNormal[[lnorm.meths.2001.2P]]$Model,tot.2001.2P.LogNormal[[lnorm.meths.2001.2P]]$Distr,tot.2001.2P.LogNormal[[lnorm.meths.2001.2P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(P1.2P,P2.2P))

#####################
#3 Perturbation Model
M             <- 0.00001
N0.ini        <- 355.5
P1.ini        <- 215.6
P2.ini        <- 125.5
P3.ini        <- 95.5
k.ini         <- 0.55e-5
alpha.ini     <- 0.95
beta.ini      <- 0.80
pars.ini.3P   <- c(log(M),
                   log(N0.ini),
                   log(P1.ini),
                   log(P2.ini),
                   log(P3.ini),
                   log(k.ini),
                   log(alpha.ini),
                   log(beta.ini))
#Dates
P1.3P             <- 20
P2.3P             <- 31
P3.3P             <- 41
dates.2001.tot.3P <- c(head(SeasonData.2001.tot$period,1),
                       P1.3P,
                       P2.3P,
                       P2.3P,
                       tail(SeasonData.2001.tot$period,1))
sealen.2001.tot   <- tail(dates.2001.tot.3P,1)-head(dates.2001.tot.3P,1)+1
##Catch Dynamics Matrix
tot.2001.3P.ini   <- CDMN3P(par=pars.ini.3P,
                            dates=dates.2001.tot.3P,
                            obscat=SeasonData.2001.tot$obscat,
                            obseff=SeasonData.2001.tot$obseff1,
                            obsmbm=SeasonData.2001.tot$obsmbm,
                            M.fixed=FALSE,
                            distr='normal')
#AIC
AIC.2001.tot.3P.ini <- 2*length(pars.ini.3P) -
                       2*(-((sealen.2001.tot-2)/2)*
                       log(sum(tot.2001.3P.ini$resids^2)))
#Plot ini
plot(x=tot.2001.3P.ini,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.3P.ini$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(AIC.2001.tot.3P.ini,1),
     top.text="Raised Catch-Effort 2001 - 1P Model - Normal",
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(P1.3P,P2.3P))
#Fit 1 - Normal
t.start             <- Sys.time()
tot.2001.3P.Normal  <- catdyn(p=3,
                              par=pars.ini.3P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.3P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="normal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 
#Selected method
norm.meths.2001.3P <- 
#Predict, plot
tot.2001.3P.Normal.fit   <- CDMN3P(par=log(tot.2001.3P.Normal[[norm.meths.2001.3P]]$bt.par),
                                   dates=dates.2001.tot.3P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='normal')

plot(x=tot.2001.3P.Normal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.3P.Normal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.3P.Normal[[norm.meths.2001.3P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.3P.Normal[[norm.meths.2001.3P]]$Model,tot.2001.3P.Normal[[norm.meths.2001.3P]]$Distr,tot.2001.3P.Normal[[norm.meths.2001.3P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(P1.3P,P2.3P,P3.3P))
#Fit 2 - LogNormal
t.start             <- Sys.time()
tot.2001.3P.LogNormal  <- catdyn(p=3,
                              par=pars.ini.3P,
                              itnmax=50000,
                              method=c("spg", "CG", "BFGS"),
                              hessian=TRUE,
                              dates=dates.2001.tot.3P,
                              obscat=SeasonData.2001.tot$obscat,
                              obseff=SeasonData.2001.tot$obseff1,
                              M.fixed=FALSE,
                              M=NULL,
                              distr="lognormal")
t.end               <- Sys.time()
(t.process <- t.end-t.start)
#Time difference of 
#Selected method
lnorm.meths.2001.3P <- 
#Predict, plot
tot.2001.3P.LogNormal.fit   <- CDMN3P(par=log(tot.2001.3P.LogNormal[[lnorm.meths.2001.3P]]$bt.par),
                                   dates=dates.2001.tot.3P,
                                   obscat=SeasonData.2001.tot$obscat,
                                   obseff=SeasonData.2001.tot$obseff1,
                                   obsmbm=SeasonData.2001.tot$obsmbm,
                                   M.fixed=FALSE,
                                   distr='lognormal')
plot(x=tot.2001.3P.LogNormal.fit,
     tstep='Week',
     mult='Millions',
     Biom=round(tail(tot.2001.3P.LogNormal.fit$npred,1)*1e6
                *mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6),
     AIC=round(tot.2001.3P.LogNormal[[lnorm.meths.2001.3P]]$AIC,1),
     top.text=paste("Raised Catch-Effort 2001",tot.2001.3P.LogNormal[[lnorm.meths.2001.3P]]$Model,tot.2001.3P.LogNormal[[lnorm.meths.2001.3P]]$Distr,tot.2001.3P.LogNormal[[lnorm.meths.2001.3P]]$method,sep="-"),
     leg.pos='topright',
     AIC.xpos=0.28,
     AIC.ypos=0.1,
     Biom.tstep="fin",
     Biom.xpos=0.28,
     Biom.ypos=0,
     p.dates=c(P1.3P,P2.3P,P3.3P))

################################################################################

norm.meths.2001 <- c(norm.meths.2001.0P,norm.meths.2001.1P,norm.meths.2001.2P,-1,-1)
lnorm.meths.2001 <- c(lnorm.meths.2001.0P,lnorm.meths.2001.1P,lnorm.meths.2001.2P,-1,-1)

Lox.SSA.MLE.2001 <- data.frame(
    Year=    rep(2001,10),
    Model=   rep(c("0P","1P","2P","3P","4P"),2),
    Distr=   c(rep("Normal",5),rep("LogNormal",5)),
    Method=  c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,tot.2001.0P.Normal[[norm.meths.2001[1]]]$method)),                  
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,tot.2001.1P.Normal[[norm.meths.2001[2]]]$method)),                  
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,tot.2001.2P.Normal[[norm.meths.2001[3]]]$method)),                  
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,tot.2001.3P.Normal[[norm.meths.2001[4]]]$method)),                  
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,tot.2001.3P.Normal[[norm.meths.2001[5]]]$method)),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$method)),              
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$method)),              
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$method)),              
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$method)),              
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$method))),
    AIC=     c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$AIC,1))),            
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$AIC,1))),            
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$AIC,1))),            
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$AIC,1))),            
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$AIC,1))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$AIC,1))),        
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$AIC,1))),        
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$AIC,1))),        
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$AIC,1))),        
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$AIC,1)))),
    M=       c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.par[1],12))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.par[1],12))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[1],12))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[1],12))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.par[1],12))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.par[1],12))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.par[1],12))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[1],12))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[1],12))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[1],12)))),
    SE.M=    c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.stdev[1],12))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.stdev[1],12))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[1],12))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[1],12))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.stdev[1],12))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.stdev[1],12))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.stdev[1],12))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[1],12))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[1],12))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[1],12)))),
    N0=      c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.par[2],1))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.par[2],1))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[2],1))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[2],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.par[2],1))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.par[2],1))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.par[2],1))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[2],1))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[2],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[2],1)))),
    SE.N0=   c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.stdev[2],1))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.stdev[2],1))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[2],1))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[2],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.stdev[2],1))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.stdev[2],1))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.stdev[2],1))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[2],1))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[2],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[2],1)))),
    P1=      c(0,                                                                
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.par[3],1))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[3],1))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[3],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.par[3],1))),
               0,                                                                
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.par[3],1))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[3],1))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[3],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[3],1)))),
    SE.P1=   c(0,                                                                
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.stdev[3],1))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[3],1))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[3],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.stdev[3],1))),
               0,                                                                
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.stdev[3],1))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[3],1))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[3],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[3],1)))),
    TS.P1=   c(0,                                                                
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,dates.2001.tot.1P[2])),                                             
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,dates.2001.tot.2P[2])),                                             
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,dates.2001.tot.3P[2])),                                             
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[2])),
               0,                                                                
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,dates.2001.tot.1P[2])),                                             
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,dates.2001.tot.2P[2])),                                             
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,dates.2001.tot.3P[2])),                                             
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[2]))),
    P2=      c(0,                                                                
               0,                                                                
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[4],1))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[4],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.par[4],1))),
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[4],1))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[4],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[4],1)))),
    SE.P2=   c(0,                                                                
               0,                                                                
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[4],1))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[4],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.stdev[4],1))),
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[4],1))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[4],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[4],1)))),
    TS.P2=   c(0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,dates.2001.tot.2P[3])),                                             
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,dates.2001.tot.3P[3])),                                             
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[3])),
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,dates.2001.tot.2P[3])),                                             
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,dates.2001.tot.3P[3])),                                             
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[3]))),
    P3=      c(0,                                                                
               0,                                                                
               0,                                                                
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[5],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.par[5],1))),
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[5],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[5],1)))),
    SE.P3=   c(0,                                                                
               0,                                                                
               0,                                                                
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[5],1))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.stdev[5],1))),
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[5],1))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[5],1)))),
    TS.P3=   c(0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,dates.2001.tot.3P[4])),                                             
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[4])), 
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,dates.2001.tot.3P[4])),                                             
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[4]))),
    P4=      c(0,                                                                
               0,                                                                
               0,                                                                
               0,                                                                
              ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.par[6],1))),
               0,                                                                
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[6],1)))),
    SE.P4=   c(0,                                                                
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[5]]]$bt.stdev[6],1))),
               0,                                                                
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[6],1)))),
    TS.P4=   c(0,                                                                
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[5])), 
               0,                                                                
               0,                                                                
               0,                                                                
               0,                                                                
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,dates.2001.tot.4P[5]))),
    k=       c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.par[3],7))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.par[4],7))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[5],7))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[6],7))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[7],7))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.par[3],7))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.par[4],7))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[5],7))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[6],7))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[7],7)))),
    SE.k=    c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.stdev[3],7))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.stdev[4],7))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[5],7))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[6],7))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[7],7))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.stdev[3],7))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.stdev[4],7))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[5],7))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[6],7))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[7],7)))),
    alpha=   c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.par[4],3))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.par[5],3))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[6],3))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[7],3))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[8],3))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.par[4],3))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.par[5],3))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[6],3))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[7],3))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[8],3)))),
    SE.alpha=c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.stdev[4],3))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.stdev[5],3))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[6],3))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[7],3))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[8],3))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.stdev[4],3))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.stdev[5],3))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[6],3))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[7],3))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[8],3)))),
    beta=    c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.par[5],3))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.par[6],3))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par[7],3))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[8],3))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.par[9],3))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.par[5],3))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.par[6],3))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.par[7],3))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.par[8],3))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.par[9],3)))),
    SE.beta= c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.Normal[[norm.meths.2001[1]]]$bt.stdev[5],3))),      
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.Normal[[norm.meths.2001[2]]]$bt.stdev[6],3))),      
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev[7],3))),      
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[8],3))),      
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.Normal[[norm.meths.2001[4]]]$bt.stdev[9],3))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tot.2001.0P.LogNormal[[lnorm.meths.2001[1]]]$bt.stdev[5],3))),  
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tot.2001.1P.LogNormal[[lnorm.meths.2001[2]]]$bt.stdev[6],3))),  
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tot.2001.2P.LogNormal[[lnorm.meths.2001[3]]]$bt.stdev[7],3))),  
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[4]]]$bt.stdev[8],3))),  
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tot.2001.3P.LogNormal[[lnorm.meths.2001[5]]]$bt.stdev[9],3)))),
    B.esc=   c(ifelse(norm.meths.2001[[1]]==0,0,ifelse(norm.meths.2001[[1]]==-1,-1,round(tail(tot.2001.0P.Normal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(norm.meths.2001[[2]]==0,0,ifelse(norm.meths.2001[[2]]==-1,-1,round(tail(tot.2001.1P.Normal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(norm.meths.2001[[3]]==0,0,ifelse(norm.meths.2001[[3]]==-1,-1,round(tail(tot.2001.2P.Normal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(norm.meths.2001[[4]]==0,0,ifelse(norm.meths.2001[[4]]==-1,-1,round(tail(tot.2001.3P.Normal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(norm.meths.2001[[5]]==0,0,ifelse(norm.meths.2001[[5]]==-1,-1,round(tail(tot.2001.4P.Normal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(lnorm.meths.2001[[1]]==0,0,ifelse(lnorm.meths.2001[[1]]==-1,-1,round(tail(tot.2001.0P.LogNormal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(lnorm.meths.2001[[2]]==0,0,ifelse(lnorm.meths.2001[[2]]==-1,-1,round(tail(tot.2001.1P.LogNormal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(lnorm.meths.2001[[3]]==0,0,ifelse(lnorm.meths.2001[[3]]==-1,-1,round(tail(tot.2001.2P.LogNormal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(lnorm.meths.2001[[4]]==0,0,ifelse(lnorm.meths.2001[[4]]==-1,-1,round(tail(tot.2001.3P.LogNormal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6))),
               ifelse(lnorm.meths.2001[[5]]==0,0,ifelse(lnorm.meths.2001[[5]]==-1,-1,round(tail(tot.2001.4P.LogNormal.fit$npred,1)*1e6*mean(tail(SeasonData.2001.tot$obsmbm,7))*1e-6)))),
    Sel.Mod= c(0,0,1,0,0,0,0,0,0,0))
                                                                                                                                                                                                                             
Lox.SSA.Cor.2001 <- tot.2001.2P.Normal[[norm.meths.2001[3]]]$Cor                                                                                                                                                                  
Lox.SSA.SD.2001  <- tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.stdev
Lox.SSA.Par.2001 <- tot.2001.2P.Normal[[norm.meths.2001[3]]]$bt.par
                                                                                                                                                                                                                           
save.image("CDMN.2001.Lox_RR.RData")

