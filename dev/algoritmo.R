###############################################################################
# EJ(20111003)
# Statistical Catch at Age with FLR
# a=ages
# y=years
# N=population
# R=recruitment
# F=fishing mortality
# M=natural mortality
# C=catch
# S=gear selectivity
# mu=exploitation fraction
# f=net fecundity
# I=abundance index
###############################################################################

library(FLCore)
library(optimx)
library(lattice)
data(ple4)
data(ple4.index)
source("a4a3.R")

#==============================================================================
# Do
#==============================================================================

dnms <- dimnames(catch.n(NSH))
pnms <- c(paste("S", dnms[[1]], sep=""), paste("F", dnms[[2]], sep=""), paste("R", dnms[[2]], sep=""), paste("q", 1:5, sep=""))

# fit 
system.time(
#Rprof()
ple4.ox2 <- a4a4(ple4, ple4.index)
#Rprof(NULL)
)

res <- data.frame(stat=substr(pnms,1,1), val=unlist(nsh3$par), x=as.numeric(substr(pnms,2,5)))
res[res$stat=="R","val"] <- exp(res[res$stat=="R","val"])
res[res$stat=="q","val"] <- log(res[res$stat=="q","val"])

xyplot(val~x|stat, data=res, scales=list(relation="free"), type="l", auto.key=list(points=FALSE, lines=TRUE, space="right"), layout=c(1,4))


out <- benchmark(a4a3(ple4, ple4.index),a4a4(ple4, ple4.index))


#==============================================================================
# apply to herring data and compare to SAM
#==============================================================================
library(FLSAM)
data(NSH)
stck <- NSH
tun  <- NSH.tun
data(NSH.sam)

#-------------------------------------------------------------------------------
# SAM
#-------------------------------------------------------------------------------
ctrl <- NSH.ctrl
sam <- FLSAM(stck,tun,ctrl)
stck.sam <- stck
harvest(stck.sam) <- sam@harvest[,as.character(1960:2010)]
stock.n(stck.sam) <- sam@stock.n[,as.character(1960:2010)]

#-------------------------------------------------------------------------------
# YAS (Yet Another SCAM)
#-------------------------------------------------------------------------------

nsh2 <- a4a4(stck,tun[[2]])

dnms <- dimnames(catch.n(NSH))
pnms <- c(paste("S", dnms[[1]], sep=""), paste("F", dnms[[2]], sep=""), paste("R", dnms[[2]], sep=""), paste("q", dimnames(index(tun[[2]]))[[1]], sep=""))
res <- data.frame(stat=substr(pnms,1,1), val=unlist(nsh2$par), x=as.numeric(substr(pnms,2,5)))
stck2.sca <- stck
harvest(stck2.sca)[] <- t(res[res$stat=="F","val"]%*%t(res[res$stat=="S","val"]))
units(harvest(stck2.sca)) <- 'f'
stock.n(stck2.sca)[1] <- exp(res[res$stat=="R","val"])

nsh3 <- a4a4(stck,tun[[3]])
nsh3.mod5 <- a4a5(stck,tun[[3]])
nsh3.mod6 <- a4a6(stck,tun[[3]])
nsh3.mod7 <- a4a7(stck,tun[[3]])
nsh3.mod8 <- a4a8(stck,tun[[3]])


dnms <- dimnames(catch.n(NSH))
pnms <- c(paste("S", dnms[[1]], sep=""), paste("F", dnms[[2]], sep=""), paste("R", dnms[[2]], sep=""), paste("q", dimnames(index(tun[[3]]))[[1]], sep=""))
res <- data.frame(stat=substr(pnms,1,1), val=unlist(nsh3$par), x=as.numeric(substr(pnms,2,5)))
stck3.sca <- stck
harvest(stck3.sca)[] <- t(res[res$stat=="F","val"]%*%t(res[res$stat=="S","val"]))
units(harvest(stck3.sca)) <- 'f'
stock.n(stck3.sca)[1] <- exp(res[res$stat=="R","val"])

nsh4 <- a4a4(stck,tun[[4]])

dnms <- dimnames(catch.n(NSH))
pnms <- c(paste("S", dnms[[1]], sep=""), paste("F", dnms[[2]], sep=""), paste("R", dnms[[2]], sep=""), paste("q", dimnames(index(tun[[4]]))[[1]], sep=""))
res <- data.frame(stat=substr(pnms,1,1), val=unlist(nsh4$par), x=as.numeric(substr(pnms,2,5)))
stck4.sca <- stck
harvest(stck4.sca)[] <- t(res[res$stat=="F","val"]%*%t(res[res$stat=="S","val"]))
units(harvest(stck4.sca)) <- 'f'
stock.n(stck4.sca)[1] <- exp(res[res$stat=="R","val"])

stks <- FLStocks(sam=stck.sam, sca2=stck2.sca, sca3=stck3.sca, sca4=stck4.sca)


#-------------------------------------------------------------------------------
# a4a9
#-------------------------------------------------------------------------------

obs <- list(catch = Data( stck @ catch.n ),
            index = Data(tun[[3]] @ index)) # 4 has more ages .... but missing age 1s at start of series
M <- Data( stck @ m )

get_stk_funs <-
function(pars, data_summaries)
{
  with(data_summaries,
  list(
    #S = function(a) ifelse(a < pars $ aE, 0.05 + 0.95 * (a - Cminage)/(pars $ aE - Cminage), 1),
    #S = function(a) iglogit( 1 * glogit(a, Cminage-1e-9, Cmaxage+1e-9) + log(99) - glogit(pars $ aE, Cminage-1e-9, Cmaxage+1e-9) * 1, 0.05, 1) ,
    S = function(a) iglogit( (exp(glogit(pars $ aE, Cminage, Cmaxage)) + .5) * 
                             glogit(a, Cminage-1e-9, Cmaxage+1e-9) + 
                             log(99) - glogit(Cmaxage - .5, Cminage, Cmaxage) * (exp(glogit(pars $ aE, Cminage, Cmaxage)) + .5),
                              0.05, 1) ,
    q = function(a) exp(pars $ logq) * exp(- pars $ lambda * a)
  ))
}

fit3 <- SCA.fit(obs, M, trace = 50)






