###############################################################################
# EJ(20111212)
# A4A assessment model 7 (same as 6 but fixed selectivities): 
#	Type: Separable SCAA 
#	Data:
#		catch-at-age
#		1 index of abundance, 
#	Error:
#		catch-at-age
#		abundance-at-age  
#	Fit: minimum squares	
#
#	Sa = linear up to full exploited age, 1 afterwards
#	Isel = q*sel
#	selq = sel0*exp(-beta*t) - exponential decay
###############################################################################

#==============================================================================
# Generic
#==============================================================================
setGeneric('a4a8', function(stk, idx, ...){standardGeneric('a4a8')})

#==============================================================================
# fobj
#==============================================================================

#------------------------------------------------------------------------------
# attempt ...
#------------------------------------------------------------------------------

fobj8 <- function(par, flqini, flcini, Cayobs, Iayobs, May, Cminage, Cmaxage, Cnages, Cnyrs, Cageidx, Cyrsidx, Iminage, Imaxage, Iageidx, Ichtidx) {
		aE <- par[1] # fully exploited age
		Fy <- exp(par[1+Cyrsidx])
		Ry <- exp(par[1+Cnyrs+Cyrsidx])
		lambda <- log(par[1+2*Cnyrs+1])
		q <- exp(-10*par[1+2*Cnyrs+2])

		Sa <- c(seq(0.05,1,length=aE-Cminage), rep(1,Cmaxage-aE+1))
		flcini[] <- Sa
		Sa <- flcini

		flqini[1,] <- Fy
		Fy <- flqini[rep(1,Cnages)]
		Fy <- FLCohort(Fy)@.Data[,,1,1,1,1, drop=TRUE] # <<<=== INEFICIENT
		
		flcini[1,1:(Cnages-1)] <- median(Ry)
		flcini[1,Cnages+Cyrsidx-1] <- Ry
		Ry <- flcini[rep(1,Cnages),]

		qa <- Iminage:Imaxage
		qa <- q*exp(-lambda*qa)
		flcini <- flcini[Iageidx, Ichtidx]
		flcini[] <- qa
		qa <- flcini

		Fay <- Fy*Sa 
		Zay <- Fay+May
		Zay[is.na(Zay)] <- 0 # trick to use cumsum which doesn't accept na.rm=T
		cumZ <- apply(Zay, 2, cumsum) # <<<=== INEFICIENT 	
		muay <- Fay/Zay*(1-exp(-Zay))
		Cayhat <- Ry[drop=T]*muay[drop=T]*exp(-cumZ)
		Iayhat <- qa[drop=T]*((Ry[drop=T]*exp(-cumZ))[Iageidx, Ichtidx])
		e1 <- (log(Cayobs[drop=T])-log(Cayhat))^2
		e1[is.infinite(e1)] <- NA
		e2 <- (log(Iayobs[drop=T])-log(Iayhat))^2
		e2[is.infinite(e2)] <- NA
		sum(e1, na.rm=T) + sum(e2, na.rm=T)  	
	}
#)

#==============================================================================
# a4a
#==============================================================================

setMethod("a4a8", c("FLStock", "FLIndex"), 
	function(stk, idx, hessian=FALSE) {
		# data
		cth <- catch.n(stk)
		flqini <- FLQuant(NA, dimnames=dimnames(cth))
		Cmaxage <- stk@range['max']
		Cminage <- stk@range['min'] 
		Cnages <- Cmaxage-Cminage+1
		Cnyrs <- stk@range['maxyear']-stk@range['minyear']+1
		Cageidx <- 1:Cnages
		Cyrsidx <- 1:Cnyrs
		Cayobs <- FLCohort(cth)
		flcini <- FLCohort(NA, dimnames=dimnames(Cayobs))
		May <- FLCohort(m(stk))
		Iayobs <- FLCohort(index(idx))
		Imaxage <- idx@range['max'] 
		Iminage <- idx@range['min']
		# recruitment indices make a mess with nbeta 0/0 so
		if(Imaxage==0 & Iminage==0) Imaxage <- Iminage <- 0.1 
		Inages <- idx@range['max']-idx@range['min']+1
		# To subset the catch matrix matching the index
		Cdnms <- dimnames(Cayobs)
		Idnms <- dimnames(Iayobs)
		Iageidx <- match(Idnms[[1]], Cdnms[[1]])
		Ichtidx <- match(Idnms[[2]], Cdnms[[2]])
		# starting values
		aE <- floor(Cnages/3)
		Fyini <- rep(0.5, Cnyrs)
		Ryini <- log(cth[1]/(0.1*(1-exp(-1)))) #M0=0.9 F0=0.1
		lambda <- 0.5
		# transformation for q y=-0.1*log(q); q=exp(-10*y)
		q <- exp(10e-10)
		par <- c(aE, Fyini, Ryini, lambda, q)
		# bounds - WARNING transformation on pars may require adjustement
		low <- c(floor(Cnages/3-1), sqrt(rep(1e-2, Cnyrs)), log(rep(min(Ryini)*0.1, Cnyrs)), exp(0.4), -0.1*log(1e-5))
		upp <- c(floor(Cnages/3+1), sqrt(rep(2, Cnyrs)), log(rep(max(Ryini)*10, Cnyrs)), exp(0.6), -0.1*log(1e-15))
		#ox <- optimx(par, fobj, flqini=flqini, flcini=flcini, Cayobs=Cayobs, Iayobs=Iayobs, May=May, Cnages=Cnages, Cnyrs=Cnyrs, Cageidx=Cageidx, Cyrsidx=Cyrsidx, Iageidx=Iageidx, Ichtidx=Ichtidx, parscale=rep(1,4), method='nlminb', lower=low, upper=upp, hessian=FALSE, itnmax=1000) 
		flcini=flcini@.Data[,,1,1,1,1, drop=TRUE]
		Cayobs=Cayobs@.Data[,,1,1,1,1, drop=TRUE]
		Iayobs=Iayobs@.Data[,,1,1,1,1, drop=TRUE]
		May=May@.Data[,,1,1,1,1, drop=TRUE]
		
		#nlminb(par, fobj7, hessian=FALSE, flqini=flqini, flcini=flcini, Cayobs=Cayobs, Iayobs=Iayobs, May=May, Cminage=Cminage, Cmaxage=Cmaxage, Cnages=Cnages, Cnyrs=Cnyrs, Cageidx=Cageidx, Cyrsidx=Cyrsidx, Imaxage=Imaxage, Iminage=Iminage, Iageidx=Iageidx, Ichtidx=Ichtidx, lower=low, upper=upp, control=list(iter.max=10000, eval.max=10000))
		DEoptim(fobj7, flqini=flqini, flcini=flcini, Cayobs=Cayobs, Iayobs=Iayobs, May=May, Cminage=Cminage, Cmaxage=Cmaxage, Cnages=Cnages, Cnyrs=Cnyrs, Cageidx=Cageidx, Cyrsidx=Cyrsidx, Imaxage=Imaxage, Iminage=Iminage, Iageidx=Iageidx, Ichtidx=Ichtidx, lower=low, upper=upp, control=DEoptim.control(strategy=6, c=0.4, trace=FALSE)) 
	}
)


