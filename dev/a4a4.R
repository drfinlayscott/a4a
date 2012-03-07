###############################################################################
# EJ(20111212)
# A4A assessment model 4 (same as 2 but without FLR objects): 
#	Type: Separable SCAA 
#	Data:
#		catch-at-age
#		1 index of abundance, 
#	Error:
#		catch-at-age
#		abundance-at-age  
#	Fit: minimum squares	
#
###############################################################################

#==============================================================================
# Generic
#==============================================================================
setGeneric('a4a4', function(stk, idx, ...){standardGeneric('a4a4')})

#==============================================================================
# fobj
#==============================================================================

#------------------------------------------------------------------------------
# second attempt ...
#------------------------------------------------------------------------------

fobj2 <- function(par, flqini, flcini, Cayobs, Iayobs, May, Cnages, Cnyrs, Cageidx, Cyrsidx, Iageidx, Ichtidx, parscale=rep(1,4)) {
		Sa <- par[Cageidx]*parscale[1] # scaling
		Fy <- par[Cnages+Cyrsidx]*parscale[2] # scaling
		Ry <- exp(par[Cnages+Cnyrs+Cyrsidx])*parscale[3] # scaling
		qa <- log(par[Cnages+2*Cnyrs+Iageidx])*parscale[4] # scaling

		flcini[] <- Sa
		Sa <- flcini

		flqini[1,] <- Fy
		Fy <- flqini[rep(1,Cnages)]
		Fy <- FLCohort(Fy)@.Data[,,1,1,1,1, drop=TRUE] # <<<=== INEFICIENT
		
		flcini[1,1:(Cnages-1)] <- median(Ry)
		flcini[1,Cnages+Cyrsidx-1] <- Ry
		Ry <- flcini[rep(1,Cnages),]

		flcini <- flcini[Iageidx, Ichtidx] # could be a FLCohort
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
# a4a3
#==============================================================================

setMethod("a4a4", c("FLStock", "FLIndex"), 
	function(stk, idx, hessian=FALSE) {
		# data
		cth <- catch.n(stk)
		flqini <- FLQuant(NA, dimnames=dimnames(cth))
		Cnages <- stk@range['max']-stk@range['min']+1
		Cnyrs <- stk@range['maxyear']-stk@range['minyear']+1
		Cageidx <- 1:Cnages
		Cyrsidx <- 1:Cnyrs
		Cayobs <- FLCohort(cth)
		flcini <- FLCohort(NA, dimnames=dimnames(Cayobs))
		May <- FLCohort(m(stk))
		Iayobs <- FLCohort(index(idx))
		Inages <- idx@range['max']-idx@range['min']+1
		# To subset the catch matrix matching the index
		Cdnms <- dimnames(Cayobs)
		Idnms <- dimnames(Iayobs)
		Iageidx <- match(Idnms[[1]], Cdnms[[1]])
		Ichtidx <- match(Idnms[[2]], Cdnms[[2]])
		# starting values
		Saini <- c(seq(0.1,1,len=floor(Cnages/3)), rep(1,Cnages-floor(Cnages/3)))
		Fyini <- rep(0.5, Cnyrs)
		Ryini <- log(cth[1]/(0.1*(1-exp(-1)))) #M0=0.9 F0=0.1
		qaini <- exp(rep(10e-6,Inages))
		par <- c(Saini, Fyini, Ryini, qaini)
		# bounds - WARNING transformation on pars may require adjustement
		low <- c(rep(1e-10, Cnages), rep(1e-10, Cnyrs), log(rep(1e-10, Cnyrs)), exp(rep(1e-10, Inages)))
		upp <- c(rep(1, Cnages), rep(2, Cnyrs), log(rep(Inf, Cnyrs)), exp(rep(1e-3, Inages)))
		#ox <- optimx(par, fobj, flqini=flqini, flcini=flcini, Cayobs=Cayobs, Iayobs=Iayobs, May=May, Cnages=Cnages, Cnyrs=Cnyrs, Cageidx=Cageidx, Cyrsidx=Cyrsidx, Iageidx=Iageidx, Ichtidx=Ichtidx, parscale=rep(1,4), method='nlminb', lower=low, upper=upp, hessian=FALSE, itnmax=1000) 
		flcini=flcini@.Data[,,1,1,1,1, drop=TRUE]
		Cayobs=Cayobs@.Data[,,1,1,1,1, drop=TRUE]
		Iayobs=Iayobs@.Data[,,1,1,1,1, drop=TRUE]
		May=May@.Data[,,1,1,1,1, drop=TRUE]
		
		nlminb(par, fobj2, hessian=FALSE, flqini=flqini, flcini=flcini, Cayobs=Cayobs, Iayobs=Iayobs, May=May, Cnages=Cnages, Cnyrs=Cnyrs, Cageidx=Cageidx, Cyrsidx=Cyrsidx, Iageidx=Iageidx, Ichtidx=Ichtidx, parscale=rep(1,4), lower=low, upper=upp, control=list(iter.max=10000, eval.max=10000)) 
	}
)


