###############################################################################
# EJ(20111117)
# A4A assessment model 2 (same as 1 but with improved code): 
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
setGeneric('a4a2', function(pars, stk, idx, ...){standardGeneric('a4a2')})
#setGeneric('fobj', function(par, flqini, flcini, Cayobs, Iayobs, May, Cnages, Cnyrs, Cageidx, Cyrsidx, Inages, Inyrs, Iageidx, Iyrsidx, parscale){standardGeneric('fobj')})

#==============================================================================
# fobj
#==============================================================================

#------------------------------------------------------------------------------
# second attempt ...
#------------------------------------------------------------------------------

fobj <- function(par, flqini, flcini, Cayobs, Iayobs, May, Cnages, Cnyrs, Cageidx, Cyrsidx, Iageidx, Ichtidx, parscale=rep(1,4)) {
		Sa <- par[Cageidx]*parscale[1] # scaling
		Fy <- par[Cnages+Cyrsidx]*parscale[2] # scaling
		Ry <- exp(par[Cnages+Cnyrs+Cyrsidx]*parscale[3]) # scaling
		qa <- log(par[Cnages+2*Cnyrs+Iageidx])*parscale[4] # scaling

		#flqini[] <- Sa # could be a FLCohort
		#Sa <- flqini
		flcini[] <- Sa
		Sa <- flcini

		flqini[1,] <- Fy
		Fy <- flqini[rep(1,Cnages)]
		Fy <- FLCohort(Fy)
		
		flcini[1,1:(Cnages-1)] <- median(Ry)
		flcini[1,Cnages+Cyrsidx-1] <- Ry
		Ry <- flcini[rep(1,Cnages)]

		#flqini <- flqini[Iageidx, Iyrsidx] # could be a FLCohort
		#flqini[] <- qa
		#qa <- FLCohort(flqini)
		flcini <- flcini[Iageidx, Ichtidx] # could be a FLCohort
		flcini[] <- qa
		qa <- flcini

		#Fay <- FLCohort(Fy*Sa) 
		Fay <- Fy*Sa 
		Zay <- Fay+May
		Zay[is.na(Zay)] <- 0 # trick to use cumsum which doesn't accept na.rm=T
		cumZ <- apply(Zay, 2, cumsum)	
		muay <- Fay/Zay*(1-exp(-Zay))
		Cayhat <- Ry[drop=T]*muay[drop=T]*exp(-cumZ)
		Iayhat <- qa[drop=T]*((Ry[drop=T]*exp(-cumZ))[Iageidx, Ichtidx])
		sum((log(Cayobs[drop=T])-log(Cayhat))^2, na.rm=T) + sum((log(Iayobs[drop=T])-log(Iayhat))^2, na.rm=T)  	
	}
#)

#==============================================================================
# a4a2
#==============================================================================

setMethod("a4a2", c("FLPar", "FLStock", "FLIndex"), 
	function(pars, stk, idx, lower, upper, hessian=FALSE) {
		par <- c(pars)
		Cayobs <- catch.n(stk)
		flqini <- FLQuant(NA, dimnames=dimnames(Cayobs))
		Cnages <- stk@range['max']-stk@range['min']+1
		Cnyrs <- stk@range['maxyear']-stk@range['minyear']+1
		Cageidx <- 1:Cnages
		Cyrsidx <- 1:Cnyrs
		Cayobs <- FLCohort(Cayobs)
		flcini <- FLCohort(NA, dimnames=dimnames(Cayobs))
		May <- FLCohort(m(stk))
		Iayobs <- FLCohort(index(idx))
		sdim <- dim(Iayobs)
		# To subset the catch matrix matching the index
		Iageidx <- Cageidx[dimnames(Cayobs)[[1]] %in% dimnames(Iayobs)[[1]]]
		Ichtidx <- Cyrsidx[dimnames(Cayobs)[[2]] %in% dimnames(Iayobs)[[2]]]
		ox <- optimx(par, fobj, flqini=flqini, flcini=flcini, Cayobs=Cayobs, Iayobs=Iayobs, May=May, Cnages=Cnages, Cnyrs=Cnyrs, Cageidx=Cageidx, Cyrsidx=Cyrsidx, Iageidx=Iageidx, Ichtidx=Ichtidx, parscale=rep(1,4)) 
	}
)


