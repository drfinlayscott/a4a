###############################################################################
# EJ(20111117)
# A4A assessment model 1: 
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
setGeneric('a4a1', function(pars, stk, idx, ...){standardGeneric('a4a1')})
setGeneric('fobj', function(par, Cayobs, Iayobs, May, nages, nyrs){standardGeneric('fobj')})

#==============================================================================
# fobj
#==============================================================================

setMethod("fobj", c("vector", 'FLQuant', 'FLQuant', "FLQuant", 'vector', 'vector'), 
	function(par, Cayobs, Iayobs, May, nages, nyrs) {
		flq <- Cayobs
		flc <- FLCohort(flq)
		Sa <- par[1:nages]
		Fy <- par[nages+1:nyrs]
		Ry <- exp(par[nages+nyrs+1:nyrs]*10)
		qa <- log(par[nages+2*nyrs+1:(nages-2)])*1e-5 # survey only has 8 years

		flq[] <- Sa
		Sa <- flq
		flq[1,] <- Fy
		Fy <- flq[rep(1,nages)]
		flc[1,1:(nages-1)] <- median(Ry)
		flc[1,nages+1:nyrs-1] <- Ry
		Ry <- flc[rep(1,nages)]
		flq[1:8] <- qa
		qa <- flq[1:8]

		Fay <- Fy*Sa 
		Zay <- Fay+May
		Zay <- FLCohort(Zay)
		Zay[is.na(Zay)] <- 0 # trick to use cumsum which doesn't accept na.rm=T
		cumZ <- apply(Zay, 2, cumsum)	
		muay <- FLCohort(Fay)/Zay*(1-exp(-Zay))
		Cayobs <- FLCohort(Cayobs)				
		Cayhat <- Ry[drop=T]*muay[drop=T]*exp(-cumZ)
		Iayobs <- FLCohort(Iayobs)
		Iayhat <- FLCohort(qa)[drop=T]*(Ry[drop=T]*exp(-cumZ))[1:8]
		Iayhat <- Iayhat[,dimnames(Iayobs)[[2]]]

		sum((log(Cayobs[drop=T])-log(Cayhat))^2, na.rm=T) + sum((log(Iayobs[drop=T])-log(Iayhat))^2, na.rm=T)  	
	}
)

#==============================================================================
# a4a1
#==============================================================================

setMethod("a4a1", c("FLPar", "FLStock", "FLIndex"), 
	function(pars, stk, idx, lower, upper, hessian=FALSE) {
		pars <- c(pars)
		Cayobs <- catch.n(stk)
		sdim <- dim(Cayobs)
		nyrs <- sdim[2]
		nages <- sdim[1]
		May <- m(stk)
		Iayobs <- index(idx)
		ox <- optimx(pars, fobj, Cayobs=Cayobs, Iayobs=Iayobs, May=May, nyrs=nyrs, nages=nages, method="nlminb", lower=lower, upper=upper)
	}
)


