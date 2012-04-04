###############################################################################
# EJ CM (20120403)
# A4A assessment model 9 (same as 6 but fixed selectivities): 
#	Type: Separable SCAA 
#	Data:
#		catch-at-age
#		1 index of abundance, 
#	Error:
#		catch-at-age
#		abundance-at-age  
#	Fit: minimum squares with difference observation variance for catch and survey	
#
#	Sa = linear up to full exploited (continuous) age, 1 afterwards
#	Isel = q*sel
#	selq = sel0*exp(-beta*t) - exponential decay
###############################################################################

#==============================================================================
# Generic
#==============================================================================
setGeneric('a4a9', function(stock, index, ...) standardGeneric('a4a9') )

#==============================================================================
# Functions
#==============================================================================

# some utility functions
iglogit <- function(x, min = 0, max = 1) 1/(1 + exp(-x)) * (max - min) + min
glogit <- function(p, min = 0, max = 1) {p <- (p - min) / (max - min); log(p / (1-p))}


# the main fitting function/algorithm
setMethod("a4a9", c("FLStock", "FLIndex"), 
function(stock, index, control.model, trace = 10)
{
  time0 <- proc.time()

  # extract data from stock and index FL objects
  obs <- list(catch = drop( catch.n(stock) @ .Data ),
              index = drop( index(index) @ .Data ))

  M <- drop( m(stock) @ .Data )

  # collect useful data summaries
  data_summaries <- getDataSummaries(obs)
  
  # CHECK CONTROL.MODEL
  
  # get starting values (improve this along the lines of EJ code)
  par.init <- getIntialParams(obs, M, control.model, data_summaries)

  opt <- nlminb(par.init, neg.llik, 
                control.model = control.model, 
                data_summaries = data_summaries, M = M, obs = obs, 
                control = list(trace = trace, iter.max = 1e5, eval.max = 1e5))

  out <- predict_pop(opt $ par, control.model, data_summaries, M, full = TRUE)
  out $ par $ sigma <- exp(opt $ par[1:2])
  out $ opt <- opt
  out $ llik <- -1 * opt $ objective
  out $ aic <- -2 * out $ llik + 2 * length(par.init)
  out $ ellapsed <- proc.time() - time0
  
  out
})

# calculate useful data summaries
getDataSummaries <-
function(obs)
{
  ## code runs at ~4,600 per sec
  # but its only called once
  out <-
  with(obs,
  list(
    Iages = as.numeric(rownames(index)),
    Iyrs = as.numeric(colnames(index)),
    Cages = as.numeric(rownames(catch)),
    Cyrs = as.numeric(colnames(catch))
  )) 
  within(out,
  {
    Iminage = min(Iages)
    Imaxage = max(Iages)
    Inyrs = length(Iyrs)
    Inages = length(Iages)
    Cminage = min(Cages)
    Cmaxage = max(Cages)
    Cnyrs = length(Cyrs)
    Cnages = length(Cages)
    Cyrsidx = 1:length(Cyrs)
    CinIages = Cages %in% Iages
    CinIyrs = Cyrs %in% Iyrs
  })
}

# good starting values... (NOT YET!)
getInitialParams <-
function(obs, M, control.model, data_summaries)
{
  ## code runs at ~40,000 per sec (because it doen't do anything yet!!)
  npars <- 2*data_summaries $ Cnyrs + 5
  if (control.model $ recruitment != "none") npars <- npars + 3
  
  # now do something clever
  with(data_summaries, rep(0, npars))
}


# the negative log likelihood
neg.llik <-
function(par, control.model, data_summaries, M, obs)
{

  pred <- predict_pop(par, control.model, data_summaries, M)
  
  # add smallest observed catch for zero catch
  obs $ catch[obs $ catch == 0] <- min(obs $ catch[obs $ catch > 0], na.rm = TRUE)
  
  llik.catch <- sum(dnorm(log(obs $ catch),  log(pred $ catch), exp(par[1]), log = TRUE), na.rm = TRUE)
  llik.index <- sum(dnorm(log(obs $ index),  log(pred $ index), exp(par[2]), log = TRUE), na.rm = TRUE)
   
  -1 * (llik.catch + llik.index)
}


# Predict population -- CONVERT TO USE Forward... 
#                       ONLY REQUIRED IF Recruit model is deterministic
predict_pop <-
function(par, control.model, data_summaries, M, full = FALSE)
{
  model_par <- getModelParams(par, control.model, data_summaries)
  
  Sfun <- getSelModel(control.model $ fishery.selectivity)
  qfun <- getSelModel(control.model $ survey.selectivity)
  Rfun <- getRecModel(control.model $ recruitment)
  
  # survey selectivity
  qa <- exp(model_par $ logq) * qfun(data_summaries $ Iages, model_par $ lambda, data_summaries)
  q <- qa %*% t(rep(1, data_summaries $ Inyrs))

  # F and Z
  S <- Sfun(data_summaries $ Cages, model_par $ aE, data_summaries)
  F <- S %*% t(model_par $ Fy)
  Z <- F + M
  
  # cohortify Z
  Zc <- 
    with(data_summaries,
      do.call(rbind, 
        lapply(1:Cnages, 
               function(i, nr) c(rep(NA,nr-i), Z[i,], rep(NA,i-1)), 
               nr = Cnages)))
  Zc <- Zc[,-(1:(data_summaries $ Cnages - 1))]      

  # fill out N with recruits
  Nc <- matrix(model_par $ Ry, data_summaries $ Cnages, data_summaries $ Cnyrs, byrow = TRUE)      
  # do cumulative Zs down cohorts
  Zc.cum <- rbind(0, Zc[-nrow(Zc),])
  for (i in 3:data_summaries $ Cnages) Zc.cum[i,] <- Zc.cum[i-1,] + Zc.cum[i,]
  
  # generate population
  Nc <- Nc * exp(-Zc.cum)

  # un-cohortify N
  N <-
    do.call(rbind,
      lapply(1:nrow(Nc), 
        function(i, nr) c(rep(NA, i-1), Nc[i,], rep(NA, nr-i+1)), 
        nr = nrow(Nc)))
  N <- N[,-(2:nrow(N) + ncol(F))]

  catch <- F / (F + M) * (1 - exp(-Z)) * N[,-ncol(N)]
  index <-  q * with(data_summaries, N[CinIages, c(CinIyrs, TRUE)])

  if (full)
    list(catch = catch, index = index, 
         N = N, F = F, S = S, qa = qa, 
         Sfun = Sfun, qfun = qfun, Rfun = Rfun, 
         par = model_par)
  else
    list(catch = catch, index = index)
}


# convert model parameters from vector to list 
# and convert supports
getModelParams <-
function(par, control.model, data_summaries)
{
  ## code runs at ~8,000 per sec
  with(data_summaries,
  {
    nsig <- 2
    out <-
    list(
      Fy     = exp(par[nsig + Cyrsidx]), 
	    Ry     = exp(par[nsig + Cnyrs+Cyrsidx]),
      logq   = par[nsig + 2*Cnyrs+1])
      
    # fishery selection model  
    out $ aE <- 
      getSelModelParams( par[nsig + 2*Cnyrs+2], 
                         control.model $ fishery.selectivity, 
                         data_summaries)
   
    # survey selection model      
    out $ lambda <- 
      getSelModelParams( par[nsig + 2*Cnyrs+3], 
                         control.model $ survey.selectivity,
                         data_summaries)
    
    # recruitment model (can be 'none')          
    out $ recruits <-
      getRecModelParams( par[nsig + 2*Cnyrs+4:5], 
                         control.model $ recruitment,
                         data_summaries)
                           
    out
  })
}

# convert params from real line to the desired domain
getSelModelParams <-
function(par, model, data_summaries)
{
  ## code runs at ~50,000 per sec 
  with(data_summaries,
    switch(model,
      linear = iglogit(par, min = Cminage, max = Cmaxage),
      # for the logit models theta is the age at 99% selection
      # so it makes no sense for theta to be able to reach 9
      logit1 = iglogit(par, min = Cminage, max = Cmaxage - 1e-9), 
      logit2 = iglogit(par, min = Cminage, max = Cmaxage - 1e-9),
      exponential = exp(par),
      stop("this shouldn't happen")
    ))
}    

      
# get a 1 parameter selection function
getSelModel <-
function(model)
{
  ## code runs at ~1000,000 per sec
  switch(model,
    linear = 
      function(a, theta, data_summaries) 
        with(data_summaries, 
          ifelse(a < theta, 0.05 + 0.95 * (a - Cminage)/(theta - Cminage), 1)),
    # for the logit models theta is the age at 99% selection
    # so it makes no sense for theta to be able to reach 9
    logit1 = 
      function(a, theta, data_summaries) 
        with(data_summaries,
          iglogit( glogit(a, Cminage, Cmaxage) + 
                     log(99) - glogit(theta, Cminage, Cmaxage), 
                   0.05, 1)),
    logit2 = 
      function(a, theta, data_summaries) 
        with(data_summaries, 
          iglogit( (exp(glogit(theta, Cminage, Cmaxage)) + .5) * 
                     glogit(a, Cminage, Cmaxage) + log(99) - 
                     glogit(Cmaxage - .5, Cminage, Cmaxage) * 
                     (exp(glogit(theta, Cminage, Cmaxage)) + .5),
                   0.05, 1)),
    exponential = 
      function(a, theta, data_summaries)
        exp(-theta * a),
    stop("this shouldn't happen")
  )      
}

# convert params from real line to the desired domain
getRecModelParams <-
function(par, model, data_summaries)
{
  ## code runs at ~XX per sec 
  with(data_summaries,
    switch(model,
      none = numeric(0),
      ricker = par,
      stop("this shouldn't happen")
    ))
}    

      
# get a 2 parameter recruitment function ...
getRecModel <-
function(model)
{
  ## code runs at ~1000,000 per sec
  switch(model,
    none = NULL, 
    ricker = 
      function(ssb, theta, data_summaries)
        rep( theta[1], length(ssb)), # constant recruitment
    stop("this shouldn't happen")
  )      
}

