# see README
require(rcolgem) 

# MODEL ORGANISATION AND DEMES
UMDEMES <- paste( 'UM', 0:4 , sep = '')
#~ UMDEMES <- c('Em', 'Cm', 'Am')
DMDEMES <- paste(sep='', 'DM', 0:4) 
TMDEMES <- paste(sep='', 'TM', 0:4) 

UFDEMES <- paste( 'UF', 0:4 , sep = '')
#~ UFDEMES <- c('Ef', 'Cf', 'Af')
DFDEMES <- paste(sep='', 'DF', 0:4) 
TFDEMES <- paste(sep='', 'TF', 0:4) 

UMSMDEMES <- paste( 'Umsm', 0:4 , sep = '')
#~ UMSMDEMES <- c('Emsm', 'Cmsm', 'Amsm')
DMSMDEMES <- paste(sep='', 'Dmsm', 0:4) 
TMSMDEMES <- paste(sep='', 'Tmsm', 0:4) 

MDEMES <- c( UMDEMES, DMDEMES, TMDEMES )
FDEMES <- c( UFDEMES, DFDEMES, TFDEMES )
MSMDEMES <- c( UMSMDEMES, DMSMDEMES, TMSMDEMES )

E_DEMES <- c(UMDEMES[1], DMDEMES[1], TMDEMES[1], UFDEMES[1], DFDEMES[1], TFDEMES[1], UMSMDEMES[1], DMSMDEMES[1], TMSMDEMES[1] )
C_DEMES <- c(UMDEMES[2:4], DMDEMES[2:4], TMDEMES[2:4], UFDEMES[2:4], DFDEMES[2:4], TFDEMES[2:4], UMSMDEMES[2:4], DMSMDEMES[2:4], TMSMDEMES[2:4] )
#~ A_DEMES <- c(UMDEMES[3], DMDEMES[3], TMDEMES[3], UFDEMES[3], DFDEMES[3], TFDEMES[3], UMSMDEMES[3], DMSMDEMES[3], TMSMDEMES[3] )
A_DEMES <- c(UMDEMES[5], DMDEMES[5], TMDEMES[5], UFDEMES[5], DFDEMES[5], TFDEMES[5], UMSMDEMES[5], DMSMDEMES[5], TMSMDEMES[5] )
AnT_DEMES <- c(UMDEMES[5], DMDEMES[5], UFDEMES[5], DFDEMES[5], UMSMDEMES[5], DMSMDEMES[5])
AT_DEMES <- c(TMDEMES[5], TFDEMES[5],  TMSMDEMES[5] )

U_DEMES <- c( UMDEMES, UFDEMES, UMSMDEMES )
D_DEMES <- c( DMDEMES, DFDEMES, DMSMDEMES )
T_DEMES <- c( TMDEMES, TFDEMES, TMSMDEMES )

DEMENAMES <- c( UMDEMES, DMDEMES, TMDEMES
  , UFDEMES, DFDEMES, TFDEMES
  , UMSMDEMES, DMSMDEMES, TMSMDEMES
  ,   'Source') #28 demes 
NONSOURCEDEMES <- DEMENAMES[-length(DEMENAMES) ]
STATENAMES <- c( DEMENAMES[-length(DEMENAMES)],  'Sm', 'Sf', 'Smsm')
m <- length(DEMENAMES)
# /MODEL ORGANISATION

# DATA
TREES <- read.tree( 'mergeTrees0.nwk' )
ssdf  <- read.csv( 'sampleStates0.csv')
ss    <- ssdf[,2:ncol(ssdf)]
rownames(ss) <- ssdf[,1]
ss <- as.matrix ( ss )
st <- setNames( sapply( strsplit( rownames(ss), split='_'), function(splits) as.numeric(splits[length(splits)] ) ) , rownames(ss))
model5sampleStates.to.model6sampleStates <- function(ss) 
{
	ss1 <- matrix(0, nrow=nrow(ss), ncol = length(DEMENAMES ) )
	colnames(ss1) <- DEMENAMES
	rownames(ss1) <- rownames(ss) 
	wgamma <- 	c( gamma1 = 1/0.157 
	, gamma2 = 1/0.350  
	, gamma3 = 1/0.282 )
	wgamma <- unname(wgamma) / sum(wgamma )
	for (i in 1:nrow(ss))
	{
		ss1[i, 'Source'] <- ss[i, 'Source']
		ss1[i, c('UM0', 'UF0', 'Umsm0')] <- unname( ss[i, c('Em', 'Ef', 'Emsm')] )
		ss1[i, c('UM4', 'UF4', 'Umsm4')] <- unname( ss[i, c('Am', 'Af', 'Amsm')] )
		ss1[i,  paste( 'UM', 1:3, sep='' ) ] <- unname( ss[i, 'Cm'] * wgamma )
		ss1[i,  paste( 'UF', 1:3, sep='' ) ] <- unname( ss[i, 'Cf'] * wgamma )
		ss1[i,  paste( 'Umsm', 1:3, sep='' ) ] <- unname( ss[i, 'Cmsm'] * wgamma )
	}
	ss1
}
sampleStates <- model5sampleStates.to.model6sampleStates( ss )
BDTS <- lapply( TREES, function(tre)
{
	binaryDatedTree( tre, st , sampleStates )
})
#/DATA 




#TIMES
#~  dates relative to birthday
BIRTHDAY <- as.Date('1980-11-11')
TSTART <- as.Date('1930-1-1')  
MAXSAMPLETIME <- 12395
MAXHEIGHT <- as.numeric( BIRTHDAY + MAXSAMPLETIME - TSTART ) # corresponds to 1930 
DELTAT <- 31 #61 # days
TIMES <- seq(MAXSAMPLETIME - MAXHEIGHT, MAXSAMPLETIME+DELTAT, by=DELTAT) #note first time point is negative, corresponds to 1930

#~ functions to map day time axis relative to tstart onto calendar years; 
tdays2calyears <- function(tdays)
{
	1980 + as.numeric( BIRTHDAY - as.Date('1980-1-1') ) / 365  + tdays / 365 
}
calyrs2days <- function(yr)
{
	dBIRTHDAY <- 0
	yrBIRTHDAY <- 1980 + as.numeric( BIRTHDAY - as.Date('1980-1-1') ) / 365
	(yr - yrBIRTHDAY ) * 365  + dBIRTHDAY
}
TIMES_CALYEARS <-  tdays2calyears( TIMES )
EXPANSIONTIME <-  as.numeric( as.Date('1975-1-1') - BIRTHDAY ) # 1975 
# /TIMES




#GLOBAL PARMS
NCORES <- 20
SOURCESIZE <- 1e7
# default MSM parameters  
# 4.310345e-05 #1.436782e-05 # proportion of men MSM --- stef 1pc  --- considering just abuja 3e6*1e-2 / 2 = 15k --- 15k / 174e6/2 
#~  DHS: 25476 total msm in 2012, all nigeria; 
#~ > 25476 / 45.24e6 # dhs report / adult male nigeria 
#~ [1] 0.00056313
#~ > 25476 / 16.5e6 # dhs report / urban adult male nigeria
#~ [1] 0.001544
#~ PMSM <- .005 
MAX_PMSM <- .01 
ASSORTPROB <- 0.90 #MSM vs general 
RETENTION_PROB <- .78 # retention on art within 1 year; .78 = exp(-r*365 ) ; 
T2D_RATE <-  -log(RETENTION_PROB) / 365  # 
MU <- 1 / (365*(58-13)) # natural mortality, approximately life expectancy minus sexual debut
#~ Population of 'Federal Capital Territory' (2006 census)
#~  Total	1,405,201
#~ N1 <- 3e6 # 175*1e6 # nigeria
rdemographic <- 0.02527781 / 365  # growth rate 2.56% per year 
N1 <- 1405201 * exp(rdemographic * (2014-2007) *365 ) # note adjustment from 2007 census estimate 
N_NIGERIA <- 175e6 # 2014
PROP_ABUJA <- N1 / N_NIGERIA
N.t <- function(t) { #exponential growth before 1970
	N1 * exp(-rdemographic * (MAXSAMPLETIME - t))
}
TRE_BREAK <- .25 # how much disease progression slowed by treatment
# parameter names and boundaries 
EPSTLB <- .01 #applied in simulate 
EPSTUB <- .27 
BETA    <- .85 / 365 #1 * .58035640 / 365
BETAMSM <- 1 * BETA 
RATIOPOPSIZENEWOLD <-  .5 # ? calibrated to give realistic prevalence # 3e6 / 174e6
DEFAULTPARAMETERS <- c( y0 = .06973729
	 , beta = BETA
	 ,  w1 = 1/26.04, w2 = 7/ 26.04
	 , eps = .75 
	 , epsilonT = .5*(EPSTUB -EPSTLB )
	 , epsilonD = .5
	 , immigrationRate = 1 / 10 / 365
	 , pmsm = .00154 
	 , assortprob = ASSORTPROB
	 , assortprob_msm = ASSORTPROB
	 , betamsm = BETAMSM
	 , gamma0 = 1/365 # EHI
	, gamma1 = 1/ ( (1-.75/12) /(0.157 / 365 )  )
	, gamma2 = 1/ ( (1-.75/12) /(0.350 / 365 )  )
	, gamma3 = 1/ ( (1-.75/12) /(0.282 / 365 )  )
	, gamma4 = 1/ ( (1-.75/12) /(.434 / 365 )  ) #aids 
	, pstartstage1 = 0.58
	, pstartstage2 = 0.23  #reg4 combine these 
	, pstartstage3 = 0.16
	, pstartstage4 = 0.03
	, treBreak = TRE_BREAK
) 
#/GLOBAL PARMS



source.size <- function(t, parms)
{ #model 5 - expon growth of source
#~ note that sourcesize is fixed, w/ will reduce accuracy of this approximation 
	unname( SOURCESIZE * exp(- parms$sourceGrowthRate * (MAXSAMPLETIME - t)) )
}





# force of infection scaling 
################
# scale transmission rates with generalised logistic function of time; f(1980 = 1): 
f.t <- function(t, ns_parmlist )
{ 
	tt <- tdays2calyears(t) 
	if (tt < 1980) return(1) 
	with ( ns_parmlist , 
		min( 1, 1 + ( 1 + exp(-fB*(1980-fM) ) )^(-1/fV)  - ( 1 + exp(-fB*(tt-fM) ) )^(-1/fV))
	)
}
f.t.msm <- function(t, ns_parmlist )
{ 
	tt <- tdays2calyears( t )
	#if (tt < 1980) return(f.t.msm( , ns_parmlist ) ) 
	with ( ns_parmlist , 
		min( 1, ( 1 + exp(-fmsmB*(tt-fmsmM) ) )^(-1)  ) 
	)
}

################
diagnosis.rate <<- function(t, ns_parmlist) 
{ # assumes this applies to all undiagnosed demes equally 
	tt <- tdays2calyears( t )
	with ( ns_parmlist , 
		dK*( 1 + exp(-dB*(tt-dM) ) )^(-1)
	)
}
diagnosis.rate.msm <<- diagnosis.rate

TINTERVENTIONSTART <- Inf


treatment.rate.vector <<- function(t, ns_parmlist)
{ # treatment rates for each stage of infection (5-vector )
	tt <- tdays2calyears( t) 
	if (tt >= TINTERVENTIONSTART )  { # universal treatment
		return( rep( ns_parmlist$treRateIntervention, 5 ) ) 
	} else if (tt > 2011) #incl period after intervention
	{ #cd4 < 350
		return( c( 0,0,0, rep(ns_parmlist$treRate, 2 ) ) )
	} else if (tt > 2004 & tt <= 2011)
	{ #cd4 < 200 only 
		return( c( 0,0,0,0, ns_parmlist$treRate ) )
	} else{
		return ( rep(0, 5))
	}
}
treatment.rate.vector.msm <<- treatment.rate.vector

non.adherence.rate <<- function(t) T2D_RATE
non.adherence.rate.msm <<- function(t) T2D_RATE

.F.from.state <- function(t, y, parms)
{
	FF <- matrix(0, nrow=m, ncol=m)
	colnames(FF) = rownames(FF) <- DEMENAMES
	
	Im <- sum(y[c(MDEMES)])
	Imsm <- sum(y[c(MSMDEMES)])
	If <- sum(y[FDEMES])
	
	y <- setNames( pmax(0, y), names(y) )
	with(as.list(y), 
	{
		Nf  <- Sf + If
		Nm  <- Sm + Im 
		Nmsm <- Smsm + Imsm
		N   <- Nf + Nm + Nmsm
		I <- If + Im + Imsm
		
		foi <- f.t( t, parms ) * parms$beta
		foi_msm <- f.t.msm(t, parms ) * parms$betamsm
		
		w <- rep( c(1, rep(parms$w1,3), parms$w2 ), 3 )
		
		f_m2f <-   foi * w * y[MDEMES]  * (Sf  / Nf)  
		
		f_f2m <-  parms$eps * foi * w * y[FDEMES] * parms$assortprob * Sm  / Nm
		f_f2msm <-  parms$eps * foi * w * y[FDEMES] * (1-parms$assortprob) *  Smsm  / Nmsm
		
		f_msm2f   <- foi_msm * w * y[MSMDEMES] * (1-parms$assortprob_msm) * Sf / Nf
		f_msm2msm <- foi_msm * w * y[MSMDEMES] * parms$assortprob_msm * Smsm / Nmsm
		
		FF[MDEMES, 'UF0'] <- f_m2f
		
		FF[FDEMES, 'UM0'] <-   f_f2m
		FF[FDEMES, 'Umsm0'] <-  f_f2msm
		
		FF[MSMDEMES, 'Umsm0'] <- f_msm2msm
		FF[MSMDEMES, 'UF0']   <- f_msm2f
		
		# treatment
		FF[T_DEMES,] <- FF[T_DEMES,] * parms$epsilonT
		# diagnosis
		FF[D_DEMES,] <- FF[D_DEMES,] * parms$epsilonD
		
		FF['Source', 'Source'] <-  source.size(t, parms) 
#~ if (t > 12e3) browser()
		FF
	})
}


.G.from.state <- function(t, y, parms)
{
	GG <- matrix(0, nrow=m, ncol=m)
	colnames(GG) = rownames(GG) <- DEMENAMES
	
	p <- unname( unlist(parms[ c(paste(sep='', 'pstartstage', 1:4)) ] )  )
	gammas13 <-  unname( unlist(  parms[c('gamma1', 'gamma2', 'gamma3') ]  ))
	
		#stage progression 0
		GG[ 'UM0', paste(sep='', 'UM', 1:4) ] <- y['UM0'] * parms$gamma0 * p 
		GG[ 'DM0', paste(sep='', 'DM', 1:4) ] <- y['DM0'] * parms$gamma0 * p 
		GG[ 'TM0', paste(sep='', 'TM', 1:4) ] <- y['TM0'] * parms$gamma0 * p 
		
		GG[ 'UF0', paste(sep='', 'UF', 1:4) ] <- y['UF0'] * parms$gamma0 * p 
		GG[ 'DF0', paste(sep='', 'DF', 1:4) ] <- y['DF0'] * parms$gamma0 * p 
		GG[ 'TF0', paste(sep='', 'TF', 1:4) ] <- y['TF0'] * parms$gamma0 * p 
		
		GG[ 'Umsm0', paste(sep='', 'Umsm', 1:4) ] <- y['Umsm0'] * parms$gamma0 * p 
		GG[ 'Dmsm0', paste(sep='', 'Dmsm', 1:4) ] <- y['Dmsm0'] * parms$gamma0 * p 
		GG[ 'Tmsm0', paste(sep='', 'Tmsm', 1:4) ] <- y['Tmsm0'] * parms$gamma0 * p 
		
		# stage progression 1-3
		GG[ cbind( UMDEMES[c(-1,-5)], UMDEMES[c(-1,-2)] ) ] <- gammas13 * y[ UMDEMES[c(-1,-5)] ]
		GG[ cbind(DMDEMES[c(-1,-5)], DMDEMES[c(-1,-2)] ) ] <- gammas13 * y[ DMDEMES[c(-1,-5)] ]
		GG[ cbind( TMDEMES[c(-1,-5)], TMDEMES[c(-1,-2)]  )] <- gammas13 * y[ TMDEMES[c(-1,-5)] ] * parms$treBreak

		GG[ cbind( UFDEMES[c(-1,-5)], UFDEMES[c(-1,-2)] )] <- gammas13 * y[ UFDEMES[c(-1,-5)] ]
		GG[ cbind( DFDEMES[c(-1,-5)], DFDEMES[c(-1,-2)] ) ] <- gammas13 * y[ DFDEMES[c(-1,-5)] ]
		GG[ cbind( TFDEMES[c(-1,-5)], TFDEMES[c(-1,-2)] ) ] <- gammas13 * y[ TFDEMES[c(-1,-5)] ] * parms$treBreak
		
		GG[ cbind( UMSMDEMES[c(-1,-5)], UMSMDEMES[c(-1,-2)] )] <- gammas13 * y[ UMSMDEMES[c(-1,-5)] ]
		GG[ cbind( DMSMDEMES[c(-1,-5)], DMSMDEMES[c(-1,-2)] ) ] <- gammas13 * y[ DMSMDEMES[c(-1,-5)] ]
		GG[ cbind( TMSMDEMES[c(-1,-5)], TMSMDEMES[c(-1,-2)] ) ] <- gammas13 * y[ TMSMDEMES[c(-1,-5)] ] * parms$treBreak
		
		# diagnosis
		dr <- diagnosis.rate(t, parms )
		dr_msm <- diagnosis.rate.msm(t, parms )
		GG[cbind(UMDEMES, DMDEMES)] <-  dr * y[UMDEMES] 
		GG[cbind(UFDEMES, DFDEMES)] <-  dr * y[UFDEMES] 
		GG[cbind(UMSMDEMES, DMSMDEMES)] <-  dr_msm * y[UMSMDEMES] 
		
		# treatment
		trv <- treatment.rate.vector(t, parms )
		trv_msm <- treatment.rate.vector.msm(t, parms )
		#trv <- rep(trv, 3 )
		GG[ cbind( DMDEMES, TMDEMES) ] <- trv * y[DMDEMES] 
		GG[ cbind( DFDEMES, TFDEMES) ] <- trv * y[DFDEMES] 
		GG[ cbind( DMSMDEMES, TMSMDEMES) ] <- trv_msm * y[DMSMDEMES] 
		
		# non adherence 
		t2drate <- non.adherence.rate( t )
		t2drate_msm <- non.adherence.rate.msm( t )
		GG[ cbind( TMDEMES, DMDEMES) ] <- t2drate * y[TMDEMES] 
		GG[ cbind( TFDEMES, DFDEMES) ] <- t2drate * y[TFDEMES] 
		GG[ cbind( TMSMDEMES, DMSMDEMES) ] <- t2drate_msm * y[TMSMDEMES] 
		
		GG['Source',] <- y[DEMENAMES] * parms$immigrationRate
		if (t < EXPANSIONTIME) GG['Source',] <-  GG['Source',] * 10
		GG['Source', 'Source'] <- 0
#~ if (t > 12e3) browser()
		GG
}
.Y.from.state <- function(t, y, parms)
{
	yy <- c(y, Source=source.size(t, parms )) 
	yy <- setNames( pmax(0, yy) , names(yy) )
	yy[DEMENAMES]
}


dy <- function(t, y, parms, ...)
{
	y <- setNames( pmax(0, y), names(y) )
	Im <- sum(y[c(MDEMES)])
	Imsm <- sum(y[c(MSMDEMES)])
	If <- sum(y[FDEMES])
	interventionEffect <- 1 #  
	if (t < EXPANSIONTIME) return(list(  setNames( rep(0, length(y)) , STATENAMES  )) )
	FF <- .F.from.state(t, y, parms)
	GG <- .G.from.state(t, y, parms)
	csFF <- colSums(FF)
	rsGG <- rowSums(GG)
	csGG <- colSums(GG) 
	names(csFF) = names(rsGG) = names(csGG)  <- rownames(FF) 
	d <-  setNames( rep(0, length(y) ), STATENAMES) 
	d[NONSOURCEDEMES]  <- csFF[NONSOURCEDEMES] + csGG[NONSOURCEDEMES] - rsGG[NONSOURCEDEMES]
	#with(as.list(y), 
	{
		Nf  <- y['Sf'] + If
		Nm  <- y['Sm'] + Im 
		Nmsm <- y['Smsm'] + Imsm
		N   <- Nf + Nm + Nmsm
		b   <- log( N.t(t+1) / N.t(t) ) +
		  ( parms$gamma4 * sum(y[AnT_DEMES]) + parms$gamma4 * parms$treBreak * sum(y[AT_DEMES]) ) /N + 
		  MU 
		d['Sm'] <- b * (1-parms$pmsm) * N / 2 - sum(csFF[c(MDEMES)])   - MU*y['Sm']
		d['Smsm'] <- b * (parms$pmsm) * N / 2 - sum(csFF[c(MSMDEMES)]) - MU*y['Smsm']
		d['Sf'] <- b * N / 2 - sum(csFF[c(FDEMES)]) - MU*y['Sf']
	}#)
	
	d[NONSOURCEDEMES] <- d[NONSOURCEDEMES] - MU * y[NONSOURCEDEMES]
	
	# aids deaths 
	d[AnT_DEMES] <- d[AnT_DEMES] - parms$gamma4 * y[AnT_DEMES] 
	d[AT_DEMES]  <- d[AT_DEMES]  - parms$gamma4 * parms$treBreak * y[AT_DEMES] 
	list(d)
	
}


simulate <- function(theta)
{ # theta : vector 
	params <- DEFAULTPARAMETERS
	params[names(theta)] <- theta 
	params['pmsm'] <- min( MAX_PMSM, params['pmsm'] )
	
	I0 <- params['y0'] #* (MAXY-MINY) + MINY
	S0 <- N.t(MAXSAMPLETIME - MAXHEIGHT)  #* exp(-rdemographic*(MAXHEIGHT))
	times <- seq(MAXSAMPLETIME - MAXHEIGHT, MAXSAMPLETIME+DELTAT, by=DELTAT)
	y0 <- setNames( rep(0, length(STATENAMES)), STATENAMES)
	y0['Sf'] <- S0 / 2
	y0['Smsm'] <- params['pmsm'] * S0/2
	y0['Sm'] <- (1 - params['pmsm']) * S0/2
	y0[UFDEMES] = y0[UMDEMES] = y0[UMSMDEMES]  <- I0 #/ 1e2 #
	#y0['UM0'] = y0['UF0'] <- I0
	odeparms <- as.list(params)
	odeparms$t0 <- MAXSAMPLETIME - MAXHEIGHT
	
	# to speed things up, will use rk4 and solve ODEs only after expansiontime
	times2 <- times[ times > EXPANSIONTIME] 
	y0['Sm'] = y0['Sf'] <-  N.t(EXPANSIONTIME) / 2
	
	o0 <- cbind( times[ times <= EXPANSIONTIME],  matrix(y0,  byrow = TRUE, nrow = length(times) - length(times2), ncol = length(STATENAMES)) )
	o1 <- ode( y=y0, times=times2, func = dy, parms = odeparms, method='rk4')
	o <- rbind(o0, o1)
	#class( o) <- 'deSolve'
	colnames(o)[2:ncol(o)] <- STATENAMES
	oo <- o[,2:ncol(o)]
	Fs <- lapply( 1:nrow(o), function(k) .F.from.state(o[k,1], oo[k,], odeparms) ) #\times DELTAT  ? ?  No. 
	Gs <- lapply( 1:nrow(o), function(k) .G.from.state(o[k,1], oo[k,], odeparms) )
	Ys <- lapply( 1:nrow(o), function(k) .Y.from.state(o[k,1], oo[k,], odeparms) )
	list(times, Fs, Gs, Ys, o[,'Sm']+o[,'Sf']+o[,'Smsm'], o)
}



LOG_SCALE_PARAMETERS <- c( 'y0', 'immigrationRate', 'sourceGrowthRate'
  ,  'beta',  'fB', 'fM', 'fV'
  , 'betamsm',  'fmsmM' #'fmsmB'can be negative
  , 'dK', 'dB', 'dM'
  , 'treRate'
) #
LOGISTIC_SCALE_PARAMETERS <- c(  'assortprob', 'assortprob_msm', 'w1', 'w2', 'pmsm')
.natural2transformed.scale  = natscale2transformed <- function(params)
{
	logsp <- names(params)[ names(params) %in% LOG_SCALE_PARAMETERS ] 
	lcsp <- names(params)[ names(params) %in% LOGISTIC_SCALE_PARAMETERS ]
	params[logsp] <- log(params[logsp])
	params[lcsp] <- log(params[lcsp]/ (1-params[lcsp]) )
	params
}

.transformed2natural.scale = transformed2natscale <- function(params)
{
	logsp <- names(params)[ names(params) %in% LOG_SCALE_PARAMETERS ] 
	lcsp <- names(params)[ names(params) %in% LOGISTIC_SCALE_PARAMETERS ]
	params[logsp] <- exp(params[logsp])
	params[lcsp] <- 1 / (1 + exp(-params[lcsp]) )
	params
}





#~ INFERENCE
rprior = rprior.vec <- function()
{ #
	c( assortprob =  rbeta(1,  5,1 )
	  , assortprob_msm = rbeta(1,  5,1 )
	  , w1 =  rbeta(1, 1, 2 ) 
	  , w2 = rbeta( 1, 1, 2 )
	  , y0 = runif(1, 0, 10 )
	  , immigrationRate = 1 / (365 * runif( 1, 5, 30) )
	  , sourceGrowthRate = 1 / (365 * runif(1, 10/12, 3.5) ) #1/x prior based on empirical observations of doubling time
	  
	  , pmsm = rlnorm( 1,  log(.001544), 1)
	  	
	  	, beta = rlnorm(1,  log( BETA ) , 1/2) 
		, fB = rlnorm( 1, log(.25), 1 )
		, fM = rnorm(1,  mean = 1990, sd = 5)
		, fV = rlnorm(1, log(1), .5)
		
		, betamsm = rlnorm(1,  log( BETAMSM ) , abs(log(2* .58))) ## good
		, fmsmB = rnorm(1, 0, sd = .25 ) #rlnorm( 1, log(.25), 1 )
		, fmsmM = rnorm(1,  mean = 1995, sd = 5)

		, dK =  rlnorm( 1, log(1/10/365), .6 ) #median 4yrs
		, dB = rlnorm( 1, log(.5), .4  )  
		, dM = rnorm(1, 2006, sd= 1) 
		
		, treRate = rlnorm( 1, log(1/(2*365)), .4 ) 
	)
}


dprior <- function( theta, sim = NA, log = TRUE)
{
	d <- with(as.list(theta), {
		ps <- c( dbeta(assortprob, 5,1, log = log)
		  , dbeta(assortprob_msm, 5,1, log = log)
		  , dbeta( w1, 1, 2, log = log)
		  , dbeta( w2, 1, 2 , log = log) 
		  , dunif( y0, 0, 10, log = log)
		  , dunif(  1 / immigrationRate, 5*365, 30*365, log = log )
		  , dunif(  1 / sourceGrowthRate, 365*10/12, 365*3.5, log = log )
		  , dlnorm( pmsm, log(.001544), 1, log = log )
		  
			, beta = dlnorm(beta,  log( BETA ) , 1/2, log = log ) 
			, fB = dlnorm( fB, log(.25), 1, log = log  ) #dunif(fB, 0, 1, log = log)
			, fM = dnorm(fM,  mean = 1990, sd = 5, log = log)
			, fV = dlnorm(fV, log(1), .5, log = log)
			
			, betamsm = dlnorm(betamsm,  log( BETAMSM ) , abs(log(2* .58)), log = log) 
			, fmsmB = dnorm(fmsmB, 0, sd = .25, log = log  )
			, fmsmM = dnorm(fM,  mean = 1990, sd = 5, log = log)
			#, fmsmV = dlnorm(fV, log(1), .5, log = log)
			
			, dK = dlnorm( dK, log(1/10/365), .6, log = log  ) # dunif( dK, 0, 1/8, log = log)
			, dB = dlnorm( dB, log(.5), .4 , log = log ) 
			, dM = dnorm(dM, 2006, sd= 1, log = log)
			
			, treRate = dlnorm( treRate, log(1/(2*365)), .4 , log = log ) 
		)
		ifelse( log, sum(ps), prod(ps) )
	})
	if ( length(sim) > 1 )
	{
		i2001 <- which.min( abs(TIMES_CALYEARS - 2001)  )
		i2008 <- which.min( abs(TIMES_CALYEARS - 2008)  ) 
		i2012 <- which.min( abs(TIMES_CALYEARS - 2012)  )
		
		I2001 <- sum( sim[[4]][[i2001]][NONSOURCEDEMES] ) 
		I2012 <- sum( sim[[4]][[i2012]][NONSOURCEDEMES] ) 
		
		p2001 <- (I2001) / ( sim[[5]][ i2001 ] +  (I2001) )
		p2012 <- (I2012) / ( sim[[5]][ i2012 ] +  (I2012) )
		
		# unaids dhs report 
		d <- ifelse(log
		  ,  d + dnorm( p2001, mean=.035, sd=.01 /(2*1.96), log = log ) + dnorm( p2012, mean=.031, sd=.007 /(2*1.96), log = log )
		  ,  d *dnorm( p2001, mean=.035, sd=.01 /(2*1.96), log = log ) * dnorm( p2012, mean=.031, sd=.007 /(2*1.96), log = log )
		)
		
		#number on treatment
		#http://www.unaids.org/sites/default/files/media_asset/20130630_treatment_report_en_0.pdf
		Ntreated2012 <-  491021 # unaids 2012 
		Ntreat_abuja2012 <- PROP_ABUJA * Ntreated2012
		Ntreat <- sum( sim[[4]][[i2012]][T_DEMES] )  
		# note sd used here is arbitrary and small; simply used for smoothly docking sims to observed values
		d <- ifelse(log
		  ,  d + dnorm( Ntreat, mean=Ntreat_abuja2012, sd=500, log = log )
		  ,  d * dnorm( Ntreat, mean=Ntreat_abuja2012, sd=500, log = log ) 
		)
		
		# testing rates
		G2008 <- sim[[3]][[i2008]]
		Y2008 <- sim[[4]][[i2008]]
		tr <- 365 * sum( G2008[U_DEMES, D_DEMES] ) / sum(Y2008[U_DEMES] )
		# note sd used here is arbitrary and small; simply used for smoothly docking sims to observed values
		d <- ifelse(log
		  ,  d + dnorm( tr, mean=.066, sd=(.04) /(2*1.96), log = log )
		  ,  d * dnorm( tr, mean=.066, sd=(.04) /(2*1.96), log = log ) 
		)
		
		#MSM treatment 
		#1.       For those already diagnosed at enrollment, ~30-32% were on ART. (source : bridging trust)
		TELIG_DMSM_DEMES <- c(DMSMDEMES[4:5] ) 
		Nelig_msm <- sum( sim[[4]][[i2012]][TELIG_DMSM_DEMES] ) 
		Ntreat_msm <- sum( sim[[4]][[i2012]][TMSMDEMES] ) 
		ptreat_msm <- Ntreat_msm / (Ntreat_msm + Nelig_msm )
		d <- ifelse(log
		  ,  d + dnorm( ptreat_msm, mean=.31, sd=(.40-.23) /(2*1.96), log = log )
		  ,  d * dnorm( ptreat_msm, mean=.31, sd=(.40-.23) /(2*1.96), log = log ) 
		)
		
		#MSM prop diagnosed 
		#(source : bridging trust)
		Ndiag_msm <- sum( sim[[4]][[i2012]][c(DMSMDEMES,TMSMDEMES)] )  
		pdiag_msm <- Ndiag_msm / sum( sim[[4]][[i2012]][MSMDEMES] )  
		p <- .33 
		n <- 220 
		s <- sqrt(n*p*(1-p) / n^2)
		d <- ifelse(log
		  ,  d + dnorm( pdiag_msm, mean=p, sd=s, log = log )
		  ,  d * dnorm( pdiag_msm, mean=p, sd=s, log = log ) 
		)
	}
	d
}



logx2meanlogx <- function(x)
{
	log( mean( exp( x-max(x) ) ) ) + max(x)
}

co.log.lik_sim.dompi <- function( theta  ) 
{
	tfgy <- simulate(theta) 
	clls <- foreach( bdt = iter(BDTS), .inorder=F) %dopar% 
	{ #.export=ls(envir=globalenv()) 
		cll <-  coalescent.log.likelihood.fgy(  bdt
		 , times=tfgy[[1]]
		 , births=tfgy[[2]]
		 , migrations=tfgy[[3]]
		 , demeSizes=tfgy[[4]]
		 , integrationMethod = 'rk4'
		 , censorAtHeight=30 * 365
		 , forgiveAgtY=.2
		 , returnTree=FALSE
		)
		cll
	}
	cll <- logx2meanlogx( unlist( clls )  )
	print( c( date(), theta, cll) )
	list( cll, tfgy )
}




# DIMENSION REDUCTION
DEMEMAP0 <- list(
	UM0='UM0'
	, UM13 = paste(sep='', 'UM', 1:3)
	, UM4 = 'UM4'
	, UF0='UF0'
	, UF13 = paste(sep='', 'UF', 1:3)
	, UF4 = 'UF4'
	, Umsm0='Umsm0'
	, Umsm13 = paste(sep='', 'Umsm', 1:3)
	, Umsm4 = 'Umsm4'
	
	, DM0= c('DM0', 'TM0')
	, DM13 = c(paste(sep='', 'DM', 1:3), paste(sep='', 'TM', 1:3) )
	, DM4 = c('DM4', 'TM4')
	
	, DF0= c('DF0', 'TF0')
	, DF13 = c(paste(sep='', 'DF', 1:3), paste(sep='', 'TF', 1:3) )
	, DF4 = c('DF4', 'TF4')
	
	, Dmsm0= c('Dmsm0', 'Tmsm0')
	, Dmsm13 = c(paste(sep='', 'Dmsm', 1:3), paste(sep='', 'Tmsm', 1:3) )
	, Dmsm4 = c('Dmsm4', 'Tmsm4')
	
	, Source = 'Source'
)
DEMES0 <- names(DEMEMAP0)
dim.reduce.matrix0 <- function(M) 
{
	MM0 <- matrix(0, nrow = length(DEMES0), ncol=ncol(M))
	colnames(MM0) <- colnames(M)
	rownames(MM0) <- DEMES0
	MM <- matrix(0, nrow = length(DEMES0), ncol=length(DEMES0))
	colnames(MM)=rownames(MM) <- DEMES0
	for (n in (DEMES0))
	{
		if (length( DEMEMAP0[[n]] ) > 1 )
		{
			MM0[n, ] <- colSums(M[ DEMEMAP0[[n]] , ] ) 
		} else
		{
			MM0[n,] <- M[ DEMEMAP0[[n]],  ] 
		}
	}
	for (n in (DEMES0))
	{
		if (length( DEMEMAP0[[n]] ) > 1 )
		{
			MM[, n] <- rowSums(MM0[  ,DEMEMAP0[[n]] ] ) 
		}else
		{
			MM[, n] <- MM0[, DEMEMAP0[[n]] ] 
		}
	}
	MM
}
dim.reduce.vector0 <- function(V) 
{
	setNames( sapply( DEMES0, function(n)
	{
		sum( V[DEMEMAP0[[n]] ] ) 
	}), DEMES0)
}
dim.reduce.sim0 <- function(sim)
{
# combine chronic stages and (D,T)
	list(
		sim[[1]]
		, lapply( sim[[2]], dim.reduce.matrix0)
		, lapply( sim[[3]], function(M) { MM <- dim.reduce.matrix0(M); diag(MM)<- 0; MM } ) 
		, lapply( sim[[4]], dim.reduce.vector0)
		, sim[[5]]
		, sim[[6]]
	)
}
sampleStates_dr0 <- t(sapply( 1:nrow(sampleStates), function(i)
{
	dim.reduce.vector0(sampleStates[i, ] )
}))
rownames( sampleStates_dr0 ) <- rownames( sampleStates )
colnames(sampleStates_dr0) <- DEMES0
BDTS_dr0  <- lapply( TREES, function(tre)
{
	binaryDatedTree( tre, st , sampleStates_dr0 )
})
# /DIMENSION REDUCTION


co.log.lik_sim.serial <- function( theta  ) 
{
	tfgy <- simulate(theta) 
	require(foreach)
	require(iterators)
	{
		sim <- tfgy 
		#plot( sim[[6]] ) 
		I <- matrix( unlist(sim[[4]]), byrow=TRUE, ncol = length(DEMENAMES) )
		colnames(I) <- DEMENAMES 
		I_t <- rowSums( I[, which(colnames(I) != "Source")] )
		I_msm_t <- rowSums( I[, MSMDEMES] )
		I <- I[nrow(I),which(colnames(I) != "Source")]
		p <- sum(I) / ( sim[[5]][length(sim[[5]])] +  sum(I) )
		print(p) 
		t <- 1930 + (sim[[1]] - sim[[1]][1]) / 365
		par(mfcol=c(1,2))
		plot( t, I_t , xlim=c(1975,2014) ) 
		abline(v = 1975, col = 'red' )
		plot( t, I_msm_t , xlim=c(1975,2014) ) 
		abline(v = 1975, col = 'red' )
	}
	print('dim reduce sim time ' )
	print(system.time( {
		# likelihood goes from 6 sec to 1 sec
		tfgy <- dim.reduce.sim0( tfgy )
	}))
	clls <- foreach( bdt = iter(BDTS_dr0), .inorder=F) %do% 
	{
		cll <-  coalescent.log.likelihood.fgy(  bdt
		 , times=tfgy[[1]]
		 , births=tfgy[[2]]
		 , migrations=tfgy[[3]]
		 , demeSizes=tfgy[[4]]
		 , integrationMethod = 'adams'
		 , censorAtHeight=25 * 365
		 , forgiveAgtY=1 #.2
		 , returnTree=FALSE
		)
		print(date())
		print(cll)
		cll
	}
	cll <- logx2meanlogx( unlist( clls )  )
	print( c( date(), theta, cll) )
	list( cll, tfgy )
}



#  OBJFUNS
co.log.lik_sim.dr0.mcwm <- function( theta  ) 
{
	tfgy0 <- simulate(theta)
	bdt <- sample(BDTS_dr0, size = 1)[[1]]
	{ 
			coalescent.log.likelihood.fgy(  bdt
			 , times=tfgy[[1]]
			 , births=tfgy[[2]]
			 , migrations=tfgy[[3]]
			 , demeSizes=tfgy[[4]]
			 , integrationMethod = 'adams'
			 , censorAtHeight = 20 * 365
			 , forgiveAgtY = .25
			 , returnTree=FALSE
			) -> cll 
	}
	list( cll, tfgy0 )
}

lnlik.dr0.mcwm <- function(theta) { #
	ll_sim <- tryCatch({
	  co.log.lik_sim.dr0.mcwm( theta ) #}))
	 }, error = function(e) list(-Inf, NA) )
	
	  dp <- dprior( theta , sim = ll_sim[[2]], log = TRUE) 
	dp +  ll_sim[[1]]
}

co.log.lik_sim.dr0 <- function( theta  ) 
{ #NOTE only 10 trees
	tfgy0 <- simulate(theta) 
	tfgy <- dim.reduce.sim0( tfgy0 )
	#clls <- sapply(BDTS_dr0[1:10], function(bdt ) 
	clls <- sapply(BDTS_dr0[1:10], function(bdt ) #
	{
		cll <-  tryCatch({
			coalescent.log.likelihood.fgy(  bdt
			 , times=tfgy[[1]]
			 , births=tfgy[[2]]
			 , migrations=tfgy[[3]]
			 , demeSizes=tfgy[[4]]
			 , integrationMethod = 'adams'
			 , censorAtHeight=20 * 365
			 , forgiveAgtY=.25 #.2
			 , returnTree=FALSE
			)
		}, error = function(e) -Inf)
		cll
	})
	cll <- logx2meanlogx( unlist( clls )  )
	list( cll, tfgy0 )
}
lnlik.dr0 <- function(theta) { #
	ll_sim <- tryCatch({
	  co.log.lik_sim.dr0( theta )
	 }, error = function(e) list(-Inf, NA) )	
	dp <- dprior( theta , sim = ll_sim[[2]], log = TRUE) 
	dp +  ll_sim[[1]]
}

lnlik.nongen0 <- function ( theta )
{# just dprior.sim
	tryCatch({
		dprior( theta , sim = simulate(theta), log = TRUE) 
	}, error = function(e) -Inf)
}


MODEL5FITS <<- NULL

#~ lnlik <- lnlik.dr0.mcwm
lnlik <- lnlik.dr0
# /OBJFUNS


