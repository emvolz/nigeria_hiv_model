# see README 

require(rcolgem) 
source('rcolgem.R') #includes coratemodifier


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

U0_DEMES <- c('UM0', 'UF0', 'Umsm0')

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
DELTAT <- 31 #93 #31 #61 # days
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
SOURCESIZE <- 1e7
# default MSM parameters  
# 4.310345e-05 #1.436782e-05 # proportion of men MSM --- stef 1pc  --- considering just abuja 3e6*1e-2 / 2 = 15k --- 15k / 174e6/2 
#~  DHS: 25476 total msm in 2012, all nigeria; 
#~ > 25476 / 45.24e6 # dhs report / adult male nigeria 
#~ [1] 0.00056313
#~ > 25476 / 16.5e6 # dhs report / urban adult male nigeria
#~ [1] 0.001544
#~ PMSM <- .005 
MIN_PMSM <- .001544
MAX_PMSM <- .01 
IMMIGRATION_LB <- 1 / (30*365)
IMMIGRATION_UB <- 1 / (5 * 365 )
SOURCEGROWTHRATE_LB <- 1/(365*3.5)
SOURCEGROWTHRATE_UB <- 1/( 365*10/12 )
MAX_W <- 0.30 
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
BETA    <- .85 / 365 / 2 #1 * .58035640 / 365
BETAMSM <- 1 * BETA 
RATIOPOPSIZENEWOLD <-  .5 # ? calibrated to give realistic prevalence # 3e6 / 174e6
DEFAULTPARAMETERS <- c( y0 = .06973729
	 , beta = BETA
	 ,  w1 = 1/26.04, w2 = 7/ 26.04
	 , eps = .75 
	 , epsilonT = .5*(EPSTUB -EPSTLB )
	 , immigrationRate = 1 / 10 / 365
	 , assortprob = ASSORTPROB
	 , assortprob_msm = ASSORTPROB
	 , betamsm = BETAMSM
	 , gamma0 = 1/365 # EHI
	 # shrink chronc&aids periods by 9 months, or factor of 1-.75/12  
	 # 1/ ( (1-.75/12) /(0.157 / 365 )  )
	, gamma1 = 1/ ( (1-.75/12) /(0.157 / 365 )  )
	, gamma2 = 1/ ( (1-.75/12) /(0.350 / 365 )  )
	, gamma3 = 1/ ( (1-.75/12) /(0.282 / 365 )  )
	, gamma4 = 1/ ( (1-.75/12) /(.434 / 365 )  ) #aids 
	#~ 	, gamma1 = 0.157 / 365 
	#~ 	, gamma2 = 0.350 / 365 
	#~ 	, gamma3 = 0.282 / 365
	#~ 	, gamma4 = .434 / 365 #aids 
	, pstartstage1 = 0.58
	, pstartstage2 = 0.23  #reg4 combine these 
	, pstartstage3 = 0.16
	, pstartstage4 = 0.03
	, treBreak = TRE_BREAK
	, treRateIntervention = 0
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
	tt <- min(2015, tdays2calyears(t) )
	if (tt < 1980) return(ns_parmlist$beta) 
	with ( ns_parmlist , 
		#min( 1, 1 + ( 1 + exp(-fB*(1980-fM) ) )^(-1/fV)  - ( 1 + exp(-fB*(tt-fM) ) )^(-1/fV))
		max(0, beta + betaSlope * (tt - 1980))
	)
}

f.t.msm <- function(t, ns_parmlist )
{ 
	tt <- min(2015, tdays2calyears(t) )
	if (tt < 1980) return(ns_parmlist$betamsm) 
	with ( ns_parmlist , 
		max(0, betamsm + betamsmSlope * (tt - 1980))
	)
}

#~ p.t.msm <- function(t, ns_parmlist)
#~ {
#~ 	tt <- tdays2calyears( t )
#~ 	if (tt < 1980) return(1e-6)
#~ 	with (ns_parmlist, 
#~ 		max(1e-6, min(pmsmFin,  pmsmFin * (tt - 1980) / (2014 - 1980)) )
#~ 	)
#~ }

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
		
		foi <- f.t( t, parms )
		foi_msm <- f.t.msm(t, parms )
		
		w <- rep( c(1, rep(parms$w1,3), parms$w2 ), 3 )
		
		f_m2f <-   foi * w * y[MDEMES]  * (Sf  / Nf)  
		
		#note in following that assortprob contacts are reserved for within-type; remainder are allocated randomly
		f_f2m <-  parms$eps * foi * w * y[FDEMES] * 
		  (parms$assortprob * Sm  / (Nm) +  (1-parms$assortprob)*Sm/(Nm+Nmsm) )
		f_f2msm <-  parms$eps * foi * w * y[FDEMES] * (1-parms$assortprob) *  Smsm  / (Nmsm+Nm)
		
		f_msm2msm <- foi_msm * w * y[MSMDEMES] * 
		  (parms$assortprob_msm * Smsm / Nmsm + (1-parms$assortprob_msm) * Smsm / (Nmsm + Nf ) )
		f_msm2f   <- foi_msm * w * y[MSMDEMES] * (1-parms$assortprob_msm) * Sf / (Nmsm + Nf )
		
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
	#yy <- setNames( pmax(0, yy) , names(yy) )
#NOTE min pop size 
yy <- setNames( pmax(1e-1, yy) , names(yy) )
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
	pmsm <- parms$pmsmFin #p.t.msm( t, parms) 
	#with(as.list(y), 
	{
		Nf  <- y['Sf'] + If
		Nm  <- y['Sm'] + Im 
		Nmsm <- y['Smsm'] + Imsm
		N   <- Nf + Nm + Nmsm
		b   <- log( N.t(t+1) / N.t(t) ) +
		  ( parms$gamma4 * sum(y[AnT_DEMES]) + parms$gamma4 * parms$treBreak * sum(y[AT_DEMES]) ) /N + 
		  MU 
		d['Sm'] <- b * (1-pmsm) * N / 2 - sum(csFF[c(MDEMES)])   - MU*y['Sm']
		d['Smsm'] <- b * (pmsm) * N / 2 - sum(csFF[c(MSMDEMES)]) - MU*y['Smsm']
		d['Sf'] <- b * N / 2 - sum(csFF[c(FDEMES)]) - MU*y['Sf']
	}#)
	
	d[NONSOURCEDEMES] <- d[NONSOURCEDEMES] - MU * y[NONSOURCEDEMES]
	
	# aids deaths 
	d[AnT_DEMES] <- d[AnT_DEMES] - parms$gamma4 * y[AnT_DEMES] 
	d[AT_DEMES]  <- d[AT_DEMES]  - parms$gamma4 * parms$treBreak * y[AT_DEMES] 
#~ if (t > TIMES[197]) browser()
	list(d)
	
}


simulate <- function(theta)
{ # theta : vector 
	params <- DEFAULTPARAMETERS
	params[names(theta)] <- theta 
	params['pmsmFin'] <- MIN_PMSM + params['pmsmFin'] * (MAX_PMSM - MIN_PMSM) 
	params['w1'] <-  0 + params['w1'] * (MAX_W - 0) 
	params['w2'] <-  0 + params['w2'] * (MAX_W - 0) 
	params['immigrationRate'] <- IMMIGRATION_LB + (IMMIGRATION_UB - IMMIGRATION_LB) * params['immigrationRate']
	params['sourceGrowthRate'] <- SOURCEGROWTHRATE_LB + (SOURCEGROWTHRATE_UB - SOURCEGROWTHRATE_LB) * params['sourceGrowthRate'] 
	
	I0 <- params['y0'] #* (MAXY-MINY) + MINY
	S0 <- N.t(MAXSAMPLETIME - MAXHEIGHT)  #* exp(-rdemographic*(MAXHEIGHT))
	times <- seq(MAXSAMPLETIME - MAXHEIGHT, MAXSAMPLETIME+DELTAT, by=DELTAT)
	y0 <- setNames( rep(0, length(STATENAMES)), STATENAMES)
	y0['Sf'] <- S0 / 2
	y0['Smsm'] <- params['pmsmFin'] * S0/2 #1e-6#
	y0['Sm'] <- S0/2 # (1 - params['pmsm']) * S0/2
	y0[UFDEMES] = y0[UMDEMES] = y0[UMSMDEMES]  <- I0 #/ 1e2 #
	#y0['UM0'] = y0['UF0'] <- I0
	odeparms <- as.list(params)
	odeparms$t0 <- MAXSAMPLETIME - MAXHEIGHT
	
	# to speed things up, will use rk4 and solve ODEs only after expansiontime
	times2 <- times[ times > EXPANSIONTIME] 
	y0['Sm'] = y0['Sf'] <-  N.t(EXPANSIONTIME) / 2
	
	o0 <- cbind( times[ times <= EXPANSIONTIME],  matrix(y0,  byrow = TRUE, nrow = length(times) - length(times2), ncol = length(STATENAMES)) )
	o1 <- ode( y=y0, times=times2, func = dy, parms = odeparms, method='adams')
	o <- rbind(o0, o1)
	#class( o) <- 'deSolve'
#~ 	o <- ode( y=y0, times=times2, func = dy, parms = odeparms, method='adams')
#~ X11(); plot(o)
	colnames(o)[2:ncol(o)] <- STATENAMES
	oo <- o[,2:ncol(o)]
	Fs <- lapply( 1:nrow(o), function(k) .F.from.state(o[k,1], oo[k,], odeparms) ) #\times DELTAT  ? ?  No. 
	Gs <- lapply( 1:nrow(o), function(k) .G.from.state(o[k,1], oo[k,], odeparms) )
	Ys <- lapply( 1:nrow(o), function(k) .Y.from.state(o[k,1], oo[k,], odeparms) )
#~ browser()
	list(times, Fs, Gs, Ys, o[,'Sm']+o[,'Sf']+o[,'Smsm'], o)
}



LOG_SCALE_PARAMETERS <- c( 'y0'
  ,  'beta',  'fB', 'fM', 'fV'
  , 'betamsm',  'fmsmM' #'fmsmB'can be negative
  , 'dK', 'dB', 'dM'
  , 'treRate'
  , 'coRateModifier'
) #
LOGISTIC_SCALE_PARAMETERS <- c(  'assortprob', 'assortprob_msm', 'w1', 'w2', 'pmsm', 'immigrationRate', 'sourceGrowthRate', 'epsilonD', 'epsilonT', 'pmsmFin')
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
	  , w1 =  runif(1, 0, 1) #rbeta(1, 1, 2 ) 
	  , w2 = runif(1, 0, 1) #rbeta( 1, 1, 2 )
	  , y0 = runif(1, 0, 10 )
	  , immigrationRate = runif(1, 0, 1) # 1 / (365 * runif( 1, 5, 30) )
	  , sourceGrowthRate =  runif(1, 0, 1) #= 1 / (365 * runif(1, 10/12, 3.5) ) #1/x prior based on empirical observations of doubling time
	  
	  	, beta = rlnorm(1,  log( BETA ) , 1/2) 
	  	, betaSlope = rnorm(1, 0, BETA / (35) /2)
		
		, betamsm = rlnorm(1,  log( BETAMSM ) , abs(log(2* .58))) ## good
		, betamsmSlope = rnorm( 1, 0 , BETAMSM/ (35) /2)
		
		, dK =  rlnorm( 1, log(1/10/365), .6 ) #median 4yrs
		, dB = rlnorm( 1, log(.5), .4  )  
		, dM = rnorm(1, 2006, sd= 1) 
		
		, treRate = rlnorm( 1, log(1/(2*365)), .4 ) 
		
		, epsilonD = runif(1)
		, coRateModifier = rlnorm(1, log(1), .1 )
		, pmsmFin = runif(1)
	)
}

rprior.1 = rprior.vec.1 <- function()
{ #
	c( assortprob =  rbeta(1,  5,1 )
	  , assortprob_msm = rbeta(1,  5,1 )
	  , w1 =  runif(1, 0, 1) #rbeta(1, 1, 2 ) 
	  , w2 = runif(1, 0, 1) #rbeta( 1, 1, 2 )
	  , y0 = runif(1, 0, 10 )
	  , immigrationRate = runif(1, 0, 1) # 1 / (365 * runif( 1, 5, 30) )
	  , sourceGrowthRate =  runif(1, 0, 1) #= 1 / (365 * runif(1, 10/12, 3.5) ) #1/x prior based on empirical observations of doubling time
	  
	  	, beta = rlnorm(1,  log( BETA ) , 1/2) 
	  	, betaSlope = rnorm(1, 0, .5 / (35) /2)
		
		, betamsm = rlnorm(1,  log( BETAMSM ) , abs(log(2* .58))) ## good
		, betamsmSlope = rnorm( 1, 0 , .5 / (35) /2)
		
		, dK =  rlnorm( 1, log(1/10/365), .6 ) #median 4yrs
		, dB = rlnorm( 1, log(.5), .4  )  
		, dM = rnorm(1, 2006, sd= 1) 
		
		, treRate = rlnorm( 1, log(1/(2*365)), .4 ) 
		
		, epsilonD = runif(1)
		, pmsmFin = runif(1)
	)
}

dprior.sim <- function(sim, log =TRUE){
	d <- ifelse( log, 0, 1)
		i2001 <- which.min( abs(TIMES_CALYEARS - 2001)  )
		i2008 <- which.min( abs(TIMES_CALYEARS - 2008)  ) 
		i2012 <- which.min( abs(TIMES_CALYEARS - 2012)  )
		
		I2001 <- sum( sim[[4]][[i2001]][NONSOURCEDEMES] ) 
		I2012 <- sum( sim[[4]][[i2012]][NONSOURCEDEMES] ) 
		
		p2001 <- (I2001) / ( sim[[5]][ i2001 ] +  (I2001) )
		p2012 <- (I2012) / ( sim[[5]][ i2012 ] +  (I2012) )
		#d <- ifelse(log,  d + dnorm( p, mean=.017, sd=.0015, log = log ), d *  dnorm( p, mean=.017, sd=.0015, log = log ) )
		# unaids dhs report 
		d <- ifelse(log
		  ,  d + dnorm( p2001, mean=.035, sd=.01 /(2*1.96), log = log ) + dnorm( p2012, mean=.031, sd=.007 /(2*1.96), log = log )
		  ,  d * dnorm( p2001, mean=.035, sd=.01 /(2*1.96), log = log ) * dnorm( p2012, mean=.031, sd=.007 /(2*1.96), log = log )
		)
		

		# sim prior on treatment 
		# proportion eligible on treatment
		# this involves model-estimated values by unaids; more reliable to use num treated 
		#~ 		TELIG_DEMES <- c(DMDEMES[4:5], DFDEMES[4:5], DMSMDEMES[4:5] ) 
		#~ 		Nelig <- sum( sim[[4]][[i2012]][TELIG_DEMES] ) 
		#~ 		Ntreat <- sum( sim[[4]][[i2012]][T_DEMES] ) 
		#~ 		ptreat <- Ntreat / (Ntreat + Nelig )
		#~ 		d <- ifelse(log
		#~ 		  ,  d + dnorm( ptreat, mean=.36, sd=(.39-.33) /(2*1.96), log = log )
		#~ 		  ,  d *dnorm( ptreat, mean=.36, sd=(.39-.33)/(2*1.96), log = log ) 
		#~ 		)
		
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
		
		#BRIDGING TRUST 
		# prevalence msm 
		i2010_end <- which.min( abs(TIMES_CALYEARS - 2011)  ) 
		i2014_end <- which.min( abs(TIMES_CALYEARS - 2015)  ) 
        #        geom_errorbar( aes(x = 2010, ymin=25.5/100,ymax=45.9/100) ) + geom_point( aes(x = 2010, y = 34.9/100) )  + 
        #        geom_errorbar( aes(x = 2014, ymin=39.9/100,ymax=49.1/100, guide=F), colour='green' ) + geom_point( aes(x = 2014, y = 44.5/100, guide=F), colour='green' ) 
		Imsm_2010_end <- sum( sim[[4]][[i2010_end]][MSMDEMES] ) 
		Imsm_2014_end <- sum( sim[[4]][[i2014_end]][MSMDEMES] ) 
		
		prevmsm2010_end <- (Imsm_2010_end) / ( sim[[5]][ i2010_end ] +  (Imsm_2010_end) )
		prevmsm2014_end <- (Imsm_2014_end) / ( sim[[5]][ i2014_end ] +  (Imsm_2014_end) )
		
		prevmsm_sd2010_end <- (45.9/100 - 25.5/100) / 1.96
		prevmsm_center2010_end <- 34.9/100
		prevmsm_sd2014_end <- (49.1/100 - 39.9/100) / 1.96
		prevmsm_center2014_end <- 44.5/100
		
		if ( prevmsm2010_end > 25.5/100 & prevmsm2010_end < 45.9/100){
			d <-  ifelse(log
			  ,  d + dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log )
			  ,  d * dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log ) 
			)
#~ 			print('dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log )')
#~ 			print( dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log ) )
		} else{
			# prevalence penalty, smaller sd
			d <-  ifelse(log
			  ,  d + dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end/1e3, log = log )
			  ,  d * dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end/1e3, log = log ) 
			)
#~ 			print('dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end/1e3, log = log )')
#~ 			print( dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end/1e3, log = log ) )
		}
		
		if (prevmsm2014_end > 39.9/100 & prevmsm2014_end < 49.1/100)
		{
			d <-  ifelse(log
			  ,  d + dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log )
			  ,  d * dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log ) 
			)
#~ 			print('dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log )')
#~ 			print(dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log ))
		} else{
			# prevalence penalty, smaller sd
			d <-  ifelse(log
			  ,  d + dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end/1e3, log = log )
			  ,  d * dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end/1e3, log = log ) 
			)
#~ 			print('dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end/1e3, log = log )')
#~ 			print( dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end/1e3, log = log ) )
		}
		
		#incidence msm 
		#geom_errorbar( aes(x = 2014, ymin=7/100,ymax=27.8/100, guide=F), colour='green' ) + geom_point( aes(x = 2014, y = 13.9/100, guide=F), colour='green' ) 
		i2014_mid <- which.min( abs(TIMES_CALYEARS - 2014.5)  ) 
		smsm_2014 <- sim[[6]][i2014_mid, 'Smsm'] 
		pcinc_msm_2014 <- sum(sim[[2]][[i2014_mid]][NONSOURCEDEMES, MSMDEMES]) * 365 / smsm_2014
		pcincmsm_sd <- (27.8/100 - 7/100) / 1.96
		pcincmsm_center <- 13.9/100
		if (pcinc_msm_2014 > 7/100 & pcinc_msm_2014 < 27.8/100) {
			d <-  ifelse(log
			  ,  d + dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log )
			  ,  d * dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log ) 
			)
#~ 			print('dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log )')
#~ 			print(dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log ))
		} else{# penalty, smaller sd
			d <-  ifelse(log
			  ,  d + dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd/1e3, log = log )
			  ,  d * dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd/1e3, log = log ) 
			)
#~ 			print( 'dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd/1e3, log = log ') 
#~ 			print( dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd/1e3, log = log ) )
		}
	d
}

dprior <- function( theta, sim = NA, log = TRUE)
{
	d <- with(as.list(theta), {
		ps <- c( dbeta(assortprob, 5,1, log = log)
		  , dbeta(assortprob_msm, 5,1, log = log)
		  , dunif(w1, 0, 1, log = log) #dbeta( w1, 1, 2, log = log)
		  , dunif(w2, 0, 1, log = log )#dbeta( w2, 1, 2 , log = log) 
		  , dunif( y0, 0, 10, log = log)
		  , dunif(immigrationRate, 0, 1 , log = log ) #dunif(  1 / immigrationRate, 5*365, 30*365, log = log )
		  , dunif(sourceGrowthRate, 0, 1, log = log )# dunif(  1 / sourceGrowthRate, 365*10/12, 365*3.5, log = log )
		  
			, beta = dlnorm(beta,  log( BETA ) , 1/2, log = log ) 
			, betaSlope = dnorm( betaSlope, 0, .5 / (35) /2, log = log)
			
			, betamsm = dlnorm(betamsm,  log( BETAMSM ) , abs(log(2* .58)), log = log) 
			, betamsmSlope = dnorm( betamsmSlope, 0 , .5 / (35) /2, log = log)
			
			, dK = dlnorm( dK, log(1/10/365), .6, log = log  ) # dunif( dK, 0, 1/8, log = log)
			, dB = dlnorm( dB, log(.5), .4 , log = log ) 
			, dM = dnorm(dM, 2006, sd= 1, log = log)
			
			, treRate = dlnorm( treRate, log(1/(2*365)), .4 , log = log ) 
			
		  , epsilonD = dunif(epsilonD, 0, 1, log = log)
		  , coRateModifier = dlnorm(coRateModifier, log(1), .1, log = log )
		  , pmsmFin = dunif(pmsmFin, 0, 1, log = log)
		)
		ifelse( log, sum(ps), prod(ps) )
	})
	if ( length(sim) > 1 )
	{
		d <- ifelse(log, d + dprior.sim( sim, log ), d * dprior.sim(sim, log))
	}
	d
}

dprior.1 <- function( theta, sim = NA, log = TRUE)
{
d <- with(as.list(theta), {
		ps <- c( dbeta(assortprob, 5,1, log = log)
		  , dbeta(assortprob_msm, 5,1, log = log)
		  , dunif(w1, 0, 1, log = log) #dbeta( w1, 1, 2, log = log)
		  , dunif(w2, 0, 1, log = log )#dbeta( w2, 1, 2 , log = log) 
		  , dunif( y0, 0, 10, log = log)
		  , dunif(immigrationRate, 0, 1 , log = log ) #dunif(  1 / immigrationRate, 5*365, 30*365, log = log )
		  , dunif(sourceGrowthRate, 0, 1, log = log )# dunif(  1 / sourceGrowthRate, 365*10/12, 365*3.5, log = log )
		  
			, beta = dlnorm(beta,  log( BETA ) , 1/2, log = log ) 
			, betaSlope = dnorm( betaSlope, 0, .5 / (35) /2, log = log)
			
			, betamsm = dlnorm(betamsm,  log( BETAMSM ) , abs(log(2* .58)), log = log) 
			, betamsmSlope = dnorm( betamsmSlope, 0 , .5 / (35) /2, log = log)
			
			, dK = dlnorm( dK, log(1/10/365), .6, log = log  ) # dunif( dK, 0, 1/8, log = log)
			, dB = dlnorm( dB, log(.5), .4 , log = log ) 
			, dM = dnorm(dM, 2006, sd= 1, log = log)
			
			, treRate = dlnorm( treRate, log(1/(2*365)), .4 , log = log ) 
			
		  , epsilonD = dunif(epsilonD, 0, 1, log = log)
		  , pmsmFin = dunif(pmsmFin, 0, 1, log = log)
		)
		ifelse( log, sum(ps), prod(ps) )
	})
	if ( length(sim) > 1 )
	{
		d <- ifelse(log, d + dprior.sim( sim, log ), d * dprior.sim(sim, log))
	}
	d
}



logx2meanlogx <- function(x)
{
	log( mean( exp( x-max(x) ) ) ) + max(x)
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
MSMDEMES0 <- DEMES0[ grepl( 'msm', DEMES0 ) ]
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


#  OBJFUNS
co.log.lik_sim.dr0.mcwm <- function( theta  ) 
{
	tfgy0 <- simulate.cpp(theta)
	tfgy <- dim.reduce.sim0( tfgy0 ) 
	bdt <- sample(BDTS_dr0, size = 1)[[1]]
	{ 
			coalescent.log.likelihood.fgy(  bdt
			 , times=tfgy[[1]]
			 , births=tfgy[[2]]
			 , migrations=tfgy[[3]]
			 , demeSizes=tfgy[[4]]
			 , integrationMethod = 'adams'
			 , censorAtHeight = 25 * 365
			 , forgiveAgtY = .25
			 , returnTree=FALSE
			 , coRateModifier = theta['coRateModifier']
			) -> cll 
	}
	#print( c( date(), theta, cll) )
	list( cll, tfgy0 )
}

lnlik.dr0.mcwm <- function(theta) { #
	ll_sim <- tryCatch({
	  co.log.lik_sim.dr0.mcwm( theta ) 
	 }, error = function(e) list(-Inf, NA) )
	
	  dp <- dprior( theta , sim = ll_sim[[2]], log = TRUE) 
	#print( c('dp', dp) )
	dp +  ll_sim[[1]]
}

# no coratemodifier,  uses func.1 priors and shorter censoring
co.log.lik_sim.dr0.mcwm.1 <- function( theta  ) 
{
	tfgy0 <- simulate.cpp(theta)
	tfgy <- dim.reduce.sim0( tfgy0 ) 
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
	#print( c( date(), theta, cll) )
	list( cll, tfgy0 )
}

lnlik.dr0.mcwm.1 <- function(theta) { #
	ll_sim <- tryCatch({
	  co.log.lik_sim.dr0.mcwm.1( theta ) 
	 }, error = function(e) list(-Inf, NA) )
	
	  dp <- dprior.1( theta , sim = ll_sim[[2]], log = TRUE) 
	#print( c('dp', dp) )
	dp +  ll_sim[[1]]
}

co.log.lik_sim.dr0 <- function( theta  ) 
{ #NOTE only x trees
	tfgy0 <- simulate.cpp(theta) 
	#tfgy <- tfgy0
	tfgy <- dim.reduce.sim0( tfgy0 )
	clls <- sapply(BDTS_dr0[1:10], function(bdt ) #
	#clls <- sapply(BDTS[1:2], function(bdt ) #
	{
		#cll <-  tryCatch({
		cll <-	coalescent.log.likelihood.fgy(  bdt
			 , times=tfgy[[1]]
			 , births=tfgy[[2]]
			 , migrations=tfgy[[3]]
			 , demeSizes=tfgy[[4]]
			 , integrationMethod = 'adams'
			 , censorAtHeight=25 * 365
			 , forgiveAgtY=0 #.2
			 , returnTree=FALSE
			 , coRateModifier = theta['coRateModifier']
			)
		#}, error = function(e) -Inf)
		ifelse( is.na(cll), -Inf, cll )
	})
	#print ('#################')
	#print(clls) 
	cll <- logx2meanlogx( unlist( clls )  )
	cll <- ifelse( is.na(cll), -Inf, cll )
	list( cll, tfgy0 )
}


lnlik.dr0 <- function(theta) { #
	ll_sim <- tryCatch({
	  co.log.lik_sim.dr0( theta )
	 }, error = function(e) list(-Inf, NA) )	
	dp <- dprior( theta , sim = ll_sim[[2]], log = TRUE) 
	dp +  ll_sim[[1]]
}




lnlik <- lnlik.dr0.mcwm.1
#~ lnlik <- lnlik.dr0
# /OBJFUNS



simulate.cpp <- function(theta)
{ # theta : vector 
	params <- DEFAULTPARAMETERS
	params[names(theta)] <- theta 
	params['pmsmFin'] <- MIN_PMSM + params['pmsmFin'] * (MAX_PMSM - MIN_PMSM) 
	params['w1'] <-  0 + params['w1'] * (MAX_W - 0) 
	params['w2'] <-  0 + params['w2'] * (MAX_W - 0) 
	params['immigrationRate'] <- IMMIGRATION_LB + (IMMIGRATION_UB - IMMIGRATION_LB) * params['immigrationRate']
	params['sourceGrowthRate'] <- SOURCEGROWTHRATE_LB + (SOURCEGROWTHRATE_UB - SOURCEGROWTHRATE_LB) * params['sourceGrowthRate'] 
	
	I0 <- params['y0'] 
	S0 <- N.t(MAXSAMPLETIME - MAXHEIGHT)  
	times <- seq(MAXSAMPLETIME - MAXHEIGHT, MAXSAMPLETIME+DELTAT, by=DELTAT)
	y0 <- setNames( rep(0, length(STATENAMES)), STATENAMES)
	y0['Sf'] <- S0 / 2
	y0['Smsm'] <- params['pmsmFin'] * S0/2 #1e-6#
	y0['Sm'] <- S0/2
	y0[UFDEMES] = y0[UMDEMES] = y0[UMSMDEMES]  <- I0 #/ 1e2 #
	
	odeparms <- as.list(params)
	odeparms$t0 <- MAXSAMPLETIME - MAXHEIGHT
	{ #add extra components to odeparms - needed for cpp code
		odeparms$m = m
		odeparms$DEMES <- DEMENAMES
		odeparms$imdemes <-  which ( DEMENAMES %in% MDEMES ) - 1
		odeparms$ifdemes <-  which ( DEMENAMES %in% FDEMES ) - 1
		odeparms$imsmdemes <-  which ( DEMENAMES %in% MSMDEMES ) - 1
		odeparms$itdemes <-  which ( DEMENAMES %in% T_DEMES ) -1 
		odeparms$iddemes <-  which ( DEMENAMES %in% D_DEMES ) - 1
		odeparms$w <- rep( c(1, rep(odeparms$w1,3), odeparms$w2 ), 3 ) 
		
		odeparms$i_um = which(DEMENAMES %in% UMDEMES) - 1
		odeparms$i_dm = which(DEMENAMES %in% DMDEMES) - 1
		odeparms$i_tm = which(DEMENAMES %in% TMDEMES) - 1
		
		odeparms$i_umsm = which(DEMENAMES %in% UMSMDEMES) - 1
		odeparms$i_dmsm = which(DEMENAMES %in% DMSMDEMES) - 1
		odeparms$i_tmsm = which(DEMENAMES %in% TMSMDEMES) - 1
		
		odeparms$i_uf = which(DEMENAMES %in% UFDEMES) - 1
		odeparms$i_df = which(DEMENAMES %in% DFDEMES) - 1
		odeparms$i_tf = which(DEMENAMES %in% TFDEMES) - 1
		
		odeparms$i_AnT = which(DEMENAMES %in% AnT_DEMES) - 1
		odeparms$i_AT = which(DEMENAMES %in% AT_DEMES) - 1
	}
	
	times2 <- times[ times > EXPANSIONTIME] 
	y0['Sm'] = y0['Sf'] <-  N.t(EXPANSIONTIME) / 2
	
	o0 <- cbind( times[ times <= EXPANSIONTIME],  matrix(y0,  byrow = TRUE, nrow = length(times) - length(times2), ncol = length(STATENAMES)) )
	o1 <- ode( y=y0, times=times2, func = dydt, parms = odeparms, method='adams')
	o <- rbind(o0, o1)
	
	colnames(o)[2:ncol(o)] <- STATENAMES
	oo <- o[,2:ncol(o)]
	Fs <- lapply( 1:nrow(o), function(k) F_from_state(o[k,1], oo[k,], odeparms) ) #\times DELTAT  ? ?  No. 
	Gs <- lapply( 1:nrow(o), function(k) G_from_state(o[k,1], oo[k,], odeparms) )
	Ys <- lapply( 1:nrow(o), function(k) .Y.from.state(o[k,1], oo[k,], odeparms) )
	list(times, Fs, Gs, Ys, o[,'Sm']+o[,'Sf']+o[,'Smsm'], o)
}



