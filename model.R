#~ author : Erik Volz <evolz@imperial.ac.uk>
#~ date : Nov 13, 2015

#~ This file provides functions for 
#~ 1) simulating HIV epidemic model incl. MSM & het males and females
#~ 2) computing priors, reducing dimension of model outputs
#~ 3) computing prob of nongenetic surveillance data (dprior.sim) 
#~ 4) computing prob of tree data using rcolgem package

#~ Requirements:
#~ 1) rcolgem package for phylodynamics: http://colgem.r-forge.r-project.org/
#~ 2) gfmodel11, compiled package for quickly simulating HIV model
#~ Other dependencies: Rcpp, RcppArmadillo, ape, deSolve

#~ NOTE 
#~ running this file unmodified will require data that is not included 
#~ 1)  'mergeTrees0.nwk' : sample of BEAST trees
#~ 2)  'sampleStates0.csv' : initial state of sampled patients
#~ if interested in using/adapting this model, or accessing full data set, contact: 
#~ Erik Volz <evolz@imperial.ac.uk>
#~ Man Charurat <MCharurat@ihv.umaryland.edu >

require(rcolgem) 
#~ source('rcolgem.R') #modified version w/ includes extra parameter; not needed for all analyses
require(gfmodel11)


# MODEL ORGANISATION AND DEMES
UMDEMES <- paste( 'UM', 0:4 , sep = '') #undiagnosed male
DMDEMES <- paste(sep='', 'DM', 0:4) #diagnosed male
TMDEMES <- paste(sep='', 'TM', 0:4) #treated male

UFDEMES <- paste( 'UF', 0:4 , sep = '') # undiag female
DFDEMES <- paste(sep='', 'DF', 0:4)  # diag female
TFDEMES <- paste(sep='', 'TF', 0:4) #treat female

UMSMDEMES <- paste( 'Umsm', 0:4 , sep = '')
DMSMDEMES <- paste(sep='', 'Dmsm', 0:4) 
TMSMDEMES <- paste(sep='', 'Tmsm', 0:4) 

MDEMES <- c( UMDEMES, DMDEMES, TMDEMES ) #male demes
FDEMES <- c( UFDEMES, DFDEMES, TFDEMES ) # fem demes
MSMDEMES <- c( UMSMDEMES, DMSMDEMES, TMSMDEMES ) #msm demes

E_DEMES <- c(UMDEMES[1], DMDEMES[1], TMDEMES[1], UFDEMES[1], DFDEMES[1], TFDEMES[1], UMSMDEMES[1], DMSMDEMES[1], TMSMDEMES[1] ) #early infection demes
C_DEMES <- c(UMDEMES[2:4], DMDEMES[2:4], TMDEMES[2:4], UFDEMES[2:4], DFDEMES[2:4], TFDEMES[2:4], UMSMDEMES[2:4], DMSMDEMES[2:4], TMSMDEMES[2:4] ) # chronic demes
A_DEMES <- c(UMDEMES[5], DMDEMES[5], TMDEMES[5], UFDEMES[5], DFDEMES[5], TFDEMES[5], UMSMDEMES[5], DMSMDEMES[5], TMSMDEMES[5] ) #aids demes
AnT_DEMES <- c(UMDEMES[5], DMDEMES[5], UFDEMES[5], DFDEMES[5], UMSMDEMES[5], DMSMDEMES[5]) # aids + nontreated
AT_DEMES <- c(TMDEMES[5], TFDEMES[5],  TMSMDEMES[5] ) # aids + treated

U_DEMES <- c( UMDEMES, UFDEMES, UMSMDEMES ) #undiag
D_DEMES <- c( DMDEMES, DFDEMES, DMSMDEMES ) #diag
T_DEMES <- c( TMDEMES, TFDEMES, TMSMDEMES ) #treated

U0_DEMES <- c('UM0', 'UF0', 'Umsm0')

DEMENAMES <- c( UMDEMES, DMDEMES, TMDEMES
  , UFDEMES, DFDEMES, TFDEMES
  , UMSMDEMES, DMSMDEMES, TMSMDEMES
  ,   'Source') #28 demes 
NONSOURCEDEMES <- DEMENAMES[-length(DEMENAMES) ]
STATENAMES <- c( DEMENAMES[-length(DEMENAMES)],  'Sm', 'Sf', 'Smsm') #incl susceptible categories
m <- length(DEMENAMES) #dimension of model
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

.Y.from.state <- function(t, y, parms)
{
	yy <- c(y, Source=source.size(t, parms )) 
	#yy <- setNames( pmax(0, yy) , names(yy) )
#NOTE min pop size 
yy <- setNames( pmax(1e-1, yy) , names(yy) )
	yy[DEMENAMES]
}


LOG_SCALE_PARAMETERS <- c( 'y0'
  ,  'beta',  'fB', 'fM', 'fV'
  , 'betamsm', 'fMmsm' , 'fVmsm' #'fBmsm'can be negative
  , 'dK', 'dB', 'dM'
  , 'treRate'
  , 'coRateModifier'
) #
LOGISTIC_SCALE_PARAMETERS <- c( 'assortprob', 'assortprob_msm', 'w1', 'w2', 'pmsm', 'immigrationRate', 'sourceGrowthRate', 'epsilonD', 'epsilonT', 'pmsmFin')
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




##################################
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
	  	, fB = rlnorm( 1, log(.25), 1 )
		, fM = rnorm(1,  mean = 1990, sd = 5)
		#, fV = rlnorm(1, log(1), .5)
		
		, betamsm = rlnorm(1,  log( BETAMSM ) , abs(log(2* .58))) ## good
		, fBmsm = rnorm( 1, (.25), .25 )
		, fMmsm = rnorm(1,  mean = 1990, sd = 5)
		#, fV = rlnorm(1, log(1), .5)
		
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
	  	, fB = rlnorm( 1, log(.25), 1 )
		, fM = rnorm(1,  mean = 1990, sd = 5)
		#, fV = rlnorm(1, log(1), .5)
		
		, betamsm = rlnorm(1,  log( BETAMSM ) , abs(log(2* .58))) ## good
		, fBmsm = rnorm( 1, (.25), .25 )
		, fMmsm = rnorm(1,  mean = 1990, sd = 5)
		#, fV = rlnorm(1, log(1), .5)
		
		, dK =  rlnorm( 1, log(1/10/365), .6 ) #median 4yrs
		, dB = rlnorm( 1, log(.5), .4  )  
		, dM = rnorm(1, 2006, sd= 1) 
		
		, treRate = rlnorm( 1, log(1/(2*365)), .4 ) 
		
		, epsilonD = runif(1)
		, pmsmFin = runif(1)
	)
}

# based on nongenetic data
dprior.sim <- function(sim, log =TRUE, debug = FALSE){
	d <- ifelse( log, 0, 1)
	i2001 <- which.min( abs(TIMES_CALYEARS - 2001)  )
	i2008 <- which.min( abs(TIMES_CALYEARS - 2008)  ) 
	i2012 <- which.min( abs(TIMES_CALYEARS - 2012)  )
	i2014m <- which.min( abs(TIMES_CALYEARS - 2013.5)  ) #midpoint of year up to '14
	
	I2001 <- sum( sim[[4]][[i2001]][NONSOURCEDEMES] ) 
	I2012 <- sum( sim[[4]][[i2012]][NONSOURCEDEMES] ) 
	
	p2001 <- (I2001) / ( sim[[5]][ i2001 ] +  (I2001) )
	p2012 <- (I2012) / ( sim[[5]][ i2012 ] +  (I2012) )
	
	
	#d <- ifelse(log,  d + dnorm( p, mean=.017, sd=.0015, log = log ), d *  dnorm( p, mean=.017, sd=.0015, log = log ) )
	# unaids dhs report 
	p2001sd <- ifelse( (p2001 > .035 - 1.96 * .01 & p2001 < .035 + 1.96*.01), .01 / (2*1.96), .01 / (2*1.96) / 1e3)
	{
		d <- ifelse(log
		 , d + dnorm( p2001, mean = .035, sd = p2001sd, log = log)
		 , d * dnorm( p2001, mean = .035, sd = p2001sd, log = log)
		)
	}
	
	p2012sd <- ifelse( (p2012 > .031 - 1.96 * .007 & p2012 < .031 + 1.96*.007), .007 / (2*1.96), (.007 / (2*1.96)) / 1e3)
	#{
		d <- ifelse(log
		 , d + dnorm( p2012, mean = .031, sd = p2012sd, log = log)
		 , d * dnorm( p2012, mean = .031, sd = p2012sd, log = log)
		)
	#} 
	
	# see written notes 13/8 -- docking to unaids incidence in '14
	if (F)
	{ #removing- gave trajectories with unrealisticly late time of peak prev in GP
		sgp_2014 <- sim[[6]][i2014m, 'Sm'] + sim[[6]][i2014m, 'Sf']
		inc2014 <- sum(sim[[2]][[i2014m]][NONSOURCEDEMES, c(MDEMES,FDEMES)]) * 365 / sgp_2014
		
		i2014m_sd <- ifelse(inc2014 > .001163
		  , (.00152-.001163)/(2*1.96)
		  , ((.00152-.001163)/(2*1.96))/1e3
		) 
		d <- ifelse(log
		 , d + dnorm( inc2014, mean = .00134, sd = i2014m_sd, log = log)
		 , d * dnorm( inc2014, mean = .00134, sd = i2014m_sd, log = log)
		)
	}
	
	if (debug ){
		print ('p2001')
		print( c(  p2001, .035, dnorm( p2001, mean = .035, sd = p2001sd, log = log) ) )
		print( 'p2012')
		print( c( p2012, .031,dnorm( p2012, mean = .031, sd = p2012sd, log = log) ) )
	}
	
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
	
	prevmsm2010_end <- (Imsm_2010_end) / ( sim[[6]][i2010_end, 'Smsm'] +  (Imsm_2010_end) )
	prevmsm2014_end <- (Imsm_2014_end) / ( sim[[6]][i2014_end, 'Smsm'] +  (Imsm_2014_end) )
	
	prevmsm_sd2010_end <- (45.9/100 - 25.5/100) / 1.96
	prevmsm_center2010_end <- 34.9/100
	prevmsm_sd2014_end <- (49.1/100 - 39.9/100) / 1.96
	prevmsm_center2014_end <- 44.5/100
	
	if (!( prevmsm2010_end > 25.5/100 & prevmsm2010_end < 45.9/100)){
		prevmsm_sd2010_end <- prevmsm_sd2010_end/ 1e3
	}
	d <-  ifelse(log
	  ,  d + dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log )
	  ,  d * dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log ) 
	)
	if (debug){
		print('prevmsm2010_end')
		print ( prevmsm2010_end)
		print( dnorm( prevmsm2010_end, mean=prevmsm_center2010_end, sd=prevmsm_sd2010_end, log = log ) )
	}
	
	if (!(prevmsm2014_end > 39.9/100 & prevmsm2014_end < 49.1/100))
	{
		prevmsm_sd2014_end <- prevmsm_sd2014_end / 1e3
	}
	d <-  ifelse(log
	  ,  d + dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log )
	  ,  d * dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log ) 
	)
	if (debug){
		print('prevmsm2014_end' )
		print( prevmsm2014_end)
		print(dnorm( prevmsm2014_end, mean=prevmsm_center2014_end, sd=prevmsm_sd2014_end, log = log ))
	}
	
	#incidence msm 
	#geom_errorbar( aes(x = 2014, ymin=7/100,ymax=27.8/100, guide=F), colour='green' ) + geom_point( aes(x = 2014, y = 13.9/100, guide=F), colour='green' ) 
	i2014_mid <- which.min( abs(TIMES_CALYEARS - 2014.5)  ) 
	smsm_2014 <- sim[[6]][i2014_mid, 'Smsm'] 
	pcinc_msm_2014 <- sum(sim[[2]][[i2014_mid]][NONSOURCEDEMES, MSMDEMES]) * 365 / smsm_2014
	pcincmsm_sd <- (27.8/100 - 7/100) / 1.96
	pcincmsm_center <- 13.9/100
	if (! (pcinc_msm_2014 > 7/100 & pcinc_msm_2014 < 27.8/100)) {
		pcincmsm_sd <- pcincmsm_sd / 1e3
	}
	d <-  ifelse(log
	  ,  d + dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log )
	  ,  d * dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log ) 
	)
	if (debug)
	{
		print( 'pcinc_msm_2014' )
		print(dnorm( pcinc_msm_2014, mean=pcincmsm_center, sd=pcincmsm_sd, log = log ))
	}
	#browser()
	unname( d )
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
			, fB = dlnorm( fB, log(.25), 1, log = log  ) #dunif(fB, 0, 1, log = log)
			, fM = dnorm(fM,  mean = 1990, sd = 5, log = log)
			
			, betamsm = dlnorm(betamsm,  log( BETAMSM ) , abs(log(2* .58)), log = log) 
			, fBmsm = dnorm( fBmsm, (.25), .25 , log = log )
			, fMmsm = dnorm(fM,  mean = 1990, sd = 5, log = log)
			
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
			, fB = dlnorm( fB, log(.25), 1, log = log  ) #dunif(fB, 0, 1, log = log)
			, fM = dnorm(fM,  mean = 1990, sd = 5, log = log)
			
			, betamsm = dlnorm(betamsm,  log( BETAMSM ) , abs(log(2* .58)), log = log) 
			, fBmsm = dnorm( fBmsm, (.25), .25 , log = log )
			, fMmsm = dnorm(fM,  mean = 1990, sd = 5, log = log)
			
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

#######################################################
# DIMENSION REDUCTION
# condensed version of model outputs used for phylodynamic inference
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






########################################################################
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
	#~ ll_sim <-  co.log.lik_sim.dr0.mcwm.1( theta ) 
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
	y0['Sf'] <- (1-params['pmsmFin']) * S0 / 2
	y0['Smsm'] <- params['pmsmFin'] * S0/2 
	y0['Sm'] <- S0/2
	y0[UFDEMES] = y0[UMDEMES]  <- I0 
	y0[UMSMDEMES] <- I0 * params['pmsmFin']
	
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


########################################################################
# example usage / debugging
if (F) {
	theta <- rprior.1()
	require(Rcpp)
	sourceCpp('model.cpp')
	s <- simulate.cpp( theta )
	print( dprior.sim( s ) )
}
