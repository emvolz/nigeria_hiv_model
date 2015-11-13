// [[Rcpp]]
#include <Rcpp.h>

using namespace Rcpp; 

static const double  SOURCESIZE = 1.e7;
static const double  MAXSAMPLETIME = 12395.;
static const double EXPANSIONTIME = -2141.;
const double N1 = 1677199.;
const double RDEMOGRAPHIC = 6.925427e-05;
const double MU = 6.08828e-05;
static const double T2D_RATE = 0.0006807161;

// intervention stuff
// NOTE intervention in this model is for both GP & MSM
// different model will handle intervention in MSM only
static const double TINTERVENTIONSTART_YRS = 2016.;
static const double TINTERVENTIONFIN_YRS  = 2021.; // scale up intervention over five year period (linear)
static const double DRFIN = 1. / 365.;
static const double TRFIN = 1. / 365.;
static const double RETENTION_PROB_FIN = .90; // # retention on art within 1 year; .78 = exp(-r*365 ) ; 
static const double T2D_RATE_FIN =  -log(RETENTION_PROB_FIN) / 365. ;

double N_t(double t){
	return N1 * exp(-RDEMOGRAPHIC * (MAXSAMPLETIME - t));
}

double tdays2calyears(double tdays){
	// note offset corresponds to birthday - 1980
	return 1980 + 0.8630137 + tdays / 365 ;
}

//[[Rcpp::export]]
double f_t( double t, List ns_parmlist){
	double tt = std::min(2015.0, tdays2calyears(t));
	//~ double tt = tdays2calyears(t);
	double beta = (double)ns_parmlist["beta"];
	double fB = (double)ns_parmlist["fB"]; 
	double fM = (double)ns_parmlist["fM"]; 
	if (t < EXPANSIONTIME) return f_t( EXPANSIONTIME, ns_parmlist);
	return beta / ( 1 + exp(fB * ( tt - fM) ) ); //note fB strictly positive, transm rate should decrease
}

//[[Rcpp::export]]
double f_t_msm( double t, List ns_parmlist){
	double tt = std::min(2015.0, tdays2calyears(t));
	//~ double tt =  tdays2calyears(t);
	double betamsm = (double)ns_parmlist["betamsm"];
	double fBmsm = (double)ns_parmlist["fBmsm"]; 
	double fMmsm = (double)ns_parmlist["fMmsm"]; 
	if (t < EXPANSIONTIME) return f_t_msm( EXPANSIONTIME, ns_parmlist);
	return betamsm / ( 1 + exp(fBmsm * ( tt - fMmsm) ) ); 
}

//[[Rcpp::export]]
double source_size( double t, List parms){
	return SOURCESIZE * exp(-(double)parms["sourceGrowthRate"] * (MAXSAMPLETIME - t));
}

double diagnosis_rate( double t, List parms){
	double tt  = tdays2calyears(t );
	if ( tt > TINTERVENTIONSTART_YRS ){
		double dr0 = diagnosis_rate( MAXSAMPLETIME, parms);
		return dr0 + (DRFIN - dr0) * std::min(1., (tt - TINTERVENTIONSTART_YRS)/ (TINTERVENTIONFIN_YRS - TINTERVENTIONSTART_YRS));
	}
	return  (double)parms["dK"]/(( 1 + exp(-(double)parms["dB"]*(tt-(double)parms["dM"]) ) ));
}

NumericVector treatment_rate_vector(double t, List parms) {
	double tt  = tdays2calyears(t );
	double treRateIntervention = (double)parms["treRateIntervention"];
	double treRate = (double)parms["treRate"];
	NumericVector trv(5);
	if (tt >= TINTERVENTIONSTART_YRS){
		double tr1 = treRate + (TRFIN - treRate) * std::min(1., (tt - TINTERVENTIONSTART_YRS) / (TINTERVENTIONFIN_YRS-TINTERVENTIONSTART_YRS));
		trv = NumericVector::create( tr1, tr1, tr1, tr1, tr1);
	} else if (tt > 2011){
		trv = NumericVector::create( 0., 0., 0., treRate, treRate);
	} else if (tt > 2004 && tt <= 2011){
		trv = NumericVector::create( 0., 0., 0., 0., treRate);
	} 
	return trv;
}

double  non_adherence_rate( double t){
	double tt  = tdays2calyears(t );
	if (tt >= TINTERVENTIONSTART_YRS){
		return T2D_RATE + (T2D_RATE_FIN - T2D_RATE) * std::min(1., (tt - TINTERVENTIONSTART_YRS) / (TINTERVENTIONFIN_YRS-TINTERVENTIONSTART_YRS));
	}
	return T2D_RATE;
}

NumericVector rowSums( NumericMatrix A){
	NumericVector rval(A.nrow());
	for (int i =0 ; i < A.nrow(); i++){
		for (int j = 0; j < A.ncol(); j++){
			rval[i]+= A(i,j);
		}
	}
	return rval;
}
NumericVector colSums( NumericMatrix A){
	NumericVector rval(A.ncol());
	for (int i =0 ; i < A.ncol(); i++){
		for (int j = 0; j < A.nrow(); j++){
			rval[i]+= A(j,i);
		}
	}
	return rval;
}


//[[Rcpp::export]]
NumericMatrix F_from_state(double t, NumericVector y, List parms){
	// note deme names need to be passed in parms
	int m = (int)parms["m"];
	CharacterVector DEMES = as<CharacterVector>( parms["DEMES"]);
	NumericVector imdemes = as<NumericVector>( parms["imdemes"]);
	NumericVector ifdemes = as<NumericVector>( parms["ifdemes"]);
	NumericVector imsmdemes = as<NumericVector>( parms["imsmdemes"]);
	NumericVector itdemes = as<NumericVector>( parms["itdemes"]);
	NumericVector iddemes = as<NumericVector>( parms["iddemes"]);
	//~ CharacterVector MDEMES = as<CharacterVector>( parms["MDEMES"]);
	//~ CharacterVector MSMDEMES = as<CharacterVector>( parms["MSMDEMES"]);
	//~ CharacterVector FDEMES = as<CharacterVector>( parms["FDEMES"]);
	double eps = (double)parms["eps"];
	double assortprob = (double)parms["assortprob"];
	double assortprob_msm = (double)parms["assortprob_msm"];
	double epsilonT = (double)parms["epsilonT"];
	double epsilonD = (double)parms["epsilonD"];
	NumericVector w = as<NumericVector>(parms["w"]);
	
	NumericMatrix F(m, m);
	List dimnms = List::create(DEMES, DEMES);
	F.attr("dimnames") = dimnms;
	double Im = 0.;
	for (int i = 0; i < imdemes.size(); i++){
		int k = (int)imdemes[i];
		Im += y[k];
	}
	double Imsm = 0.;
	for (int i = 0; i < imsmdemes.size(); i++){
		int k = (int)imsmdemes[i];
		Imsm += y[k];
	}
	double If = 0.;
	for (int i = 0; i < ifdemes.size(); i++){
		int k = (int)ifdemes[i];
		If += y[k];
	}
	y = pmax(y, 0.);
	
	double Nf = y["Sf"] + If;
	double Nm = y["Sm"] + Im;
	double Nmsm = y["Smsm"] + Imsm;
	double Sm = y["Sm"];
	double Sf = y["Sf"];
	double Smsm = y["Smsm"];
	double N = Nf + Nm + Nmsm;
	double I = If + Im + Imsm;
	
	double foi = f_t(t, parms);
	double foi_msm = f_t_msm( t, parms);
	
	
	
	NumericVector f_m2f(w.size());
	for (int i = 0; i < imdemes.size(); i++){
		int k = imdemes[i];
		f_m2f[i] =  foi * w[i] * (Sf  / Nf) * y[k] ;
	}
	
	NumericVector f_f2m(w.size());
	double x =  eps * foi *(assortprob * Sm  / (Nm) +  (1-assortprob)*Sm/(Nm+Nmsm) ); //(Sf  / Nf) *
	for (int i = 0; i < ifdemes.size(); i++){
		int k = ifdemes[i];
		f_f2m[i] =  x * w[i] *  y[k];
	}
	
	NumericVector f_f2msm(w.size());
	x =  eps * foi * (1-assortprob) *  Smsm  / (Nmsm+Nm);
	for (int i = 0; i < ifdemes.size(); i++){
		int k = ifdemes[i];
		f_f2msm[i] =  x * w[i] *  y[k];
	}
	
	NumericVector f_msm2msm(w.size());
	x =  foi_msm * (assortprob_msm * Smsm / Nmsm + (1-assortprob_msm) * Smsm / (Nmsm + Nf ) );
	for (int i = 0; i < ifdemes.size(); i++){
		int k = imsmdemes[i];
		f_msm2msm[i] =  x * w[i] *  y[k];
	}
	
	NumericVector f_msm2f(w.size());
	x =  foi_msm * (1-assortprob_msm) * Sf / (Nmsm + Nf );
	for (int i = 0; i < ifdemes.size(); i++){
		int k = imsmdemes[i];
		f_msm2f[i] =  x * w[i] *  y[k];
	}
	
	
	// fill in F values
	for (int i = 0; i < imdemes.size(); i++){
		int k = imdemes[i];
		F(k, ifdemes[0]) = f_m2f[i];
		
		k = ifdemes[i];
		F(k, imdemes[0]) = f_f2m[i];
		F(k, imsmdemes[0]) = f_f2msm[i];
		
		k = imsmdemes[i];
		F(k, imsmdemes[0]) = f_msm2msm[i];
		F(k, ifdemes[0]) = f_msm2f[i];
	}
	
	// adjust treatment and diagnosis values
	for (int i = 0; i < itdemes.size(); i++){
		int k = itdemes[i];
		F(k, imdemes[0]) *= epsilonT;
		F(k, ifdemes[0]) *= epsilonT;
		F(k, imsmdemes[0]) *= epsilonT;
	}
	for (int i = 0; i < iddemes.size(); i++){
		int k = iddemes[i];
		F(k, imdemes[0]) *= epsilonD;
		F(k, ifdemes[0]) *= epsilonD;
		F(k, imsmdemes[0]) *= epsilonD;
	}
	
	//source
	F(m-1, m-1) =  source_size(t, parms) * (double)parms["sourceGrowthRate"];
	
	return F;
}


//[[Rcpp::export]]
NumericMatrix G_from_state(double t, NumericVector y, List parms){
	int m = (int)parms["m"];
	CharacterVector DEMES = as<CharacterVector>( parms["DEMES"]);
	//~ NumericVector imdemes = as<NumericVector>( parms["imdemes"]);
	//~ NumericVector ifdemes = as<NumericVector>( parms["ifdemes"]);
	//~ NumericVector imsmdemes = as<NumericVector>( parms["imsmdemes"]);
	//~ NumericVector itdemes = as<NumericVector>( parms["itdemes"]);
	//~ NumericVector iddemes = as<NumericVector>( parms["iddemes"]);
	//~ double epsilonT = (double)parms["epsilonT"];
	//~ double epsilonD = (double)parms["epsilonD"];
	double treBreak = (double)parms["treBreak"];
	double gamma0 = (double)parms["gamma0"];
	double gamma1 = (double)parms["gamma1"];
	double gamma2 = (double)parms["gamma2"];
	double gamma3 = (double)parms["gamma3"];
	//~ gammas13 <-  unname( unlist(  parms[c('gamma1', 'gamma2', 'gamma3') ]  ))
	NumericVector gammas13 = NumericVector::create(gamma1, gamma2, gamma3);
	double pss1 = (double)parms["pstartstage1"];
	double pss2 = (double)parms["pstartstage2"];
	double pss3 = (double)parms["pstartstage3"];
	double pss4 = (double)parms["pstartstage4"];
	NumericVector p = NumericVector::create(pss1, pss2, pss3, pss4);
	
	NumericVector i_um = as<NumericVector>( parms["i_um"]);
	NumericVector i_dm = as<NumericVector>( parms["i_dm"]);
	NumericVector i_tm = as<NumericVector>( parms["i_tm"]);
	
	NumericVector i_umsm = as<NumericVector>( parms["i_umsm"]);
	NumericVector i_dmsm = as<NumericVector>( parms["i_dmsm"]);
	NumericVector i_tmsm = as<NumericVector>( parms["i_tmsm"]);
	
	NumericVector i_uf = as<NumericVector>( parms["i_uf"]);
	NumericVector i_df = as<NumericVector>( parms["i_df"]);
	NumericVector i_tf = as<NumericVector>( parms["i_tf"]);
	
	NumericMatrix G(m, m);
	List dimnms = List::create(DEMES, DEMES);
	G.attr("dimnames") = dimnms;
	
	// stage progression from zero
	for (int i = 1; i < 5; i++){
		G( i_um[0], i_um[i]) = y[i_um[0]] * p[i-1] * gamma0; // TODO double check that Y is in same order as G names
		G( i_dm[0], i_dm[i]) = y[i_dm[0]] * p[i-1] * gamma0; 
		G( i_tm[0], i_tm[i]) = y[i_tm[0]] * p[i-1] * gamma0; 
		
		G( i_uf[0], i_uf[i]) = y[i_uf[0]] * p[i-1] * gamma0; 
		G( i_df[0], i_df[i]) = y[i_df[0]] * p[i-1] * gamma0; 
		G( i_tf[0], i_tf[i]) = y[i_tf[0]] * p[i-1] * gamma0; 
		
		G( i_umsm[0], i_umsm[i]) = y[i_umsm[0]] * p[i-1] * gamma0; 
		G( i_dmsm[0], i_dmsm[i]) = y[i_dmsm[0]] * p[i-1] * gamma0; 
		G( i_tmsm[0], i_tmsm[i]) = y[i_tmsm[0]] * p[i-1] * gamma0; 
	}
	
	// stage progression 1-3
	for (int stage = 1; stage < 4; stage++){
		G( i_um[stage], i_um[stage+1]) = gammas13[stage-1] * y[i_um[stage]];
		G( i_dm[stage], i_dm[stage+1]) = gammas13[stage-1] * y[i_dm[stage]];
		G( i_tm[stage], i_tm[stage+1]) = gammas13[stage-1] * y[i_tm[stage]] * treBreak;
		
		G( i_uf[stage], i_uf[stage+1]) = gammas13[stage-1] * y[i_uf[stage]];
		G( i_df[stage], i_df[stage+1]) = gammas13[stage-1] * y[i_df[stage]];
		G( i_tf[stage], i_tf[stage+1]) = gammas13[stage-1] * y[i_tf[stage]] * treBreak;
		
		G( i_umsm[stage], i_umsm[stage+1]) = gammas13[stage-1] * y[i_umsm[stage]];
		G( i_dmsm[stage], i_dmsm[stage+1]) = gammas13[stage-1] * y[i_dmsm[stage]];
		G( i_tmsm[stage], i_tmsm[stage+1]) = gammas13[stage-1] * y[i_tmsm[stage]] * treBreak;
	}
	
	// diagnosis
	double dr = diagnosis_rate( t, parms );
	for (int stage = 0; stage <5; stage++){
		G( i_um[stage], i_dm[stage]) = dr * y[i_um[stage]];
		G( i_uf[stage], i_df[stage]) = dr * y[i_uf[stage]];
		G( i_umsm[stage], i_dmsm[stage]) = dr * y[i_umsm[stage]];
	}
	
	//treatment
	NumericVector trv = treatment_rate_vector(t, parms );
	for (int stage = 0; stage < 5; stage++){
		G( i_dm[stage], i_tm[stage]) = trv[stage] * y[i_dm[stage]];
		G( i_df[stage], i_tf[stage]) = trv[stage] * y[i_df[stage]];
		G( i_dmsm[stage], i_tmsm[stage]) = trv[stage] * y[i_dmsm[stage]];
	}
	
	//non adherence 
	double t2drate = non_adherence_rate( t);
	for (int stage = 0; stage < 5; stage++){
		G( i_tm[stage], i_dm[stage]) = t2drate * y[i_tm[stage]];
		G( i_tf[stage], i_df[stage]) = t2drate * y[i_tf[stage]];
		G( i_tmsm[stage], i_dmsm[stage]) = t2drate * y[i_tmsm[stage]];
	}
	
	//imports 
	double ir = (double)parms["immigrationRate"];
	if (t < EXPANSIONTIME){
		for (int i = 0; i < m-1; i++){
			G(m-1, i) = ir * y[i] * 10;
		}
	} else{
		for (int i = 0; i < m-1; i++){
			G(m-1, i) = ir * y[i];
		}
	}
	
	return G;
}


//[[Rcpp::export]]
List dydt(double t, NumericVector y, List parms){
	int m = (int)parms["m"];
	double pmsm = (double)parms["pmsmFin"];
	double gamma4 = (double)parms["gamma4"];
	double treBreak = (double)parms["treBreak"];
	CharacterVector DEMES = as<CharacterVector>( parms["DEMES"]);
	NumericVector imdemes = as<NumericVector>( parms["imdemes"]);
	NumericVector ifdemes = as<NumericVector>( parms["ifdemes"]);
	NumericVector imsmdemes = as<NumericVector>( parms["imsmdemes"]);
	NumericVector itdemes = as<NumericVector>( parms["itdemes"]);
	NumericVector iddemes = as<NumericVector>( parms["iddemes"]);
	NumericVector i_AnT = as<NumericVector>( parms["i_AnT"]);
	NumericVector i_AT = as<NumericVector>( parms["i_AT"]);
	
	double Im = 0.;
	for (int i = 0; i < imdemes.size(); i++){
		int k = (int)imdemes[i];
		Im += y[k];
	}
	double Imsm = 0.;
	for (int i = 0; i < imsmdemes.size(); i++){
		int k = (int)imsmdemes[i];
		Imsm += y[k];
	}
	double If = 0.;
	for (int i = 0; i < ifdemes.size(); i++){
		int k = (int)ifdemes[i];
		If += y[k];
	}
	y = pmax(y, 0.);
	
	double Nf = y["Sf"] + If;
	double Nm = y["Sm"] + Im;
	double Nmsm = y["Smsm"] + Imsm;
	double N = Nf + Nm + Nmsm;
	//~ double Sm = y["Sm"];
	//~ double Sf = y["Sf"];
	//~ double Smsm = y["Smsm"];
	//double I = If + Im + Imsm;
	
	NumericVector dy(y.size());
	dy.attr("names") = y.attr("names");
	
	if (t < EXPANSIONTIME){
		return List::create(dy);
	}
	
	NumericMatrix F = F_from_state( t, y, parms);
	NumericMatrix G = G_from_state( t, y, parms );
	
	NumericVector csF = colSums(F);
	NumericVector rsG = rowSums(G);
	NumericVector csG = colSums(G);
	
	for (int i = 0; i < m-1; i++){
		dy[i] += csF[i] + csG[i] - rsG[i];
	}
	
	double AnT = 0.;
	for (int i =0 ; i < i_AnT.size(); i++){
		AnT += y[i_AnT[i]];
	}
	double AT = 0.;
	for (int i =0 ; i < i_AT.size(); i++){
		AT += y[i_AT[i]];
	}
	double b = MU + log(N_t(t+1.) / N_t(t) ) + gamma4*AnT/N + gamma4 * treBreak * AT / N;
	
	
	//susceptibles
	double xsm  = 0.;
	for (int i = 0; i < imdemes.size(); i++){
		xsm += csF[imdemes[i]];
	}
	dy["Sm"] = b * (1 - pmsm) * N / 2. - MU * y["Sm"] - xsm ;
	
	double xsmsm  = 0.;
	for (int i = 0; i < imsmdemes.size(); i++){
		xsmsm += csF[imsmdemes[i]];
	}
	dy["Smsm"] = b * pmsm * N / 2. - MU * y["Smsm"] - xsmsm ;
	
	double xf  = 0.;
	for (int i = 0; i < ifdemes.size(); i++){
		xf += csF[ifdemes[i]];
	}
	dy["Sf"] = b * N / 2. - MU * y["Sf"] - xf ;
	
	//natural mort
	for (int i = 0; i < m - 1 ; i++){
		dy[i] -= MU * y[i]; 
	}
	
	//aids death
	for (int i = 0; i < i_AnT.size(); i++){
		dy[ i_AnT[i]] -= gamma4 * y[i_AnT[i]];
	}
	for (int i =0 ; i < i_AT.size(); i++){
		dy[ i_AT[i]] -= gamma4 * y[i_AT[i]] * treBreak;
	}
	
	return List::create(dy);
}



