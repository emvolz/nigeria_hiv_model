#Applying Phylodynamic Analysis of Linkage Between Key and General Populations to Inform Prevention Efforts in Mixed HIV Epidemic Setting.

--

Code & packages in support of submitted manuscript: 
Applying Phylodynamic Analysis of Linkage Between Key and General Populations to Inform Prevention Efforts in Mixed HIV Epidemic Setting.
Erik Volz, Nicaise Ndembi, Rebecca Nowak, Gustavo Kijak, John Idoko, Patrick Dakum, Walter Royal, Stefan Baral, Mark Dybul, William A. Blattner, and Man Charurat

This file provides functions for 
* 1. simulating HIV epidemic model incl. MSM & het males and females
* 2. computing priors, reducing dimension of model outputs
* 3. computing prob of nongenetic surveillance data (dprior.sim) 
* 4. computing prob of tree data using rcolgem package

Requirements:
* 1. rcolgem package for phylodynamics: http://colgem.r-forge.r-project.org/
* 2. gfmodel11, compiled package for quickly simulating HIV model
Other dependencies: Rcpp, RcppArmadillo, ape, deSolve

NOTE: Running this file unmodified will require data that is not included 
* 1.  'mergeTrees0.nwk' : sample of BEAST trees
* 2.  'sampleStates0.csv' : initial state of sampled patients
if interested in using/adapting this model, or accessing full data set, contact: 
* Erik Volz <evolz@imperial.ac.uk>
* Man Charurat <MCharurat@ihv.umaryland.edu>

