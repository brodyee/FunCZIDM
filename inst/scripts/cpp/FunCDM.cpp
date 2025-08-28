/*
 Author: Brody Erlandson
 
 This program implements an MCMC sampler for a functional concurrent 
 Dirichlet-multinomial (FunC-DM) regression model. 
*/

//[[Rcpp::plugins(cpp17)]]
//[[Rcpp::depends(RcppArmadillo)]] 
#include <RcppArmadillo.h>

#include "helper.h"

using namespace Rcpp;
using namespace arma; 

//[[Rcpp::export]]
List FunCDMSampler(const int ITER, const arma::umat COUNTS, const arma::mat X, 
                   const arma::uvec ID_END_IDX,
                   const arma::uvec RAND_EFF_COLS,
                   const int BURN_IN = 0,
                   const int NUM_THIN = 1,
                   const int ADJ_FREQ = 250,
                   const double PROPOSAL_CAP = 0.75,
                   const bool ADJ_PROPOSALS = true,
                   const bool RETURN_BURN_IN = false,
                   const bool CAP_PROPOSALS = false,
                   const bool PRINT_PROGRESS = true,
                   const Rcpp::CharacterVector TO_RETRUN 
                                                    = Rcpp::CharacterVector(),
                   Nullable<arma::mat> betaInitial = R_NilValue,
                   Nullable<arma::mat> rInitial = R_NilValue,
                   Nullable<List> priors = R_NilValue,
                   Nullable<List> proposalVars = R_NilValue) { 
  
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Initialization /////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  
  // ---- Data setup ----
  const int NUM_CAT = COUNTS.n_cols;
  const int NUM_OBS = COUNTS.n_rows;
  const int NUM_COEF = X.n_cols;
  const int NUM_RE_PER_ID = RAND_EFF_COLS.n_elem;
  const sp_mat Y = makeBlockMat(X, RAND_EFF_COLS, ID_END_IDX);
  const List PRIORS = (priors.isNotNull()) ? as<List>(priors) : List();
  const List PROPOSAL_VARS = (proposalVars.isNotNull()) ? as<List>(proposalVars)
                                                        : List();
  
  // ---- Saving Samples setup ----
  const int numSamples = RETURN_BURN_IN ? ITER/NUM_THIN :
                                         (ITER - BURN_IN)/NUM_THIN; 
                                         //  floor division (int/int)
  
  // Theses you must specify to save the samples
  mat uSamples, phiSamples, tauSamples, xiSamples, kappaSamples;
  cube cSamples, lambdaSamples, nuSamples, rSamples;
  bool saveC = isValueInCharVector(TO_RETRUN, "c");
  bool saveU = isValueInCharVector(TO_RETRUN, "u");
  bool saveR = isValueInCharVector(TO_RETRUN, "r");
  bool savePhi = isValueInCharVector(TO_RETRUN, "phi");
  bool saveKappa = isValueInCharVector(TO_RETRUN, "kappa");
  bool saveLambda = isValueInCharVector(TO_RETRUN, "lambda");
  bool saveNu = isValueInCharVector(TO_RETRUN, "nu");
  bool saveTau = isValueInCharVector(TO_RETRUN, "tau");
  bool saveXi = isValueInCharVector(TO_RETRUN, "xi");
  bool saveRA = isValueInCharVector(TO_RETRUN, "RA");
  if (saveC || saveRA) cSamples.set_size(NUM_OBS, NUM_CAT, numSamples); 
  if (saveU) uSamples.set_size(NUM_OBS, numSamples);
  if (saveR) rSamples.set_size(Y.n_cols, NUM_CAT, numSamples);
  if (savePhi) phiSamples.set_size(NUM_CAT, numSamples);
  if (saveKappa) kappaSamples.set_size(NUM_CAT, numSamples);
  if (saveLambda) lambdaSamples.set_size(NUM_COEF - 1, NUM_CAT, numSamples);
  if (saveNu) nuSamples.set_size(NUM_COEF - 1, NUM_CAT, numSamples);
  if (saveTau) tauSamples.set_size(numSamples, NUM_CAT);
  if (saveXi) xiSamples.set_size(numSamples, NUM_CAT);
  
  // These are always saved
  mat betaAcceptProp(NUM_COEF, NUM_CAT, fill::zeros),
      rAcceptProp(Y.n_cols, NUM_CAT, fill::zeros),
      rMeans(Y.n_cols, NUM_CAT, fill::zeros);
  cube betaSamples(NUM_COEF, NUM_CAT, numSamples);
  
  // ---- Prior Specification ----
  // Only the first beta variance, phi hyperparameters, and theta 
  // hyperparameters can be specified. 
  double firstBetaSD = 1., // constant component of the intercepts prior var
         phiShape = 3.,
         phiRate = 9., 
         kappaShape = 100.,
         kappaRate = 900.;
  if (PRIORS.containsElementNamed("first beta sd"))
    firstBetaSD = as<double>(PRIORS["first beta sd"]);
  if (PRIORS.containsElementNamed("kappaShape"))
    kappaShape = as<double>(PRIORS["kappaShape"]); 
  if (PRIORS.containsElementNamed("kappaRate"))
    kappaRate = as<double>(PRIORS["kappaRate"]); 
  // the individual and taxon specific covariate variance (phi) hyperparameters
  if (PRIORS.containsElementNamed("a"))
    phiShape = as<double>(PRIORS["a"]);
  if (PRIORS.containsElementNamed("b"))
    phiRate = as<double>(PRIORS["b"]);

  // ---- Proposal variances ----
  // The proposal variances for beta and r can be specified.
  // If proposal adjustment is on, this is the initial proposal variance for 
  // beta.
  double betaPropSD = .3,
         rProposalSD = 1.;
  if (PROPOSAL_VARS.containsElementNamed("beta proposal sd"))
    betaPropSD = as<double>(PROPOSAL_VARS["beta proposal sd"]);
  if (PROPOSAL_VARS.containsElementNamed("r proposal sd"))
    rProposalSD = as<double>(PROPOSAL_VARS["r proposal sd"]);
  
  // For proposal adjustment
  mat betaProposalSD(NUM_COEF, NUM_CAT, fill::value(betaPropSD));
  Mat<short> betaAcceptCount(NUM_COEF, NUM_CAT, fill::zeros);
  
  // ---- Initialize parameters ----
  // ordered by model hierarchy
  mat c(NUM_OBS, NUM_CAT, fill::ones);
  vec u(NUM_OBS);
  mat gamma(NUM_OBS, NUM_CAT);
  mat beta(NUM_COEF, NUM_CAT);
  mat r(Y.n_cols, NUM_CAT);
  vec phi(NUM_CAT);
  vec kappa(NUM_CAT, fill::ones);
  mat lambda(NUM_COEF - 1, NUM_CAT, fill::ones);
  mat nu(NUM_COEF - 1, NUM_CAT, fill::ones);
  rowvec tau(NUM_CAT, fill::ones);
  rowvec xi(NUM_CAT, fill::ones);
  
  // The initial values for beta and r can be specified.
  if (betaInitial.isNotNull()) {
    beta = as<mat>(betaInitial);
  } else {
    fillRandSampleMat(beta, -.75, .75);
  }
  if (rInitial.isNotNull()) {
    r = as<mat>(rInitial);
  } else {
    fillRandSampleMat(r, -.05, .05);
  }
  fillGammaVec(phi, phiShape, phiRate);
  fillLambdaNu(lambda, nu);
  
  gamma = arma::exp(X*beta + Y*r); sampleCnoZI(c, COUNTS, u, gamma);
  
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// Sampling ////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  int s = 0,
      printFreq = std::ceil(ITER/20.);
  for (int i = 0; i < ITER; i++) {
    sampleCnoZI(c, COUNTS, u, gamma);
    sampleU(u, c, COUNTS);
    sampleBetaNoZI(beta, X, c, kappa, tau, lambda, gamma, betaAcceptProp, 
                   betaProposalSD, betaAcceptCount, firstBetaSD);
    sampleKappa(kappa, beta, kappaShape, kappaRate);
    sampleRGroupedNoZI(r, Y, c, phi, gamma, NUM_RE_PER_ID, ID_END_IDX,
                       rAcceptProp, rProposalSD);
    samplePhi(phi, r, NUM_RE_PER_ID, phiShape, phiRate);
    sampleLambda(lambda, beta, nu, tau);
    sampleNu(nu, lambda);
    sampleTau(tau, beta, lambda, xi);
    sampleXi(xi, tau);
    
    // Adjusting proposal variances
    if (ADJ_PROPOSALS && (i + 1) % ADJ_FREQ == 0 && i < BURN_IN)
      adjustProposalVars(betaProposalSD, betaAcceptCount, ADJ_FREQ,
                         CAP_PROPOSALS, PROPOSAL_CAP);
    
    // Saving
    if (i == BURN_IN - 1) {
      betaAcceptProp.zeros(); rAcceptProp.zeros();
    }
    if ((i >= BURN_IN || RETURN_BURN_IN) && i % NUM_THIN == 0) {
      betaSamples.slice(s) = beta;
      rMeans += r;
      if (saveC || saveRA) cSamples.slice(s) = c;
      if (saveU) uSamples.col(s) = u;
      if (saveR) rSamples.slice(s) = r;
      if (savePhi) phiSamples.col(s) = phi;
      if (saveKappa) kappaSamples.col(s) = kappa;
      if (saveLambda) lambdaSamples.slice(s) = lambda;
      if (saveNu) nuSamples.slice(s) = nu;
      if (saveTau) tauSamples.row(s) = tau;
      if (saveXi) xiSamples.row(s) = xi;
      s++;
    }
    
    // Progress
    if (PRINT_PROGRESS && (i + 1) % printFreq == 0)
      printProgress(i+1, ITER);
  }
  
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// Output //////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  betaAcceptProp /= (ITER - BURN_IN); rAcceptProp /= (ITER - BURN_IN);
  rMeans /= numSamples;


  List output = List::create(Named("beta") = (cube) betaSamples,
                             Named("betaAcceptProp") = betaAcceptProp,
                             Named("rAcceptProp") = rAcceptProp,
                             Named("rMeans") = rMeans,
                             Named("etaMeanPropZeros") = 0.0); // no eta in this
                                                               // model
  if (saveC) output["c"] = cSamples;
  if (saveRA) {
    scaleCubeRows(cSamples);
    output["RA"] = cSamples; 
  }
  if (saveU) output["u"] = uSamples;
  if (saveR) output["r"] = rSamples;
  if (saveKappa) output["kappa"] = kappaSamples;
  if (savePhi) output["phi"] = phiSamples;
  if (saveLambda) output["lambda"] = lambdaSamples;
  if (saveNu) output["nu"] = nuSamples;
  if (saveTau) output["tau"] = tauSamples;
  if (saveXi) output["xi"] = xiSamples;
                          
  return output;
}
