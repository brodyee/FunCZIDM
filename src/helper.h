/*
 Author: Brody Erlandson
 Last Modified: 3/3/2025
 
 This is a helper file for the ZIDMwRE.cpp file. It contains the
 functions that are used within the sampler.
 */

#ifndef HELPER_H
#define HELPER_H

// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma; 

////////////////////////////////////////////////////////////////////////////////
//////////////////// Initialization & output helpers ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

sp_mat makeBlockMat(const mat& X, const uvec& RAND_EFF_COLS,
                    const uvec& ID_END_IDX) { 
  const int NUM_RE_COEF = ID_END_IDX.n_elem*RAND_EFF_COLS.n_elem,
            NUM_RE_PER_ID = RAND_EFF_COLS.n_elem,
            NUM_OBS = X.n_rows;
  unsigned int startIdx, endIdx, startCol, endCol;
  mat blockMatItems = X.cols(RAND_EFF_COLS);
  sp_mat blockMat(NUM_OBS, NUM_RE_COEF);

  startIdx = 0;
  for (int id = 0; id < ID_END_IDX.n_elem; id++) {
    endIdx = ID_END_IDX[id];
    startCol = id*NUM_RE_PER_ID;
    endCol = startCol + NUM_RE_PER_ID - 1;
    blockMat.submat(startIdx, startCol,
                    endIdx, endCol) = blockMatItems.rows(startIdx, endIdx);
    startIdx = endIdx + 1;
  }  

  return blockMat;
}

bool isValueInCharVector(Rcpp::CharacterVector vec, std::string value) {
  for (int i = 0; i < vec.size(); i++) {
      if (vec[i] == value) {
          return true;
      }
  }
  return false;
}

void scaleCubeRows(cube& toScale) {
  const int NUM_ROW = toScale.n_rows,
            NUM_COL = toScale.n_cols,
            NUM_SLICE = toScale.n_slices;
  cube rowSums = arma::sum(toScale, 1);

  for (int slice = 0; slice < NUM_SLICE; slice++) {
    for (int row = 0; row < NUM_ROW; row++) {
      for (int col = 0; col < NUM_COL; col++) {
        toScale.at(row, col, slice) = toScale.at(row, col, slice)
                                      / rowSums.at(row, 0, slice);
      }
    }
  }
}

void fillRandSampleMat(mat& mat, const double& a=-1, const double& b=1) {
  int n = mat.n_rows*mat.n_cols;

  for (int i = 0; i < n; i++) {
    mat[i] = R::runif(a, b);
  }
}

void fillGammaVec(vec& vec, const double& alpha = 3, const double& beta = 9) {
  int n = vec.n_elem;

  for (int i = 0; i < n; i++) {
    vec[i] = R::rgamma(alpha, beta);
  }
}

void fillLambdaNu(mat& lambda, mat& nu) {
  const int NUM_COEF = lambda.n_rows,
            NUM_CAT = lambda.n_cols;

  for (int coef = 0; coef < NUM_COEF; coef++) {
    for (int cat = 0; cat < NUM_CAT; cat++) {
      nu.at(coef, cat) = 1/R::rgamma(1/2, 1);
      lambda.at(coef, cat) = std::sqrt(1/R::rgamma(1/2, nu.at(coef, cat)));
    }
  }
}

void fillEta(Mat<short>& eta, const umat& COUNTS, const uvec ID_END_IDX) {
  int NUM_OBS = ID_END_IDX.n_elem,
      NUM_CAT = eta.n_cols;
  unsigned int startIdx, endIdx;
  
  for (int cat = 0; cat < NUM_CAT; cat++) {
    startIdx = 0;
    for (int obs = 0; obs < NUM_OBS; obs++) {
      endIdx = ID_END_IDX[obs];
      if (!COUNTS.col(cat).subvec(startIdx, endIdx).is_zero())
        eta.col(cat).subvec(startIdx, endIdx).fill(1);
      startIdx = endIdx + 1;
    }
  }   
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Sampling helpers ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void adjustProposalVars(mat& betaProposalSD, Mat<short>& betaAcceptCount,
                        const double freq, const bool capVar = false,
                        const double cap = .75) {
  const int NUM_COEF = betaProposalSD.n_rows,
            NUM_CAT = betaProposalSD.n_cols;
  
  double acceptProp;
  for (int coef = 0; coef < NUM_COEF; coef++) {
    for (int cat = 0; cat < NUM_CAT; cat++) {
      acceptProp = ((double) betaAcceptCount.at(coef, cat))/freq;
      if (acceptProp <= .15) {
        betaProposalSD.at(coef, cat) *= acceptProp/.25;
        betaAcceptCount.at(coef, cat) = 0;
        if (betaProposalSD.at(coef, cat) <= .001)
          betaProposalSD.at(coef, cat) = .001;
      } else if (acceptProp <= .45) {
        betaAcceptCount.at(coef, cat) = 0;
      } else {
        betaProposalSD.at(coef, cat) *= acceptProp/.25;
        betaAcceptCount.at(coef, cat) = 0;
        if (capVar && betaProposalSD.at(coef, cat) > cap)
          betaProposalSD.at(coef, cat) = cap;
      }
    }
  }
}
 
double dbetaBin(const int& x, const int& n, const double& a,
                const double& b, bool log = false) {
  double dens = 0;

  if (log) {
    dens = R::lchoose(n, x) + R::lbeta(x + a, n - x + b) - R::lbeta(a, b);
  } else {
    dens = R::choose(n, x)*R::beta(x + a, n - x + b)/R::beta(a, b);
  }

  return dens;
}

void sampleEta(Mat<short>& eta, mat& c, const mat& gamma, const umat& COUNTS,
               const double alpha, const double beta) {
  int NUM_OBS = eta.n_rows,
      NUM_CAT = eta.n_cols;
  short etaProp;
  double etaPropProb, etaPrevProb, logAcceptProb, logZratio, cProp;
  const double INF = std::numeric_limits<double>::infinity();
  vec T = arma::sum(c, 1);
  uvec zDot = arma::sum(COUNTS, 1);

  for (int cat = 0; cat < NUM_CAT; cat++) {
    for (int obs = 0; obs < NUM_OBS; obs++) { 
      // eta = 1 -> propose eta = 0 (only if COUNTS = 0)
      if ((eta.at(obs, cat) == 1) && (COUNTS.at(obs, cat) == 0)) { 
        etaProp = 0; cProp = 0;
        etaPropProb = dbetaBin(0, 1, alpha, beta, true); // p(eta = 0)
        etaPrevProb = dbetaBin(1, 1, alpha, beta, true); // p(eta = 1)
        logZratio = zDot[obs]*(std::log(T[obs]) - std::log(T[obs] 
                                                             - c.at(obs, cat)));

        // Q cancels with pC
        logAcceptProb = logZratio + etaPropProb - etaPrevProb;
      } else if ((eta.at(obs, cat) == 0) && (COUNTS.at(obs, cat) == 0)) {
        // eta = 0 -> propose eta = 1 (only if COUNTS = 0)
        etaProp = 1; cProp = R::rgamma(gamma.at(obs, cat), 1);
        etaPropProb = dbetaBin(1, 1, alpha, beta, true); // p(eta = 1)
        etaPrevProb = dbetaBin(0, 1, alpha, beta, true); // p(eta = 0)
        logZratio = zDot[obs]*(std::log(T[obs] + cProp) - std::log(T[obs])); 

        // Q cancels with pC
        logAcceptProb = logZratio + etaPropProb - etaPrevProb;
      } else { //  COUNTS != 0 
        // eta must be 1
        if (eta.at(obs, cat) == 1) {
          // Must not accept, to keep eta = 1
          logAcceptProb = -1*INF; 
        } else { // eta.at(obs, cat) == 0
          // Must make eta = 1, since COUNTS != 0
          logAcceptProb = INF;
          etaProp = 1;
          cProp = R::rgamma(gamma.at(obs, cat), 1);
        }
      }

      if (logAcceptProb > std::log(R::runif(0, 1))) {
        eta.at(obs, cat) = etaProp;
        c.at(obs, cat) = cProp;
      }
    }
  }   
}

void sampleEtaGrouped(Mat<short>& eta, mat& c, const mat& gamma, 
                      const umat& COUNTS, const uvec ID_END_IDX, 
                      const double alpha, const double beta) {
  int NUM_OBS = ID_END_IDX.n_elem,
      NUM_CAT = eta.n_cols;
  unsigned int startIdx, endIdx;
  short etaProp;
  double etaPropProb, etaPrevProb, logAcceptProb, logZratio, cProp, numMeas;
  const double INF = std::numeric_limits<double>::infinity();
  vec T = arma::sum(c, 1);
  uvec zDot = arma::sum(COUNTS, 1);
  
  auto fillcProp = [&gamma](mat& c, int cat, unsigned int startIdx,
                            unsigned int endIdx) {
    for (int i = startIdx; i <= endIdx; i++) {
      c.at(i, cat) = R::rgamma(gamma.at(i, cat), 1);
    }
  };
  
  for (int cat = 0; cat < NUM_CAT; cat++) {
    startIdx = 0;
    for (int obs = 0; obs < NUM_OBS; obs++) {
      endIdx = ID_END_IDX[obs];
      numMeas = endIdx - startIdx + 1;
      
      if ((eta.at(startIdx, cat) == 1) // all idxs are the same for each obs
            && (COUNTS.col(cat).subvec(startIdx, endIdx).is_zero())) { 
        // eta = 1 -> propose eta = 0 (only if COUNTS = 0)
        etaProp = 0; cProp = 0;
        etaPropProb = dbetaBin(0, 1, alpha, beta, true); // p(eta = 0)
        etaPrevProb = dbetaBin(1, 1, alpha, beta, true); // p(eta = 1)
        logZratio = arma::sum(zDot.subvec(startIdx, endIdx)
                                %( arma::log(T.subvec(startIdx, endIdx)) 
                                     - arma::log(T.subvec(startIdx, endIdx) 
                                     - c.col(cat).subvec(
                                         startIdx, endIdx))
                                )
        );
        
        // Q cancels with pC
        logAcceptProb = logZratio + numMeas*(etaPropProb - etaPrevProb);
      } else if ((eta.at(startIdx, cat) == 0) // all idxs are the same 
                   && COUNTS.col(cat).subvec(startIdx, endIdx).is_zero()) {
        // eta = 0 -> propose eta = 1 (only if COUNTS = 0)
        etaProp = 1;
        // fill c with props, change back to 0 if not accepted
        fillcProp(c, cat, startIdx, endIdx);
        etaPropProb = dbetaBin(1, 1, alpha, beta, true); // p(eta = 1)
        etaPrevProb = dbetaBin(0, 1, alpha, beta, true); // p(eta = 0)
        logZratio = arma::sum(zDot.subvec(startIdx, endIdx)
                                %( arma::log(T.subvec(startIdx, endIdx) 
                                               + c.col(cat).subvec(startIdx,
                                                       endIdx)) 
                                     - arma::log(T.subvec(startIdx,endIdx))
                                )
        );
        
        // Q cancels with pC
        logAcceptProb = logZratio + numMeas*(etaPropProb - etaPrevProb);
      } else { // All COUNTS != 0
        // eta must be 1
        if (eta.at(startIdx, cat) == 1) { // all idxs are the same 
          // Must not accept, to keep eta = 1
          logAcceptProb = -1*INF; 
          etaProp = 0; // to make sure we don't enter the else if
        } else { // eta.at(startIdx, cat) == 0
          // Must make eta = 1, since all COUNTS != 0
          logAcceptProb = INF;
          etaProp = 1;
          fillcProp(c, cat, startIdx, endIdx);
        }
      }
      
      if (logAcceptProb > std::log(R::runif(0, 1))) {
        eta.col(cat).subvec(startIdx, endIdx).fill(etaProp);
        if (etaProp == 0) {
          c.col(cat).subvec(startIdx, endIdx).fill(0);
        } // else c already filled
      } else if (etaProp == 1) {
        // etaProp == 1 and not accepted, change back to 0
        c.col(cat).subvec(startIdx, endIdx).fill(0);
      }
      startIdx = endIdx + 1;
    }
  }   
}
 
void sampleC(mat& c, const umat& COUNTS, const vec& u, const mat& gamma,
             const Mat<short>& eta) {
  const int NUM_OBS = COUNTS.n_rows; 
  int row, col;
  arma::uword linearIdx;
  
  c.zeros();
  uvec idx = arma::find(eta != 0);
  vec scales = 1.0 / (1.0 + u);
  mat shapes = arma::conv_to<mat>::from(COUNTS) + gamma;
  for (arma::uword i = 0; i < idx.n_elem; i++) { // only update where eta != 0
    linearIdx = idx[i];
    // Convert linear index to (row, col) indices (column-major storage)
    row = linearIdx % NUM_OBS;
    col = linearIdx / NUM_OBS;
    
    c.at(row, col) = std::max(R::rgamma(shapes.at(row, col), scales[row]),
         std::numeric_limits<double>::epsilon());
  }
}
 
void sampleU(vec& u, const mat& c, const umat& COUNTS){
  const int NUM_ROW = u.n_rows;
  
  arma::vec zDot = arma::conv_to<arma::vec>::from(arma::sum(COUNTS, 1));
  arma::vec T = arma::sum(c, 1);
  for (int row = 0; row < NUM_ROW; row++) {
    u[row] = R::rgamma(zDot[row], 1.0/T[row]);
  }
}

double ldgammaDiffScale(vec x, vec shape1, vec shape2) {
  // sum of the differences in log density of a gamma dist with rate/scale = 1
  return arma::sum((shape1 - 1)%arma::log(x) - arma::lgamma(shape1)
                   - (shape2 - 1)%arma::log(x) + arma::lgamma(shape2));
}
 
void sampleBeta(mat& beta, const mat& X, const mat& c, const vec& kappa,
                const rowvec& tau, const mat& lambda, mat& gamma, 
                mat& betaAcceptProp, const mat& betaProposalSD, 
                Mat<short>& betaAcceptCount, const Mat<short>& eta, 
                const double firstBetaSD) {
  const int NUM_COEF = beta.n_rows,
            NUM_CAT = beta.n_cols,
            NUM_OBS = c.n_rows;
  double betaProp, logAcceptCSum, logAcceptProb, betaCurr, sd, lamTilde;
  vec gammaPropVec(NUM_OBS), gammaCol(NUM_OBS), cCol(NUM_OBS);

  for (int cat = 0; cat < NUM_CAT; cat++) {
    arma::uvec obsEtaOne = arma::find(eta.col(cat) != 0);
    gammaCol = gamma.col(cat);
    cCol = c.col(cat);
    
    for (int coef = 0; coef < NUM_COEF; coef++) {
      lamTilde = ((lambda.at(coef - 1, cat)*lambda.at(coef - 1, cat)
                   *kappa[cat]*kappa[cat])
                  /(kappa[cat]*kappa[cat] 
                    + lambda.at(coef - 1, cat)*lambda.at(coef - 1, cat)
                      *tau[cat]*tau[cat]));
      sd = (coef == 0) ? firstBetaSD : tau[cat]*sqrt(lamTilde);
      betaCurr = beta.at(coef, cat);
      logAcceptCSum = 0; // reset logAcceptCSum
      betaProp = R::rnorm(betaCurr, betaProposalSD.at(coef, cat)); 

      // update gamma
      gammaPropVec = gammaCol%arma::exp(X.col(coef)*(betaProp - betaCurr));
      // acceptance ratio only for obs with eta = 1
      logAcceptCSum = ldgammaDiffScale(cCol.elem(obsEtaOne),
                                       gammaPropVec.elem(obsEtaOne),
                                       gammaCol.elem(obsEtaOne));

      logAcceptProb = logAcceptCSum // density of c with proposal
                      + R::dnorm(betaProp, 0, sd, true) // density of proposal 
                      - R::dnorm(betaCurr, 0, sd, true); // density of current  
      
      if (logAcceptProb > std::log(R::runif(0, 1))) {
        beta.at(coef, cat) = betaProp;
        gammaCol = gammaPropVec;
        betaAcceptProp.at(coef, cat)++;
        betaAcceptCount.at(coef, cat)++;
      }
    }
    gamma.col(cat) = gammaCol; // update gamma with new values
  }
}

void sampleKappa(vec& kappa, const mat& beta, const double& kappaShape,
                 const double& kappaRate) {
  const int NUM_CAT = kappa.n_elem,
            NUM_PARAM = beta.n_rows - 1; // first row is intercept
  mat betaNew = beta.rows(1, NUM_PARAM); // remove intercept
  
  double alphaPrime = kappaShape + NUM_PARAM*.5;
  arma::rowvec betaPrimes = arma::sum(arma::square(betaNew), 0)/2.0 + kappaRate;
  for (int cat = 0; cat < NUM_CAT; cat++) {
    kappa[cat] = std::sqrt(1.0/R::rgamma(alphaPrime, 1.0/betaPrimes[cat]));
  }
}

void sampleR(mat& r, const sp_mat& Y, const mat& c, const vec& phi, mat& gamma,
             const int& NUM_RE_PER_ID, const uvec& ID_END_IDX, mat& rAcceptProp,
             const double& rProposalSD, const Mat<short>& eta) {
  const int NUM_COEF = r.n_rows,
            NUM_CAT = r.n_cols,
            NUM_OBS = c.n_rows;
  unsigned int startIdx, endIdx, id;
  double rProp, logAcceptCSum, logAcceptProb, rCurr;
  vec gammaPropVec(NUM_OBS);
  
  id = 0;
  startIdx = 0;
  endIdx = ID_END_IDX[id];
  for (int coef = 0; coef < NUM_COEF; coef++) {
    vec YSub = vec(Y.col(coef).rows(startIdx, endIdx));
    mat cIdx = c.rows(startIdx, endIdx);
    mat gammaIdx = gamma.rows(startIdx, endIdx);

    for (int cat = 0; cat < NUM_CAT; cat++) {
      uvec idx = arma::find(eta.col(cat).subvec(startIdx, endIdx) != 0);

      // if all eta = 0, then don't need to update r
      if (idx.n_elem == 0) continue;
      
      vec cIdxCol = cIdx.col(cat);
      vec gammaIdxCol = gammaIdx.col(cat);
      rCurr = r.at(coef, cat);
      logAcceptCSum = 0; // reset logAcceptCSum
      rProp = R::rnorm(rCurr, rProposalSD); // get proposals
      
      // acceptance ratio only for obs with eta = 1
      logAcceptCSum = ldgammaDiffScale(cIdxCol.elem(idx),
                                       gammaPropVec.elem(idx),
                                       gammaIdxCol.elem(idx));
      
      logAcceptProb =  logAcceptCSum // density of c with proposal 
                     + R::dnorm(rProp, 0, phi[cat], true) // density of proposal 
                     - R::dnorm(rCurr, 0, phi[cat], true); // density of current 

      if (logAcceptProb > std::log(R::runif(0, 1))) {
        r.at(coef, cat) = rProp;
        gammaIdx.col(cat) = gammaPropVec.subvec(startIdx, endIdx);
        rAcceptProp.at(coef, cat)++;
      }
    }
    gamma.rows(startIdx, endIdx) = gammaIdx;

    // update id when done with all RE for id
    if ((coef + 1) % NUM_RE_PER_ID == 0) { 
      id++;
      startIdx = endIdx + 1;
      endIdx = ID_END_IDX[id];
    }
  }
}

void sampleRGrouped(mat& r, const sp_mat& Y, const mat& c, const vec& phi,
                    mat& gamma, const int& NUM_RE_PER_ID, 
                    const uvec& ID_END_IDX, mat& rAcceptProp,
                    const double& rProposalSD, const Mat<short>& eta) {
  const int NUM_COEF = r.n_rows,
            NUM_CAT = r.n_cols,
            NUM_OBS = c.n_rows;
  unsigned int startIdx, endIdx, id;
  double rProp, logAcceptCSum, logAcceptProb, rCurr;
  vec gammaPropVec(NUM_OBS);
  
  id = 0;
  startIdx = 0;
  endIdx = ID_END_IDX[id];
  for (int coef = 0; coef < NUM_COEF; coef++) {
    vec YSub = vec(Y.col(coef).rows(startIdx, endIdx));
    mat cIdx = c.rows(startIdx, endIdx);
    mat gammaIdx = gamma.rows(startIdx, endIdx);
    
    for (int cat = 0; cat < NUM_CAT; cat++) {
      // for grouped, eta the same for all obs in group
      // if eta = 0, then don't need to update r
      if (eta.at(startIdx, cat) == 0) continue;

      rCurr = r.at(coef, cat);
      logAcceptCSum = 0; // reset logAcceptCSum
      rProp = R::rnorm(rCurr, rProposalSD); // get proposals
      
      // update gamma
      gammaPropVec.subvec(startIdx, endIdx) 
                            = gammaIdx.col(cat)%arma::exp(YSub*(rProp - rCurr));
      // acceptance ratio only for obs with eta = 1
      logAcceptCSum = ldgammaDiffScale(cIdx.col(cat),
                                        gammaPropVec.subvec(startIdx, endIdx),
                                        gammaIdx.col(cat));
      
      logAcceptProb = logAcceptCSum // density of c with proposal 
                     + R::dnorm(rProp, 0, phi[cat], true) // density of proposal 
                     - R::dnorm(rCurr, 0, phi[cat], true); // density of current 

      if (logAcceptProb > std::log(R::runif(0, 1))) {
        r.at(coef, cat) = rProp;
        gammaIdx.col(cat) = gammaPropVec.subvec(startIdx, endIdx);
        rAcceptProp.at(coef, cat)++;
      }
    }
    gamma.rows(startIdx, endIdx) = gammaIdx;
    
    // update id when done with all RE for id
    if ((coef + 1) % NUM_RE_PER_ID == 0) { 
      id++;
      startIdx = endIdx + 1;
      endIdx = ID_END_IDX[id];
    }
  }
}

void samplePhi(vec& phi, const mat& r, const int& NUM_RE_PER_ID, 
                 const double& phiShape, const double& phiRate) {
  const int NUM_CAT = phi.n_elem,
            NUM_ROW = r.n_rows;
  
  double alphaPrime = phiShape + NUM_ROW*NUM_RE_PER_ID*.5; // NUM_RE_PER_ID = 1 
                                                     // until model is exteneded
  arma::rowvec betaPrimes = arma::sum(arma::square(r), 0)/2.0 + phiRate;
  for (int cat = 0; cat < NUM_CAT; cat++) {
    phi[cat] = std::sqrt(1.0/R::rgamma(alphaPrime, 1.0/betaPrimes[cat]));
  }
}

void sampleLambda(mat& lambda, const mat& beta, const mat& nu,
                  const rowvec& tau) {
  const int NUM_COEF = lambda.n_rows,
            NUM_CAT = lambda.n_cols;
  mat betaNew = beta.rows(1, NUM_COEF);
  
  mat betaSq = betaNew%betaNew;
  mat betaPrimes = (1.0/nu) + betaSq.each_row()/(2.0*(tau%tau));
  for (int cat = 0; cat < NUM_CAT; cat++) {
    for (int coef = 0; coef < NUM_COEF; coef++) {
      lambda.at(coef, cat) 
                    = std::sqrt(1.0/R::rgamma(1, 1.0/betaPrimes.at(coef, cat)));                                                       
    }
  }
}

void sampleNu(mat& nu, const mat& lambda) {
  const int NUM_COEF = lambda.n_rows,
            NUM_CAT = lambda.n_cols;
  
  mat betaPrimes = 1.0 + (1.0/(lambda%lambda));
  for (int coef = 0; coef < NUM_COEF; coef++) {
    for (int cat = 0; cat < NUM_CAT; cat++) {
      nu.at(coef, cat) = 1.0/R::rgamma(1, 1.0/betaPrimes.at(coef, cat));
    }
  }
}

void sampleTau(rowvec& tau, const mat& beta, const mat& lambda,
               const rowvec& xi) {
  const int NUM_COEF = lambda.n_rows,
            NUM_CAT = lambda.n_cols;
  mat betaNew = beta.rows(1, NUM_COEF);
  
  double alphaPrime = (NUM_COEF + 1.0)/2.0;
  mat betaLambdas = (betaNew%betaNew)/(2.0*(lambda%lambda));
  rowvec betaPrimes = arma::sum(betaLambdas, 0) + (1.0/xi);
  for (int cat = 0; cat < NUM_CAT; cat++) {
    tau[cat] = std::sqrt(1.0/R::rgamma(alphaPrime, 1.0/betaPrimes[cat]));
  }
}

void sampleXi(rowvec& xi, const rowvec& tau) {
  const int NUM_CAT = xi.n_elem;
  
  rowvec betaPrimes = 1.0 + (1.0/(tau%tau));
  for (int cat = 0; cat < NUM_CAT; cat++) {
    xi[cat] = 1.0/R::rgamma(1, 1.0/betaPrimes[cat]);
  }
}

void printProgress(int i, int ITER) {
  int barWidth = 50;
  double progress = static_cast<double>(i) / ITER;
  Rcout << "\r[";
  int pos = static_cast<int>(barWidth * progress);
  for (int j = 0; j < barWidth; ++j) {
    if (j < pos)
      Rcout << "=";
    else if (j == pos)
      Rcout << ">";
    else
      Rcout << " ";
  }
  Rcout << "] " << static_cast<int>(progress * 100) << "% - " << i 
        << "/" << ITER << " iterations" << std::flush;
  if (i == ITER) {
    Rcout << "\n";
  }
}

////////////////////////////////////////////////////////////////////////////////
///////////////////// Functions for Evaluating Samples /////////////////////////
////////////////////////////////////////////////////////////////////////////////

// ---- Fit related functions ----

//[[Rcpp::export]]
arma::cube getFit(const arma::cube& beta, const arma::mat& basis,
                  const arma::vec& covariates, const arma::vec& covIdxs) {
  const int NUM_PARAM = beta.n_rows,
            NUM_CAT = beta.n_cols,
            NUM_SAMPLE = beta.n_slices,
            NUM_TESTPOINTS = basis.n_rows,
            NUM_DF = basis.n_cols;
  mat Xtv(NUM_TESTPOINTS, NUM_PARAM),
      ones(NUM_TESTPOINTS, 1, arma::fill::ones);
  cube fit(NUM_TESTPOINTS, NUM_CAT, NUM_SAMPLE);

  int startIdx = 0;
  int endIdx = NUM_DF - 1; // intercept is always varying, so this is endIdx
  for (int cov = 0; cov < covariates.n_elem; cov++) {
    if (covIdxs[startIdx] != covIdxs[endIdx]) {
      endIdx = startIdx; // cov not varying, so endIdx is startIdx
      Xtv.cols(startIdx, endIdx) = ones*covariates[cov];
    } else {
      Xtv.cols(startIdx, endIdx) = basis*covariates[cov];
    } 
    startIdx = endIdx + 1;
    endIdx += NUM_DF;
  }

  for (int s = 0; s < NUM_SAMPLE; s++) {
    fit.slice(s) = Xtv*beta.slice(s);
  }

  return fit;
}

//[[Rcpp::export]]
arma::mat getSumExpFit(const arma::cube& fit) {
  const int NUM_ROWS = fit.n_rows,
            NUM_SAMPLE = fit.n_slices;

  cube temp = sum(exp(fit), 1);
  temp.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat sumExpFit = temp.slice(0);

  return sumExpFit;
}

//[[Rcpp::export]]
arma::cube getFitOfData(const arma::cube& beta, const arma::mat& Xtv) {
  /*
   covariates (vec): The first cov must be a 1 for the intercept, the rest will 
   be the covariate profile.
   */
  const int NUM_CAT = beta.n_cols,
            NUM_SAMPLE = beta.n_slices,
            NUM_OBS = Xtv.n_rows;
  cube fit(NUM_OBS, NUM_CAT, NUM_SAMPLE);

  for (int s = 0; s < NUM_SAMPLE; s++) {
    fit.slice(s) = Xtv*beta.slice(s);
  }
  
  return fit;
}

//[[Rcpp::export]]
arma::cube getRAFits(const arma::cube& beta, const arma::mat& basis,
                     const arma::vec& covariates, const arma::vec& covIdxs) {
  /*
   covariates (vec): The first cov must be a 1 for the intercept, the rest will 
   be the covariate profile.
   */
  const int NUM_PARAM = beta.n_rows,
            NUM_CAT = beta.n_cols,
            NUM_SAMPLE = beta.n_slices,
            NUM_TESTPOINTS = basis.n_rows,
            NUM_DF = basis.n_cols;
  cube fit(NUM_TESTPOINTS, NUM_CAT, NUM_SAMPLE);
  mat Xtv(NUM_TESTPOINTS, NUM_PARAM),
      ones(NUM_TESTPOINTS, 1, arma::fill::ones);

  int startIdx = 0;
  int endIdx = NUM_DF - 1;
  for (int cov = 0; cov < covariates.n_elem; cov++) {
    if (covIdxs[startIdx] != covIdxs[endIdx]) {
      endIdx = startIdx; // cov not varying, so endIdx is startIdx
      Xtv.cols(startIdx, endIdx) = ones*covariates[cov];
    } else {
      Xtv.cols(startIdx, endIdx) = basis*covariates[cov];
    } 
    startIdx = endIdx + 1;
    endIdx += NUM_DF;
  }

  for (int s = 0; s < NUM_SAMPLE; s++) {
    fit.slice(s) = exp(Xtv*beta.slice(s));
    fit.slice(s).each_col() /= sum(fit.slice(s), 1);
  }
  
  return fit;
}

//[[Rcpp::export]]
arma::cube getRAFitsPrecalc(const arma::cube& fit, const arma::mat& sumExpFit) {
  const int NUM_CAT = fit.n_cols,
            NUM_SAMPLE = fit.n_slices;

  cube toReturn = arma::exp(fit);

  for (int s = 0; s < NUM_SAMPLE; s++) {
    toReturn.slice(s).each_col() /= sumExpFit.col(s);
  }
  
  return toReturn;
}

//[[Rcpp::export]]
arma::cube getBetaVC(const arma::cube& beta, const arma::uvec VCRows, 
                     const arma::mat& basis) {
  const int NUM_CAT = beta.n_cols,
            NUM_SAMPLE = beta.n_slices,
            NUM_TESTPOINTS = basis.n_rows,
            START_IDX = VCRows[0] - 1, // The -1 is to convert to 0-based index
            END_IDX = VCRows[1] - 1; 
  cube betaVC(NUM_TESTPOINTS, NUM_CAT, NUM_SAMPLE);

  for (int s = 0; s < NUM_SAMPLE; s++) {
    betaVC.slice(s) = basis*beta.slice(s).rows(START_IDX, END_IDX);
  }

  return betaVC;
}  

// ---- Statistics related functions ----

//[[Rcpp::export]]
arma::mat getHillDiversity(const arma::cube& beta, const arma::mat& basis, 
                           const arma::vec& covariates,
                           const arma::vec& covIdxs, const double l = 0) {
  /*
   covariates (vec): The first cov must be a 1 for the intercept, the rest will 
   be the covariate profile.
   */
  const int NUM_PARAM = beta.n_rows,
            NUM_CAT = beta.n_cols,
            NUM_SAMPLE = beta.n_slices,
            NUM_TESTPOINTS = basis.n_rows,
            NUM_DF = basis.n_cols;
  mat fit(NUM_TESTPOINTS, NUM_CAT),
      aDiv(NUM_TESTPOINTS, NUM_SAMPLE),
      Xtv(NUM_TESTPOINTS, NUM_PARAM),
      ones(NUM_TESTPOINTS, 1, arma::fill::ones);
  
  int startIdx = 0;
  int endIdx = NUM_DF - 1;
  for (int cov = 0; cov < covariates.n_elem; cov++) {
    if (covIdxs[startIdx] != covIdxs[endIdx]) {
      endIdx = startIdx; // cov not varying, so endIdx is startIdx
      Xtv.cols(startIdx, endIdx) = ones*covariates[cov];
    } else {
      Xtv.cols(startIdx, endIdx) = basis*covariates[cov];
    } 
    startIdx = endIdx + 1;
    endIdx += NUM_DF;
  }
  
  for (int s = 0; s < NUM_SAMPLE; s++) {
    fit = exp(Xtv*beta.slice(s));
    fit.each_col() /= sum(fit, 1);
    if (l != 0)
      aDiv.col(s) = arma::pow(sum(fit%arma::pow(fit, -1*l), 1), 1/l);
    else
      aDiv.col(s) = sum(fit%log(fit), 1);
  }
  
  if (l != 0)
    return aDiv;
  
  return -1*aDiv;
}

//[[Rcpp::export]]
arma::mat getDeltaHillDiversity(const arma::cube& betaVCFits, 
                                const arma::cube& fit, 
                                const arma::mat& sumExpFit, const double l = 0,
                                double change = 1) {
  const int NUM_ROWS    = betaVCFits.n_rows,
            NUM_SAMPLE  = betaVCFits.n_slices;
  
  // fit of covariate profile you are changing from
  cube fitRAFrom = exp(fit); 
  
  // fit of covariate profile you are changing to
  cube fitRATo = exp(fit + change*betaVCFits); 
  mat sumExpFitTo = sum(exp(fit + change*betaVCFits), 1);
  for (int s = 0; s < NUM_SAMPLE; s++) {
    fitRAFrom.slice(s).each_col() /= sumExpFit.col(s);
    fitRATo.slice(s).each_col() /= sumExpFitTo.col(s);
  }

  cube aDiv = (l != 0) ? (cube) arma::pow(sum(fitRATo%pow(fitRATo, -1*l), 1)
                                /sum(fitRAFrom%pow(fitRAFrom, -1*l), 1), 1/l)
                       : (cube) sum(fitRATo%log(fitRATo), 1)
                                /sum(fitRAFrom%log(fitRAFrom), 1);
  aDiv.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat ret = aDiv.slice(0);
  
  return ret;
} // previously getHillDivMultiChange

//[[Rcpp::export]]
arma::mat getDeltaHillDivMeanAndCI(const arma::cube& betaVCFits,
                                   const arma::cube& fit, 
                                   const arma::mat& sumExpFit, 
                                   const double l = 0, double change = 1) {
  const int NUM_ROWS    = betaVCFits.n_rows,
            NUM_SAMPLE  = betaVCFits.n_slices;
  
  // fit of covariate profile you are changing from
  cube fitRAFrom = exp(fit); 
  
  // fit of covariate profile you are changing to
  cube fitRATo = exp(fit + change*betaVCFits); 
  mat sumExpFitTo = sum(exp(fit + change*betaVCFits), 1);
  for (int s = 0; s < NUM_SAMPLE; s++) {
    fitRAFrom.slice(s).each_col() /= sumExpFit.col(s);
    fitRATo.slice(s).each_col() /= sumExpFitTo.col(s);
  }

  cube aDiv = (l != 0) ? (cube) arma::pow(sum(fitRATo%pow(fitRATo, -1*l), 1)
                                /sum(fitRAFrom%pow(fitRAFrom, -1*l), 1), 1/l)
                       : (cube) sum(fitRATo%log(fitRATo), 1)
                                /sum(fitRAFrom%log(fitRAFrom), 1);
  aDiv.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat ret = aDiv.slice(0);
  
  mat out(ret.n_rows, 3);
  mat quantiles = arma::quantile(ret, arma::vec({0.025, 0.975}), 1);
  
  out.col(0) = quantiles.col(0); // 0.025 quantile
  out.col(1) = arma::mean(ret, 1); // mean
  out.col(2) = quantiles.col(1); // 0.975 quantile
  return out;
}

//[[Rcpp::export]]
arma::cube getAllCatDeltaRA(const arma::cube& betaVCFits, const arma::cube& fit,
                            const arma::mat& sumExpFit, double change = 1) {
  const int NUM_ROWS = betaVCFits.n_rows,
            NUM_CAT = betaVCFits.n_cols,
            NUM_SAMPLE  = betaVCFits.n_slices;
  cube deltaRA(NUM_ROWS, NUM_CAT, NUM_SAMPLE),
       expChange(NUM_ROWS, NUM_CAT, NUM_SAMPLE);

  // adding vB_p(t) to each respective cat, then summing to get denom of deltaRA  
  cube denominatorCube = sum(exp(fit + change*betaVCFits), 1); 
  denominatorCube.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat denominator = denominatorCube.slice(0);
  
  expChange = exp(change*betaVCFits);
  for (int cat = 0; cat < NUM_CAT; cat++) {
    deltaRA.col(cat) = expChange.col_as_mat(cat)%(sumExpFit/denominator);
  }

  return deltaRA;
} 

//[[Rcpp::export]]
arma::mat getDeltaRA(const int cat, const arma::cube& betaVCFits,
                     const arma::cube& fit, const arma::mat& sumExpFit,
                     double change = 1) {
  const int NUM_ROWS    = betaVCFits.n_rows,
            NUM_SAMPLE  = betaVCFits.n_slices;

  // adding vB_p(t) to each respective cat, then summing to get denom of deltaRA  
  cube denominatorCube = sum(exp(fit + change*betaVCFits), 1); 
  denominatorCube.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat denominator = denominatorCube.slice(0);

  mat expChange = exp(change*betaVCFits.col(cat));
  mat deltaRA = expChange%(sumExpFit/denominator);

  return deltaRA;
} // previously getMultChangeVCRAwithCovProfile

//[[Rcpp::export]]
arma::mat getDeltaRAMeanAndCI(const int cat, const arma::cube& betaVCFits, 
                              const arma::cube& fit, const arma::mat& sumExpFit,
                              double change = 1) {
  const int NUM_ROWS    = betaVCFits.n_rows,
            NUM_SAMPLE  = betaVCFits.n_slices;

  // adding vB_p(t) to each respective cat, then summing to get denom of deltaRA  
  cube denominatorCube = sum(exp(fit + change*betaVCFits), 1); 
  denominatorCube.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat denominator = denominatorCube.slice(0);

  mat expChange = exp(change*betaVCFits.col(cat));
  mat deltaRA = expChange%(sumExpFit/denominator);
  
  mat out(deltaRA.n_rows, 3);
  mat quantiles = arma::quantile(deltaRA, arma::vec({0.025, 0.975}), 1);
  
  out.col(0) = quantiles.col(0); // 0.025 quantile
  out.col(1) = arma::mean(deltaRA, 1); // mean
  out.col(2) = quantiles.col(1); // 0.975 quantile
  return out;
}

// [[Rcpp::export]]
double aitchison_distance(const arma::vec& xi, const arma::vec& xj) {
  vec log_ratio_diff = arma::log(xi) - arma::mean(arma::log(xi)) 
                       - (arma::log(xj) - arma::mean(arma::log(xj)));

  return std::sqrt(arma::sum(arma::square(log_ratio_diff)));
} // previously aitchison_distance

// ---- Finding cov for significant statistics related functions ----

//[[Rcpp::export]]
arma::uvec getSignifCovProfileRA(const arma::cube& betaVCFits, 
                                 const arma::mat& betaVCFit,
                                 const arma::cube& fit,
                                 const arma::mat& sumExpFit, 
                                 const double change = 1) {
  const int NUM_ROWS    = betaVCFits.n_rows,
            NUM_SAMPLE  = betaVCFits.n_slices;
  
  // adding vB_p(t) to each respective cat, then summing to get denom of deltaRA  
  cube denominator_cube = sum(exp(fit + change*betaVCFits), 1); 
  denominator_cube.reshape(NUM_ROWS, NUM_SAMPLE, 1);
  mat denominator = denominator_cube.slice(0);
  
  mat expChange = exp(change*betaVCFit);
  mat vcra = expChange%(sumExpFit/denominator);

  arma::vec quant025 = arma::quantile(vcra, arma::vec({0.025}), 1);
  arma::vec quant975 = arma::quantile(vcra, arma::vec({0.975}), 1);
  arma::uvec indices = arma::find((quant025 > 1) || (quant975 < 1));
  
  return indices;
}

//[[Rcpp::export]]
bool getSignifCovProfileAlpha(const arma::cube& fit, const arma::cube& fitBase,
                              const double l = 0) {
  const int NUM_SAMPLE = fit.n_slices,
            NUM_OBS = fit.n_rows;
  mat aDiff(NUM_OBS, NUM_SAMPLE);
  
  cube expFit = exp(fit);
  cube expFitBase = exp(fitBase);
  cube temp = sum(expFit, 1);
  temp.reshape(NUM_OBS, NUM_SAMPLE, 1);
  mat sumExpFit = temp.slice(0);
  temp = sum(expFitBase, 1);
  temp.reshape(NUM_OBS, NUM_SAMPLE, 1);
  mat sumExpFitBase = temp.slice(0);
  for (int s = 0; s < NUM_SAMPLE; s++) {
    mat RA = expFit.slice(s).each_col() / sumExpFit.col(s);
    mat RABase = expFitBase.slice(s).each_col() / sumExpFitBase.col(s);
    if (l != 0)
      aDiff.col(s) =   arma::pow(sum(RA%arma::pow(RA, -1*l), 1), 1/l)
                     - arma::pow(sum(RABase%arma::pow(RABase, -1*l), 1), 1/l);
    else
      aDiff.col(s) = sum(RA%log(RA), 1) - sum(RABase%log(RABase), 1);
  }
  
  arma::vec quant025 = arma::quantile(aDiff, arma::vec({0.025}), 1);
  arma::vec quant975 = arma::quantile(aDiff, arma::vec({0.975}), 1);

  return arma::any((quant025 > 0) || (quant975 < 0));
}

////////////////////////////////////////////////////////////////////////////////
////////////////// Functions for Truth in Simulation Study /////////////////////
////////////////////////////////////////////////////////////////////////////////

arma::vec computeFunction(const arma::vec& t, int funcNum) {
  vec result(t.n_elem, fill::zeros);
  
  if (funcNum == 1) {
    result = (-0.2 * arma::square(t - 5.0) + 5.0) / 7.0;
  } else if (funcNum == 2) {
    result = 1/(1.75 + arma::exp(-1.25*(t - 5.0)));
  } else if (funcNum == 3) {
    result = .07*t;
  } else if (funcNum == 4) {
    result.fill(.5); 
  } else {
    result.fill(0);
  }
  
  return result;
}

//[[Rcpp::export]]
arma::vec getTrueVCRAFromMat(const arma::mat& trueIntercepts, 
                             const arma::vec& testPoints, int catToCheck, 
                             const arma::vec& funcInfo, 
                             bool interceptTV = false) {
  const int NUM_ROWS = testPoints.n_elem,
            NUM_CAT = funcInfo.n_elem;
  double funcNum, plusMinusOne;//scale;
  double xCenter = 0;
  mat X(NUM_ROWS, trueIntercepts.n_rows, fill::ones),
      fit(NUM_ROWS, NUM_CAT, fill::zeros);
  vec vcra(NUM_ROWS, fill::zeros),
      numerator(NUM_ROWS, fill::zeros), 
      denominator(NUM_ROWS, fill::zeros);
  catToCheck = catToCheck - 1; // R to C++ indexing

  fit = trueIntercepts;
  for (int cat = 0; cat < NUM_CAT; cat++) {
    funcNum = std::round(funcInfo[cat]);
    plusMinusOne = (funcInfo[cat] - funcNum >= 0) ? 1 : -1;
    fit.col(cat) += xCenter*plusMinusOne*computeFunction(testPoints, funcNum);
  }
  numerator = arma::sum(arma::exp(fit), 1);
  for (int cat = 0; cat < NUM_CAT; cat++) {
    funcNum = std::round(funcInfo[cat]);
    plusMinusOne = (funcInfo[cat] - funcNum >= 0) ? 1 : -1;
    fit.col(cat) += plusMinusOne*computeFunction(testPoints, funcNum);
  }
  denominator = arma::sum(arma::exp(fit), 1);
  funcNum = std::round(funcInfo[catToCheck]);
  plusMinusOne = (funcInfo[catToCheck] - funcNum >= 0) ? 1 : -1;
  vcra = arma::exp(plusMinusOne*computeFunction(testPoints, funcNum))
         %(numerator/denominator);

  return vcra;
}

//[[Rcpp::export]]
arma::vec getTrueHillDivMultiChange(const arma::mat& baseFit, 
                                    const arma::vec& testPoints, 
                                    const arma::vec& funcInfo, 
                                    const double l = 0, double change = 1) {
  const int NUM_ROWS = testPoints.n_elem,
            NUM_CAT = funcInfo.n_elem;
  double funcNum, plusMinusOne;//scale;
  mat X(NUM_ROWS, baseFit.n_rows, fill::ones),
      fit(NUM_ROWS, NUM_CAT, fill::zeros);
  vec vcra(NUM_ROWS, fill::zeros),
      numerator(NUM_ROWS, fill::zeros), 
      denominator(NUM_ROWS, fill::zeros);
  
  fit = baseFit;
  for (int cat = 0; cat < NUM_CAT; cat++) {
    funcNum = std::round(funcInfo[cat]);
    plusMinusOne = (funcInfo[cat] - funcNum >= 0) ? 1 : -1;
    fit.col(cat) += change*plusMinusOne*computeFunction(testPoints, funcNum);
  }
  
  // fit of covariate profile you are changing from
  mat fitRAFrom = exp(baseFit); 
  fitRAFrom.each_col() /= sum(fitRAFrom, 1);
  // fit of covariate profile you are changing to
  mat fitRATo = exp(fit); 
  vec sumExpFitTo = sum(fitRATo, 1);
  fitRATo.each_col() /= sumExpFitTo;

  mat aDiv = (l != 0) ? 
                    (mat) arma::pow(sum(fitRAFrom%arma::pow(fitRAFrom, -1*l), 1)
                                 /sum(fitRATo%arma::pow(fitRATo, -1*l), 1), 1/l)
                       : (mat) sum(fitRAFrom%log(fitRAFrom), 1)
                         /sum(fitRATo%log(fitRATo), 1);
  aDiv.reshape(NUM_ROWS, 1);
  vec ret = aDiv.col(0);
  
  return ret;
}

#endif
