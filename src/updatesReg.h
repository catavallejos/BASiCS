#ifndef UPDATESREG_H
#define UPDATESREG_H

#include "utils.h"

// Functions for regression case of BASiCS

// Model matrix generation for regression
// [[Rcpp::export(".designMatrix")]]
arma::mat designMatrix(
    int const& k, /* Number of Gaussian radial basis functions to use for regression */
    arma::vec RBFLocations,
    arma::vec const& mu,
    double const& variance) {

  
  arma::vec x = log(mu);
  // arma::vec RBFLocations = estimateRBFLocations(x, k);
  double h = (RBFLocations(1) - RBFLocations(0)) * variance;

  // Possibly create this matrix outside
  arma::mat X = arma::ones(x.size(), k);
  X.col(1) = x;
  for (int i = 0; i < k - 2; i++) {
    X.col(i + 2) = exp(-0.5 * pow(x - RBFLocations(i), 2) / pow(h, 2));
  }
  return X;
}

// [[Rcpp::export(".estimateRBFLocations")]]
arma::vec estimateRBFLocations(
    arma::vec const& log_mu,
    int const& k,
    bool RBFMinMax) {

  double ran = log_mu.max() - log_mu.min();
  arma::vec RBFLocations = arma::vec(k - 2);
  double d = ran / (k - 1);
  if (RBFMinMax) {
    RBFLocations(0) = log_mu.min();
    for (unsigned int i = 1; i < RBFLocations.size(); i++) {
      // RBFLocations(i) = RBFLocations(i - 1) + d;
      RBFLocations(i) =  RBFLocations(i - 1) + ran / (k - 3);
    }  
  } else {
    RBFLocations(0) = log_mu.min() + d;
    for (unsigned int i = 1; i < RBFLocations.size(); i++) {
      RBFLocations(i) = RBFLocations(i - 1) + d;
      //RBFLocations(i) =  RBFLocations(i - 1) + ran / (k - 3);
    }    
  }
  return RBFLocations;
}


/* Metropolis-Hastings updates of mu 
* Updates are implemented simulateaneously for all biological genes
*/
// [[Rcpp::export(".muUpdateReg")]]
arma::mat muUpdateReg(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts, 
    arma::vec const& delta, 
    arma::vec const& phinu, 
    arma::vec const& sum_bycell_bio,
    arma::vec const& mu_mu,
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u, 
    arma::vec & ind,
    int const& k,
    arma::vec const& lambda,
    arma::vec const& beta,
    arma::mat const& X,
    double const& sigma2,
    double variance,
    bool FixLocations,
    bool RBFMinMax,
    arma::vec RBFLocations,
    double const& exponent,
    double const& mintol) {
  
  /* PROPOSAL STEP */
  mu1 = exp( arma::randn(q0) % sqrt(prop_var) + log(mu0) );
  u = arma::randu(q0);
  
  /* ACCEPT/REJECT STEP 
  * Note: there is a -1 factor coming from the log-normal prior. 
  * However, it cancels out as using log-normal proposals.
  */
  arma::vec log_aux = (log(mu1) - log(mu0)) % sum_bycell_bio; 
  log_aux -= (0.5 / s2_mu) *
    (pow(log(mu1) - mu_mu, 2) - pow(log(mu0) - mu_mu, 2)) * exponent;
  #pragma omp parallel for
  for (int i = 0; i < q0; i++) {
    for (int j = 0; j < n; j++) {
      log_aux(i) -= ( Counts(i,j) + 1/delta(i) ) *  
        log(
          (phinu(j) * mu1(i) + 1 / delta(i)) / 
          (phinu(j) * mu0(i) + 1 / delta(i))
        );
    }
  }
  
  // This is new due to regression prior on delta
  if (!FixLocations) {  
    RBFLocations = estimateRBFLocations(log(mu1), k, RBFMinMax);
  }
  arma::mat X_mu1 = designMatrix(k, RBFLocations, mu1, variance);
  
  // REGRESSION RELATED FACTOR
  // Some terms might cancel out here; check
  log_aux -= exponent * lambda % 
    (
      pow(log(delta) - X_mu1 * beta, 2) - 
      pow(log(delta) - X * beta, 2)
    ) / (2 * sigma2);
  
  /* CREATING OUTPUT VARIABLE & DEBUG 
  * Proposed values are automatically rejected in the following cases:
  * - If smaller than 1e-3
  * - If the proposed value is not finite
  * - When the acceptance rate cannot be numerally computed
  */
  ind = DegubInd(ind, q0, u, log_aux, mu1, mintol, "mu");
  for (int i = 0; i < q0; i++) {
    if (ind(i) == 0) {
      mu1(i) = mu0(i);
    }
  }
  
  /* OUTPUT */
  return join_rows(mu1, ind);
}

/* Metropolis-Hastings updates of delta
* Updates are implemented simulateaneously for all biological genes.
* delta: current value of delta
* prop_var: current value of the proposal variances for delta
* Counts: matrix of expression counts
* mu: current value of mu
* phi: current value of phi
* nu: current value of nu
* a_delta: shape prior hyper-parameter for delta (when using a gamma prior)
* b_delta: rate prior hyper-parameter for delta (when using a gamma prior)
* s2delta: prior hyper-variance for delta (when using a log-normal prior)
* q0: number of biological genes
* n: number of cells
* prior_delta: whether gamma or log-normal priors are being used
* delta: auxiliary vector for storage
* ind: auxiliary vector for storage
* lambda: current value of lambda
* X: current design matrix
* sigma2: current value of sigma2
* beta: current value of beta
*/
// [[Rcpp::export(".deltaUpdateReg")]]
arma::mat deltaUpdateReg(
    arma::vec const& delta0, 
    arma::vec const& prop_var,  
    arma::mat const& Counts, 
    arma::vec const& mu, 
    arma::vec const& phinu, 
    int const& q0,
    int const& n,
    arma::vec & delta1,
    arma::vec & u, 
    arma::vec & ind,
    arma::vec const& lambda,
    arma::mat const& X,
    double const& sigma2,
    arma::vec const& beta,
    double const& exponent,
    double const& mintol) {
  
  /* PROPOSAL STEP */
  delta1 = exp(arma::randn(q0) % sqrt(prop_var) + log(delta0));
  u = arma::randu(q0);
  
  /* ACCEPT/REJECT STEP 
  * Note: there is a -1 factor coming from the log-normal prior. 
  * However, it cancels out as using log-normal proposals.
  */
  arma::vec log_aux = - n * (lgamma_cpp(1/delta1) - lgamma_cpp(1/delta0));
  log_aux -= n * ( (log(delta1)/delta1) - (log(delta0)/delta0) );
  #pragma omp parallel for
  for (int i=0; i < q0; i++) {
    for (int j=0; j < n; j++) {
      log_aux(i) += std::lgamma(Counts(i,j) + (1/delta1(i)));
      log_aux(i) -= std::lgamma(Counts(i,j) + (1/delta0(i)));
      log_aux(i) -= ( Counts(i,j)+(1/delta1(i)) ) * log( phinu(j)*mu(i)+(1/delta1(i)) );
      log_aux(i) += ( Counts(i,j)+(1/delta0(i)) ) * log( phinu(j)*mu(i)+(1/delta0(i)) );
    }
  }
  
  // REGRESSION RELATED FACTOR
  // The next lines are equivalent; second one a is simplified version
  //  log_aux -= lambda%(pow(log(delta1)-X*beta,2) - pow(log(delta0)-X*beta,2))/(2*sigma2);
  log_aux -= exponent * lambda % 
    (
      pow(log(delta1), 2) - pow(log(delta0), 2) -
      2 * 
      (log(delta1) - log(delta0)) % 
      (X * beta)
    ) / (2 * sigma2);
    
  /* CREATING OUTPUT VARIABLE & DEBUG 
  * Proposed values are automatically rejected in the following cases:
  * - If smaller than 1e-3
  * - If the proposed value is not finite
  * - When the acceptance rate cannot be numerally computed
  */    
  ind = DegubInd(ind, q0, u, log_aux, delta1, mintol, "delta");
  for (int i = 0; i < q0; i++) {
    if (ind(i) == 0) {
      delta1(i) = delta0(i);
    }
  }
  
  // OUTPUT
  return join_rows(delta1, ind);
}
// [[Rcpp::export(".betaUpdateReg")]]
arma::vec betaUpdateReg(double const& sigma2,
                        arma::mat const& VAux,
                        arma::vec const& mAux) {

  arma::mat MVRNORM = mvrnormArma(1,mAux,sigma2 * VAux);
  arma::vec beta = MVRNORM.row(0).t();
  return beta;
}
// [[Rcpp::export(".sigma2UpdateReg")]]
double sigma2UpdateReg(arma::vec const& delta,
                       arma::vec const& beta,
                       arma::vec const& lambda, 
                       arma::mat const& V1,
                       double const& mInvVm0,
                       arma::vec const& m,
                       double const& sigma2_a0,
                       double const& sigma2_b0,
                       int const& q0,
                       double const& exponent) {
  double a, b;
  a = sigma2_a0 + (q0 + (beta.n_elem * exponent)) / 2;
  b = sigma2_b0 + 0.5 * mInvVm0;
  b += 0.5 * Rcpp::as<double>(wrap(beta.t() * V1 * beta - 2 * beta.t() * V1 * m));
  b += 0.5 * sum(lambda % pow(log(delta), 2));
  
  double sigma2 = pow(R::rgamma(a, 1.0 / b), -1);
  return sigma2; 
}
// [[Rcpp::export(".lambdaUpdateReg")]]
arma::vec lambdaUpdateReg(arma::vec const& delta,
                          arma::mat const& X,
                          arma::vec const& beta,
                          double const& sigma2, 
                          double const& eta, 
                          int const& q0,
                          arma::vec lambda1,
                          double const& exponent) {
  double a;
  arma::vec b;
  if (exponent == 1) {
    a = (eta + 1) / 2;
    b = 0.5 * (eta + (pow(log(delta) - X * beta, 2) / sigma2) );
  } else {
    a = ((eta / 2) - 1) * exponent + 1.5;
    b = 0.5 * (
      (eta * exponent) + (pow(log(delta) - (X * beta), 2) / sigma2)
    );
  }

  // Parameter calculations
  for (int i = 0; i < q0; i++) {
    lambda1(i) = R::rgamma(a, 1.0 / b(i));
  }
  return lambda1;
}

#endif
