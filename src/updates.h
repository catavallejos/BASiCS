#ifndef UPDATES_H
#define UPDATES_H

#include "utils.h"

/* Metropolis-Hastings updates of mu 
* Updates are implemented simulateaneously for all biological genes
*/
// [[Rcpp::export(".muUpdate")]]
arma::mat muUpdate(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts, 
    arma::vec const& invdelta, 
    arma::vec const& phinu, 
    arma::vec const& sum_bycell_bio,
    arma::vec const& mu_mu,
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u, 
    arma::vec & ind,
    double const& exponent,
    double const& mintol)
{
  
  /* PROPOSAL STEP */
  mu1 = exp( arma::randn(q0) % sqrt(prop_var) + log(mu0) );
  u = arma::randu(q0);
  
  /* ACCEPT/REJECT STEP 
  * Note: there is a -1 factor coming from the log-normal prior. 
  * However, it cancels out as using log-normal proposals.
  */
  arma::vec log_aux = (log(mu1) - log(mu0)) % sum_bycell_bio; 
  log_aux -= (0.5 / s2_mu) * 
    (pow(log(mu1) - mu_mu, 2) - pow(log(mu0) - mu_mu, 2))
    * exponent;
  #pragma omp parallel for
  for (int i=0; i < q0; i++) {
    for (int j=0; j < n; j++) {
      log_aux(i) -= ( Counts(i,j) + invdelta(i) ) *  
        log( ( phinu(j)*mu1(i) + invdelta(i) ) / 
        ( phinu(j)*mu0(i) + invdelta(i) ));
    }
  }
  
  /* CREATING OUTPUT & DEBUG 
  * Proposed values are automatically rejected in the following cases:
  * - If smaller than 1e-3
  * - If the proposed value is not finite
  * - When the acceptance rate cannot be numerally computed
  */
  ind = DegubInd(ind, q0, u, log_aux, mu1, mintol, "mu");
  for (int i=0; i < q0; i++) {
    if(ind(i) == 0) mu1(i) = mu0(i);  
  }
  
  /* OUTPUT */
  return join_rows(mu1, ind);
}

/* Metropolis-Hastings updates of delta
* Updates are implemented simulateaneously for all biological genes.
*/
// [[Rcpp::export(".deltaUpdate")]]
arma::mat deltaUpdate(
    arma::vec const& delta0, 
    arma::vec const& prop_var,  
    arma::mat const& Counts, 
    arma::vec const& mu, 
    arma::vec const& phinu, 
    double const& a_delta, 
    double const& b_delta, 
    double const& s2delta,
    double const& prior_delta,
    int const& q0,
    int const& n,
    arma::vec & delta1,
    arma::vec & u, 
    arma::vec & ind,
    double const& exponent,
    double const& mintol)
{
  
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
      log_aux(i) -= (Counts(i, j) + (1 / delta1(i))) * 
        log(phinu(j) * mu(i) + (1 / delta1(i)));
      log_aux(i) += (Counts(i, j) + (1 / delta0(i))) * 
        log(phinu(j) * mu(i) + (1 / delta0(i)));
    }
  }
  // Component related to the prior
  if(prior_delta == 1) {
    log_aux += (log(delta1) - log(delta0)) * 
      a_delta - b_delta * 
      (delta1 - delta0);
    log_aux *= exponent;
  } else { 
    log_aux -= (0.5 / s2delta) *
      (pow(log(delta1), 2) - pow(log(delta0), 2)) * 
      exponent;
  }
  
  /* CREATING OUTPUT VARIABLE & DEBUG 
  * Proposed values are automatically rejected in the following cases:
  * - If smaller than 1e-3
  * - If the proposed value is not finite
  * - When the acceptance rate cannot be numerally computed
  */ 
  ind = DegubInd(ind, q0, u, log_aux, delta1, mintol, "delta");
  for (int i=0; i < q0; i++) {
    if(ind(i) == 0) {
      delta1(i) = delta0(i);
    }
  }
  
  // OUTPUT
  return join_rows(delta1, ind);
}


/* Metropolis-Hastings updates of phi 
* Joint updates using Dirichlet proposals
*/
// [[Rcpp::export(".phiUpdate")]]
Rcpp::List phiUpdate(
    arma::vec const& phi0,
    double const& prop_var,
    arma::mat const& Counts,
    arma::vec const& mu,
    arma::vec const& invdelta,
    arma::vec const& nu,
    arma::vec const& aphi,
    arma::vec const& sum_bygene_bio,
    int const& q0,
    int const& n,
    arma::vec & phi1,
    double const& exponent) {

  int ind;
  
  // PROPOSAL STEP
  phi1 = n * rDirichlet(prop_var * phi0); 
  double u = R::runif(0,1);
  
  // ACCEPT/REJECT STEP (REJECT VALUES OUTSIDE VALID RANGE)  
  if(all(prop_var * phi1 < 2.5327372760800758e+305)  && 
     all(prop_var * phi0 < 2.5327372760800758e+305) &&
     all(phi1 > 0) && all(phi0 > 0)) {
    // There is an extra -1 but it cancels out with the proposal component
    double log_aux = sum(
      (sum_bygene_bio + (aphi * exponent)) % (log(phi1) - log(phi0))
    );
    
    // Loop to replace matrix operations, through genes and cells
    // There is an extra factor in the prior n^(-n); it cancels out in the ratio
    // There is an extra factor n^(-(sum(aphi) - 1));it cancels out in the ratio
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < q0; i++) {
        log_aux -= (Counts(i, j) + invdelta(i)) *
          log(
            (phi1(j) * nu(j) * mu(i) + invdelta(i)) / 
            (phi0(j) * nu(j) * mu(i) + invdelta(i))
          );
      }
    }

    // There is an extra factor 
    // n^(-(sum(prop_var * phi1))) / n^(-(sum(prop_var * phi0)));
    // it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
    // There is an extra factor 
    // gamma(sum(prop_var * phi0)) / gamma(sum(prop_var * phi1));
    // it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
    log_aux += prop_var * sum(phi1 % log(phi0) - phi0 % log(phi1));
    log_aux -= sum(
      lgamma_cpp_vec(prop_var * phi1) -
      lgamma_cpp_vec(prop_var * phi0)
    );    
    
    if(!R_IsNA(log_aux)){
      if (log(u) < log_aux) {
        ind = 1;
      } else {
        ind = 0;
        phi1 = phi0;
      }
    }
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
    else {
      Rcpp::Rcout << "Error when updating phi" << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset."); 
      ind = 0;
      phi1 = phi0;
    }
  } else {
    ind = 0;
    phi1 = phi0;
  }      
  return(
    Rcpp::List::create(
      Rcpp::Named("phi") = phi1,
      Rcpp::Named("ind") = ind
    )
  ); 
}

/* Draws for cell-specific normalising constants s[j] (batch case)
 * Metropolis-Hastings updates are not required as full conditionals 
 * have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells.
 */
// [[Rcpp::export(".sUpdateBatch")]]
arma::vec sUpdateBatch(
    arma::vec const& s0, 
    arma::vec const& nu, 
    arma::vec const& thetaBatch, 
    double const& as, 
    double const& bs, 
    arma::mat const& BatchDesign,
    int const& n,
    arma::vec & s1,
    double const& exponent) {
  
  // Calculating parameters to the passed as input to the Rgig function (common for all cells)
  arma::vec lambda;
  double psi = 2 * bs;
  if (exponent == 1) {
    lambda = as - 1 / thetaBatch;
  } else {
    lambda = ((as - 1) * exponent) - (1 / thetaBatch) + 1;
    psi *= exponent;
  }

  // GIG draws
  // Initialize s1 with s0 (return that for invalid samples)
  
  // Calculating parameter to the passed as input to the Rgig function (specific to each cell)
  arma::vec chi = 2 * nu / thetaBatch;
  for (int j = 0; j < n; j++) {
    if (!R_IsNA(lambda(j))) {
      if (!R_IsNA(chi(j)) && (chi(j) > 0)) {
        s1(j) = Rcpp::as<double>(Rgig(1, lambda(j), chi(j), psi));
        /* DEBUG: break in case of undefined values */
        if (R_IsNA(s1(j))) {
          Rcpp::Rcout << "Error when updating s" << j << std::endl;
          Rcpp::stop("Please consider additional filter of the input dataset.");
        }
      } else {
        if (!(chi(j) < 0) && (lambda(j) > 0)) {
          s1(j) = Rcpp::as<double>(Rgig(1, lambda(j), chi(j), psi));
        }
      }
    }
  }
  return s1;     
}

/* Metropolis-Hastings updates of nu (batch case)
* Updates are implemented simulateaneously for all cells.
*/
// [[Rcpp::export(".nuUpdateBatch")]]
arma::mat nuUpdateBatch(
    arma::vec const& nu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,
    double const& SumSpikeInput,
    arma::mat const& BatchDesign, 
    arma::vec const& mu, 
    arma::vec const& invdelta, 
    arma::vec const& phi, 
    arma::vec const& s,
    arma::vec const& thetaBatch, 
    arma::vec const& sum_bygene_all, 
    int const& q0,
    int const& n,
    arma::vec & nu1,
    arma::vec & u,
    arma::vec & ind,
    double const& exponent,
    double const& mintol) {

  using arma::span;
  
  // PROPOSAL STEP    
  nu1 = exp(arma::randn(n) % sqrt(prop_var) + log(nu0));
  u = arma::randu(n);
  
  // ACCEPT/REJECT STEP
  arma::vec log_aux = arma::zeros(n);
  #pragma omp parallel for
  for (int j=0; j < n; j++) {
    for (int i=0; i < q0; i++) {
      log_aux(j) -= (Counts(i, j) + invdelta(i)) *
        log(
          (phi(j) * nu1(j) * mu(i) + invdelta(i)) / 
          (phi(j) * nu0(j) * mu(i) + invdelta(i))
        );
    } 
  }
  
  log_aux += (log(nu1) - log(nu0)) % (sum_bygene_all + 1 / thetaBatch)
    * exponent;
  log_aux -= (nu1 - nu0)  % (SumSpikeInput + (exponent / (thetaBatch % s)));
  
  /* CREATING OUTPUT VARIABLE & DEBUG 
  * Proposed values are automatically rejected in the following cases:
  * - If smaller than 1e-5
  * - If the proposed value is not finite
  * - When the acceptance rate cannot be numerally computed
  */ 
  ind = DegubInd(ind, n, u, log_aux, nu1, mintol, "nu");
  for (int j=0; j < n; j++) {
    if(ind(j) == 0) {
      nu1(j) = nu0(j);
    }
  }
  
  // OUTPUT
  return join_rows(nu1, ind);
}

/* Metropolis-Hastings updates of theta 
*/
/* theta0: Current value of $\theta$ */
/* prop_var: Current value of the proposal variances for $\theta$ */
/* s: Current value of $s=(s_1,...,s_n)$' */
/* nu: Current value of $\nu=(\nu_1,...,\nu_n)'$ */
/* a_theta: Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
/* b_theta: Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
// [[Rcpp::export(".thetaUpdateBatch")]]
arma::mat thetaUpdateBatch(
    arma::vec const& theta0,
    arma::vec const& prop_var,
    arma::mat const& BatchDesign,
    arma::vec const& BatchSizes,
    arma::vec const& s,
    arma::vec const& nu,
    double const& a_theta,
    double const& b_theta,
    int const& n,
    int const& nBatch,
    double const& exponent,
    double const& mintol) {

  using arma::span;
  
  // CREATING VARIABLES WHERE TO STORE DRAWS
  arma::vec logtheta = log(theta0);
  
  // PROPOSAL STEP
  arma::vec y = arma::randn(nBatch) % sqrt(prop_var) + logtheta;
  arma::vec u = arma::randu(nBatch);
  
  arma::mat BatchDesignAux = BatchDesign.cols(0, nBatch - 1);
  BatchDesignAux.each_col() %= log(nu / s) - (nu / s);
  
  // ACCEPT/REJECT STEP
  arma::vec log_aux = (y - logtheta);
  if (exponent == 1) {
    log_aux *= a_theta;
    log_aux -= BatchSizes % (logtheta / theta0) % 
      ((y / logtheta) % exp(-y + logtheta) - 1);
  } else {
    log_aux *= ((a_theta - 1) * exponent);
    log_aux -= BatchSizes % (logtheta / theta0) %
      ((y / logtheta) % exp(-y + logtheta));
  }
  // Gamma component
  log_aux -= BatchSizes % (lgamma_cpp(exp(-y)) - lgamma_cpp(1 / theta0));
  // nu / s component
  log_aux += ((exp(-y + logtheta) - 1) / theta0) % sum(BatchDesignAux, 0).t();
  // exponential component
  log_aux -= (b_theta * exponent) * theta0 % (exp(y - logtheta) - 1);
  // log_aux -= (b_theta) * theta0 % (exp(y - logtheta) - 1);
  arma::umat ind = log(u) < log_aux;
  // DEBUG: Reject proposed values below 0.0001 (to avoid numerical innacuracies)
  ind %= mintol < exp(y);
  
  // CREATING OUTPUT VARIABLE
  arma::vec theta = ind % exp(y) + (1 - ind) % theta0;
  
  // OUTPUT
  return join_rows(theta, arma::conv_to<arma::mat>::from(ind));
}

#endif
