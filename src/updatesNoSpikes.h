#ifndef UPDATESNOSPIKES_H
#define UPDATESNOSPIKES_H

#include "utils.h"

/* Metropolis-Hastings updates of mu 
* Updates are implemented simulateaneously for all biological genes
*/
// [[Rcpp::export(".muUpdateNoSpikes")]]
arma::mat muUpdateNoSpikes(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,  
    arma::vec const& invdelta, 
    arma::vec const& nu, 
    arma::vec const& sum_bycell_all,
    arma::vec const& mu_mu,
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u,
    arma::vec & ind,
    double const& SizeTimesConstrain,
    int const& RefGene,
    arma::uvec const& ConstrainGene,
    arma::uvec const& NotConstrainGene,
    double const& exponent,
    double const& mintol) {
  using arma::span;
  
  int nConstrainGene = ConstrainGene.size();
  int nNotConstrainGene = NotConstrainGene.size();
  
  // PROPOSAL STEP    
  mu1 = exp(arma::randn(q0) % sqrt(prop_var) + log(mu0));
  u = arma::randu(q0);
  
  // INITIALIZE MU
  double aux; double iAux;
  double sumAux = sum(log(mu0.elem(ConstrainGene))) - log(mu0(RefGene));
  
  // ACCEPT/REJECT STEP
  
  // Step 1: Computing the likelihood contribution of the acceptance rate 
  // Calculated in the same way for all genes, 
  // but the reference one (no need to be sequential)
  arma::vec log_aux = (log(mu1) - log(mu0)) % sum_bycell_all;
  #pragma omp parallel for
  for (int i = 0; i < q0; i++) {
    if (i != RefGene) {
      for (int j = 0; j < n; j++) {
        log_aux(i) -= (Counts(i, j) + invdelta(i)) * 
          log(
            (nu(j)*mu1(i) + invdelta(i)) / 
            (nu(j)*mu0(i) + invdelta(i))
          );
      }
    }
  }
  
  // Step 2: Computing prior component of the acceptance rate 
  
  // Step 2.1: For genes that are under the constrain (excluding the reference one)
  for (int i=0; i < nConstrainGene; i++) {
    iAux = ConstrainGene(i);
    if(iAux != RefGene) {
      aux = 0.5 * (SizeTimesConstrain - (sumAux - log(mu0(iAux))));
      aux += 0.5 * (mu_mu(iAux) - mu_mu(RefGene));
      log_aux(iAux) -= (0.5 * 2 /s2_mu) * 
        (pow(log(mu1(iAux)) - aux,2)) * exponent; 
      log_aux(iAux) += (0.5 * 2 /s2_mu) * 
        (pow(log(mu0(iAux)) - aux,2)) * exponent;
      
      // ACCEPT REJECT
      if ((log(u(iAux)) < log_aux(iAux)) && (mu1(iAux) > mintol)) {
        ind(iAux) = 1;
        sumAux += log(mu1(iAux)) - log(mu0(iAux)); 
      } else {
        ind(iAux) = 0;
        mu1(iAux) = mu0(iAux);
      }
    }
  }
  
  // Step 2.2: For the reference gene 
  ind(RefGene) = 1;
  mu1(RefGene) = exp(SizeTimesConstrain - sumAux);
  
  // Step 2.3: For genes that are *not* under the constrain
  // Only relevant for a trimmed constrain
  #pragma omp parallel for
  for (int i=0; i < nNotConstrainGene; i++) {
    iAux = NotConstrainGene(i);
    log_aux(iAux) -= (0.5 / s2_mu) * 
      (
          pow(log(mu1(iAux)) - mu_mu(iAux), 2) - 
          pow(log(mu0(iAux)) - mu_mu(iAux), 2)
      ) * exponent;

    // ACCEPT REJECT
    if ((log(u(iAux)) < log_aux(iAux)) && (mu1(iAux) > mintol)) {
      ind(iAux) = 1;
    } else {
      ind(iAux) = 0;
      mu1(iAux) = mu0(iAux);
    }
  }
  
  // OUTPUT
  return join_rows(mu1, ind);
}

/* Metropolis-Hastings updates of nu (batch case)
* Updates are implemented simulateaneously for all cells.
*/
// [[Rcpp::export(".nuUpdateBatchNoSpikes")]]
arma::mat nuUpdateBatchNoSpikes(
    arma::vec const& nu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,
    arma::mat const& BatchDesign, 
    arma::vec const& mu, 
    arma::vec const& invdelta, 
    arma::vec const& s, 
    arma::vec const& thetaBatch, 
    arma::vec const& sum_bygene_all, 
    int const& q0,
    int const& n,
    arma::vec & nu1,
    arma::vec & u,
    arma::vec & ind,
    double const& exponent,
    double const& mintol)
{
  using arma::span;
  
  // PROPOSAL STEP    
  nu1 = exp(arma::randn(n) % sqrt(prop_var) + log(nu0));
  u = arma::randu(n);
  
  // ACCEPT/REJECT STEP
  arma::vec log_aux = (log(nu1) - log(nu0)) % 
    (sum_bygene_all + (1 / thetaBatch * exponent));
  log_aux -= (nu1 - nu0) % (exponent / (thetaBatch % s));
  
  #pragma omp parallel for
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < q0; i++) {
      log_aux(j) -= (Counts(i,j) + invdelta(i)) *
        log(
          (nu1(j) * mu(i) + invdelta(i)) / 
          (nu0(j) * mu(i) + invdelta(i)));
    } 
  }
  
  /* CREATING OUTPUT VARIABLE & DEBUG 
  * Proposed values are automatically rejected in the following cases:
  * - If smaller than 1e-5
  * - If the proposed value is not finite
  * - When the acceptance rate cannot be numerally computed
  */  
  ind = DegubInd(ind, n, u, log_aux, nu1, mintol, "nu");
  for (int j = 0; j < n; j++) {
    if (ind(j) == 0) {
      nu1(j) = nu0(j);
    }
  }
  
  // OUTPUT
  return join_rows(nu1, ind);
}

#endif
