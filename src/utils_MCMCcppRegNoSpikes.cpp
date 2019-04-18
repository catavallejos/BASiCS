#include "utils.h"

/* Metropolis-Hastings updates of mu
 * Updates are implemented simulateaneously for all biological genes
 */
arma::mat muUpdateRegNoSpikes(
    arma::vec const& mu0,
    arma::vec const& prop_var,
    arma::mat const& Counts,
    arma::vec const& delta,
    arma::vec const& invdelta,
    arma::vec const& nu,
    arma::vec const& sum_bycell_all,
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u,
    arma::vec & ind,
    double const& Constrain, /* No-spikes arguments from here */
    int const& RefGene,
    arma::uvec const& ConstrainGene,
    arma::uvec const& NotConstrainGene,
    int const& ConstrainType,
    int const& k, /* Regression arguments from here */
    arma::vec const& lambda,
    arma::vec const& beta,
    arma::mat const& X,
    double const& sigma2,
    double variance,
    double exponent)
{
  using arma::span;

  int nConstrainGene = ConstrainGene.size();
  int nNotConstrainGene = NotConstrainGene.size();

  // PROPOSAL STEP
  mu1 = exp(arma::randn(q0) % sqrt(prop_var) + log(mu0));
  u = arma::randu(q0);

  // INITIALIZE MU
  double aux;
  int iAux;
  double sumAux = sum(log(mu0.elem(ConstrainGene))) - log(mu0(RefGene));

  // ACCEPT/REJECT STEP

  // Step 1: Computing the likelihood contribution of the acceptance rate
  // Calculated in the same way for all genes,
  // but the reference one (no need to be sequential)
  arma::vec log_aux = (log(mu1) - log(mu0)) % sum_bycell_all;
  for (int i = 0; i < q0; i++) {
    if (i != RefGene) {
      for (int j = 0; j < n; j++) {
        log_aux(i) -= ( Counts(i, j) + invdelta(i) ) *
          log(
            (nu(j) * mu1(i) + invdelta(i)) /
            (nu(j) * mu0(i) + invdelta(i)));
      }
    }
  }
  // Revise this part
  // This is new due to regression prior on delta
  arma::mat X_mu1 = designMatrix(k, mu1, variance);

  // REGRESSION RELATED FACTOR
  log_aux -= exponent * lambda %
    (pow(log(delta) - X_mu1 * beta, 2) -
      pow(log(delta) - X * beta, 2)) /
    (2 * sigma2);

  // Step 2: Computing prior component of the acceptance rate

  // Step 2.1: For genes that are under the constrain (excluding the reference one)
  for (int i = 0; i < nConstrainGene; i++) {
    iAux = ConstrainGene(i);
    if (iAux != RefGene) {
      aux = 0.5 * (ConstrainGene.size() * Constrain - (sumAux - log(mu0(iAux))));
      log_aux(iAux) -= (0.5 * 2 /s2_mu) * (pow(log(mu1(iAux)) - aux, 2));
      log_aux(iAux) += (0.5 * 2 /s2_mu) * (pow(log(mu0(iAux)) - aux, 2));
      // ACCEPT REJECT
      if ((log(u(iAux)) < log_aux(iAux)) & (mu1(iAux) > 1e-3)) {
        ind(iAux) = 1;
        sumAux += log(mu1(iAux)) - log(mu0(iAux));
      } else{
        ind(iAux) = 0;
        mu1(iAux) = mu0(iAux);
      }
    }
  }

  // Step 2.2: For the reference gene
  ind(RefGene) = 1;
  mu1(RefGene) = exp(ConstrainGene.size() * Constrain - sumAux);

  // Step 2.3: For genes that are *not* under the constrain
  // Only relevant for a trimmed constrain
  if (ConstrainType == 2) {
    for (int i=0; i < nNotConstrainGene; i++) {
      iAux = NotConstrainGene(i);
      log_aux(iAux) -= (0.5 / s2_mu) * (pow(log(mu1(iAux)), 2) - pow(log(mu0(iAux)), 2));
      // ACCEPT REJECT
      if ((log(u(iAux)) < log_aux(iAux)) & (mu1(iAux) > 1e-3)) {
        ind(iAux) = 1;
      } else {
        ind(iAux) = 0;
        mu1(iAux) = mu0(iAux);
      }
    }
  }
  // OUTPUT
  return join_rows(mu1, ind);
}


/* Metropolis-Hastings updates of delta
* Updates are implemented simulateaneously for all biological genes
*/
arma::mat deltaUpdateRegNoSpikes(
    arma::vec const& delta0,
    arma::vec const& prop_var,
    arma::mat const& Counts,
    arma::vec const& mu,
    arma::vec const& nu,
    int const& q0,
    int const& n,
    arma::vec & delta1,
    arma::vec & u,
    arma::vec & ind,
    arma::vec const& lambda,
    arma::mat const& X,
    double const& sigma2,
    arma::vec const& beta,
    double exponent)
{
  using arma::span;

  // PROPOSAL STEP
  delta1 = exp(arma::randn(q0) % sqrt(prop_var) + log(delta0));
  u = arma::randu(q0);

  // ACCEPT/REJECT STEP
  // +1 should appear because we update log(delta) not delta.
  // However, it cancels out with the prior.
  arma::vec log_aux = - n * (lgamma_cpp(1 / delta1) - lgamma_cpp(1 / delta0));
  log_aux -= n * ((log(delta1) / delta1) - (log(delta0) / delta0));
  for (int i = 0; i < q0; i++) {
    for (int j = 0; j < n; j++) {
      log_aux(i) += std::lgamma(Counts(i, j) + (1 / delta1(i)));
      log_aux(i) -= std::lgamma(Counts(i, j) + (1 / delta0(i)));
      log_aux(i) -= (Counts(i, j) +
        (1 / delta1(i))) *
        log(nu(j) * mu(i) + (1 / delta1(i)));
      log_aux(i) += (Counts(i, j) +
        (1 / delta0(i))) *
        log(nu(j) * mu(i) + (1 / delta0(i)));
    }
  }


  // from MCMC_cppReg.cpp
  // The next operations are equivalent; second one a is simplified version
  //  log_aux -= lambda % (pow(log(delta1) - X * beta, 2) -
  //    pow(log(delta0) - X * beta, 2)) / (2 * sigma2);
  log_aux -= exponent * lambda %
    (pow(log(delta1), 2) - pow(log(delta0), 2) -
    2 * (log(delta1) - log(delta0)) % (X * beta)) / (2 * sigma2);

  // CREATING OUTPUT VARIABLE & DEBUG
  ind = DegubInd(ind, q0, u, log_aux, delta1, 1e-3, "delta");
  for (int i = 0; i < q0; i++) {
    if (ind(i) == 0) {
      delta1(i) = delta0(i);
    }
  }

  // OUTPUT
  return join_rows(delta1, ind);
}
