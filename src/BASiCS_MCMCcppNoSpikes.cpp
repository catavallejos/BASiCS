#include "utils.h"

/* MCMC sampler for the non-spike case
 * N: Total number of MCMC draws
 * Thin: Thinning period for MCMC chain
 * Burn: Burning period for MCMC chain
 * Counts: $q \times n$ matrix of expression counts
 * BatchDesign: Design matrix representing batch information
 * (number of columns must be equal to number of batches)
 * mu0: Starting value of $\mu=(\mu_1,...,\mu_q_0)'$
 * delta0: Starting value of $\delta=(\delta_1,...,\delta_{q_0})'$
 * phi0: Starting value of $\phi=(\phi_1,...,\phi_n)$'
 * nu0: Starting value of $\nu=(\nu_1,...,\nu_n)$'
 * theta0: Starting value of $\theta$
 * s2mu: Prior variance for log-Normal(0, $s^2_{\mu}$) assigned to the
 * unconstrained mean expression parameters
 * adelta: Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$)
 * prior assigned to each $\delta_i$
 * bdelta: Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$)
 * prior assigned to each $\delta_i$
 * s2delta: Prior variance for log-Normal(0, $s^2_{\delta}$) assigned to each
 * $\delta_i$
 * prior_delta: (as in HiddenBASiCS_MCMCcpp)
 * atheta: Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$)
 * prior assigned to $\theta$
 * btheta: Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$)
 * prior assigned to $\theta$
 * ar: Optimal acceptance rate for adaptive Metropolis-Hastings updates
 * LSmu0: Starting value of adaptive proposal variance of
 * $\mu=(\mu_1,...,\mu_q_0)'$ (log-scale)
 * LSdelta0: Starting value of adaptive proposal variance of
 * $\delta=(\delta_1,...,\delta_{q_0})'$ (log-scale)
 * LSnu0: Starting value of adaptive proposal variance of
 * $\nu=(\nu_1,...,\nu_n)'$ (log-scale)
 * LStheta0: Starting value of adaptive proposal variance of $\theta$ (log-scale)
 * sumByCellAll: Sum of expression counts by cell (all genes)
 * sumByGeneAll: Sum of expression counts by gene (all genes)
 * StoreAdapt: (as in HiddenBASiCS_MCMCcpp)
 * EndAdapt: (as in HiddenBASiCS_MCMCcpp)
 * PrintProgress: (as in HiddenBASiCS_MCMCcpp)
 * Constrain:
 * Index:
 * RefGene:
 * RefGenes:
 * ConstrainGene:
 * NotConstrainGene:
 * ConstrainType:
 * StochasticRef:
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcppNoSpikes(
    int N,
    int Thin,
    int Burn,
    NumericMatrix Counts,
    NumericMatrix BatchDesign,
    NumericVector mu0,
    NumericVector delta0,
    NumericVector s0,
    NumericVector nu0,
    NumericVector theta0,
    double s2mu,
    double adelta,
    double bdelta,
    double s2delta,
    double prior_delta,
    double as,
    double bs,
    double atheta,
    double btheta,
    double ar,
    NumericVector LSmu0,
    NumericVector LSdelta0,
    NumericVector LSnu0,
    NumericVector LStheta0,
    NumericVector sumByCellAll,
    NumericVector sumByGeneAll,
    int StoreAdapt,
    int EndAdapt,
    int PrintProgress,
    double Constrain,
    NumericVector Index,
    int RefGene,
    NumericVector RefGenes,
    NumericVector ConstrainGene,
    NumericVector NotConstrainGene,
    int ConstrainType,
    int StochasticRef,
    double geneExponent,
    double cellExponent)
{
  using arma::ones;
  using arma::zeros;
  using Rcpp::Rcout;

  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol();
  int q0 = Counts.nrow();
  int Naux = N/Thin - Burn/Thin;
  int nBatch = BatchDesign.ncol();

  // TRANSFORMATION TO ARMA ELEMENTS
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll);
  arma::mat Counts_arma = as_arma(Counts);
  arma::mat BatchDesign_arma = as_arma(BatchDesign);
  arma::vec Index_arma = as_arma(Index);
  arma::uvec ConstrainGene_arma = Rcpp::as<arma::uvec>(ConstrainGene);
  arma::uvec NotConstrainGene_arma = Rcpp::as<arma::uvec>(NotConstrainGene);
  arma::vec RefGenes_arma = as_arma(RefGenes);

  // OTHER GLOBAL QUANTITIES
  arma::vec BatchSizes = sum(BatchDesign_arma,0).t();

  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = zeros(Naux, q0);
  arma::mat delta = zeros(Naux, q0);
  arma::mat s = ones(Naux, n);
  arma::mat nu = zeros(Naux, n);
  arma::mat theta = zeros(Naux, nBatch);
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::mat LSnu;
  arma::mat LStheta;

  // LOG-PROPOSAL VARIANCES
  if (StoreAdapt == 1) {
    LSmu = zeros(Naux, q0);
    LSdelta = zeros(Naux, q0);
    LSnu = zeros(Naux, n);
    LStheta = zeros(Naux, nBatch);
  }

  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = zeros(q0);
  arma::vec deltaAccept = zeros(q0);
  arma::vec nuAccept = zeros(n);
  arma::vec thetaAccept = zeros(nBatch);

  arma::vec PmuAux = zeros(q0);
  arma::vec PdeltaAux = zeros(q0);
  arma::vec PnuAux = zeros(n);
  arma::vec PthetaAux = zeros(nBatch);

  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = zeros(q0,2);
  muAux.col(0) = as_arma(mu0);
  arma::mat deltaAux = zeros(q0,2);
  deltaAux.col(0) = as_arma(delta0);
  arma::vec sAux = as_arma(s0);
  arma::mat nuAux = zeros(n,2);
  nuAux.col(0) = as_arma(nu0);
  arma::mat thetaAux = zeros(nBatch, 2);
  thetaAux.col(0) = as_arma(theta0);
  arma::vec thetaBatch = BatchDesign_arma * as_arma(theta0);

  // INITIALIZATION OF ADAPTIVE VARIANCES
  arma::vec LSmuAux = as_arma(LSmu0);
  arma::vec LSdeltaAux = as_arma(LSdelta0);
  arma::vec LSnuAux = as_arma(LSnu0);
  arma::vec LSthetaAux = as_arma(LStheta0);

  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q0);
  arma::vec PdeltaAux0 = arma::zeros(q0);
  arma::vec PnuAux0 = arma::zeros(n);
  arma::vec PthetaAux0 = arma::ones(nBatch);

  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES
  // (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;

  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  // To avoid repeated initialisation
  arma::vec y_q0 = arma::ones(q0);
  arma::vec y_n = arma::ones(n);
  arma::vec ind_q0 = arma::zeros(q0);
  arma::vec ind_n = arma::zeros(n);
  arma::vec u_q0 = arma::zeros(q0);
  arma::vec u_n = arma::zeros(n);

  // INITIALIZATION OF PARAMETERS RELATED TO STOCHASTIC REF
  arma::vec RefFreq = arma::zeros(q0);
  int RefAux;

  StartSampler(N);

  // START OF MCMC LOOP
  for (i = 0; i < N; i++) {

    Rcpp::checkUserInterrupt();

    if (i == Burn) EndBurn();

    Ibatch++;

    // UPDATE OF PHI
    // WE CAN RECYCLE THE SAME FULL CONDITIONAL AS IMPLEMENTED FOR S (BATCH CASE)
    sAux = sUpdateBatch(sAux,
                        nuAux.col(0),
                        thetaBatch,
                        as,
                        bs,
                        BatchDesign_arma,
                        n,
                        y_n,
                        cellExponent);

    // UPDATE OF THETA:
    // 1st ELEMENT IS THE UPDATE,
    // 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdateBatch(thetaAux.col(0),
                                exp(LSthetaAux),
                                BatchDesign_arma,
                                BatchSizes,
                                sAux,
                                nuAux.col(0),
                                atheta,
                                btheta,
                                n,
                                nBatch,
                                cellExponent);
    PthetaAux += thetaAux.col(1);
    if (i >= Burn) {
      thetaAccept += thetaAux.col(1);
    }
    thetaBatch = BatchDesign_arma * thetaAux.col(0);

    // UPDATE OF MU:
    // 1st COLUMN IS THE UPDATE,
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    // If using stochastic reference, randomly select 1 ref gene
    if (StochasticRef == 1) {
      RefAux = as_scalar(
        arma::randi(1, arma::distr_param(0, RefGenes_arma.size() - 1))
      );
      RefGene = RefGenes(RefAux);
      if (i >= Burn) RefFreq(RefGene) += 1;
    }
    muAux = muUpdateNoSpikes(muAux.col(0),
                             exp(LSmuAux),
                             Counts_arma,
                             1 / deltaAux.col(0),
                             nuAux.col(0),
                             sumByCellAll_arma,
                             s2mu,
                             q0,
                             n,
                             y_q0,
                             u_q0,
                             ind_q0,
                             Constrain,
                             RefGene,
                             ConstrainGene_arma,
                             NotConstrainGene_arma,
                             ConstrainType);

    PmuAux += muAux.col(1);
    if (i >= Burn) {
      muAccept += muAux.col(1);
    }

    // UPDATE OF DELTA:
    // 1st COLUMN IS THE UPDATE,
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0),
                           exp(LSdeltaAux),
                           Counts_arma,
                           muAux.col(0),
                           nuAux.col(0),
                           adelta,
                           bdelta,
                           s2delta,
                           prior_delta,
                           q0,
                           n,
                           y_q0,
                           u_q0,
                           ind_q0,
                           geneExponent);
    PdeltaAux += deltaAux.col(1);
    if (i>=Burn) {
      deltaAccept += deltaAux.col(1);
    }

    // UPDATE OF NU:
    // 1st COLUMN IS THE UPDATE,
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateBatchNoSpikes(nuAux.col(0),
                                  exp(LSnuAux),
                                  Counts_arma,
                                  BatchDesign_arma,
                                  muAux.col(0),
                                  1 / deltaAux.col(0),
                                  sAux,
                                  thetaBatch,
                                  sumByGeneAll_arma,
                                  q0,
                                  n,
                                  y_n,
                                  u_n,
                                  ind_n);

    PnuAux += nuAux.col(1);
    if (i >= Burn) {
      nuAccept += nuAux.col(1);
    }

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if (i < EndAdapt) {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if (Ibatch == 50) {
        PmuAux = PmuAux/(50-RefFreq);
        PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux.elem(find(Index_arma != RefGene)) = LSmuAux.elem(find(Index_arma != RefGene)) + PmuAux.elem(find(Index_arma != RefGene))*0.1;
        PdeltaAux = PdeltaAux/50;
        PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux = LSdeltaAux+PdeltaAux*0.1;
        PnuAux = PnuAux/50;
        PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux = LSnuAux+PnuAux*0.1;
        PthetaAux = PthetaAux/50;
        PthetaAux = -1+2*arma::conv_to<arma::mat>::from(PthetaAux>ar);
        LSthetaAux = LSthetaAux + PthetaAux*0.1;

        Ibatch = 0;
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0;
        PnuAux = PnuAux0; PthetaAux = PthetaAux0;
      }
    }

    // STORAGE OF DRAWS
    if ((i % Thin == 0) & (i >= Burn)) {
      int ind = i / Thin - Burn / Thin;
      mu.row(ind) = muAux.col(0).t();
      delta.row(ind) = deltaAux.col(0).t();
      s.row(ind) = sAux.t();
      nu.row(ind) = nuAux.col(0).t();
      theta.row(ind) = thetaAux.col(0).t();

      if (StoreAdapt == 1) {
        LSmu.row(ind) = LSmuAux.t();
        LSdelta.row(ind) = LSdeltaAux.t();
        LSnu.row(ind) = LSnuAux.t();
        LStheta.row(ind) = LSthetaAux.t();
      }
    }

    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if ((i % (2 * Thin) == 0) & (PrintProgress == 1)) {
      CurrentIter(i, N);
      Rcpp::Rcout << "mu (gene 1): " << muAux(0,0) << std::endl;
      Rcpp::Rcout << "delta (gene 1): " << deltaAux(0,0) << std::endl;
      Rcpp::Rcout << "s (cell 1): " << sAux(0) << std::endl;
      Rcpp::Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
      Rcpp::Rcout << "theta (batch 1): " << thetaAux(0,0) << std::endl;
      Rcout << "-----------------------------------------------------" << std::endl;
      Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
      Rcpp::Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
      Rcpp::Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl;
      Rcpp::Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
      Rcpp::Rcout << "LStheta (batch 1): " << LSthetaAux(0) << std::endl;
    }
  }

  // END OF MCMC SAMPLER AND ACCEPTANCE RATE CONSOLE OUTPUT
  EndSampler(N);
  ReportAR(muAccept/(N-Burn), "mu[i]'s");
  ReportAR(deltaAccept/(N-Burn), "delta[i]'s");
  ReportAR(nuAccept/(N-Burn), "nu[jk]'s");
  ReportAR(thetaAccept/(N-Burn), "theta[k]'s");

  Rcout << "-----------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;

  if (StoreAdapt == 1) {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu") = mu,
        Rcpp::Named("delta") = delta,
        Rcpp::Named("s") = s,
        Rcpp::Named("nu") = nu,
        Rcpp::Named("theta") = theta,
        Rcpp::Named("ls.mu") = LSmu,
        Rcpp::Named("ls.delta") = LSdelta,
        Rcpp::Named("ls.nu") = LSnu,
        Rcpp::Named("ls.theta") = LStheta,
        Rcpp::Named("RefFreq") = RefFreq / (N - Burn)));
  }
  else {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu") = mu,
        Rcpp::Named("delta") = delta,
        Rcpp::Named("s") = s,
        Rcpp::Named("nu") = nu,
        Rcpp::Named("theta") = theta,
        Rcpp::Named("RefFreq") = RefFreq / (N - Burn)));
  }
}
