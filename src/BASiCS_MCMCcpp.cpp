#include "utils.h"

/* MCMC sampler
 * N: Total number of MCMC draws
 * Thin: Thinning period for MCMC chain
 * Burn: Burning period for MCMC chain
 * Counts: Matrix of expression counts
 * BatchDesign: Design matrix representing batch information
 * muSpikes: mu values for spike-in genes
 * mu0: Starting value for mu
 * delta0: Starting value for delta
 * phi0: Starting value for phi
 * s0: Starting value for s
 * nu0: Starting value for nu
 * theta0: Starting value for theta
 * s2mu: prior variance for mu
 * adelta: prior shape for delta (when using a gamma prior)
 * bdelta: prior rate for delta (when using a gamma prior)
 * s2delta: prior variance for delta (when using a log-normal prior)
 * prior_delta: whether gamma or log-normal prior is used for delta
 * aphi: Dirichlet hyper-parameter for phi
 * as: prior shape for s
 * bs: prior rate for s
 * atheta: prior shape for theta
 * btheta: prior rate for theta
 * ar: Optimal acceptance rate for adaptive Metropolis-Hastings updates
 * LSmu0: Starting value of adaptive proposal variance of mu (log-scale)
 * LSdelta0: Starting value of adaptive proposal variance of delta (log-scale)
 * LSphi0: Starting value of adaptive proposal precision of phi (log-scale)
 * LSnu0: Starting value of adaptive proposal variance of nu (log-scale)
 * LStheta0: Starting value of adaptive proposal variance of theta (log-scale)
 * sumByCellAll: Sum of expression counts by cell (all genes)
 * sumByCellBio:  Sum of expression counts by cell (biological genes only)
 * sumByGeneAll: Sum of expression counts by gene (all genes)
 * sumByGeneBio: Sum of expression counts by gene (biological genes only)
 * StoreAdapt: whether to store adaptive variances
 * EndAdapt: when to stop the adaptation
 * PrintProgress: whether to print progress report
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcpp(
    int N,
    int Thin,
    int Burn,
    NumericMatrix Counts,
    NumericMatrix BatchDesign,
    NumericVector muSpikes,
    NumericVector mu0,
    NumericVector delta0,
    NumericVector phi0,
    NumericVector s0,
    NumericVector nu0,
    NumericVector theta0,
    double s2mu,
    double adelta,
    double bdelta,
    double s2delta,
    double prior_delta,
    NumericVector aphi,
    double as,
    double bs,
    double atheta,
    double btheta,
    double ar,
    NumericVector LSmu0,
    NumericVector LSdelta0,
    double LSphi0,
    NumericVector LSnu0,
    NumericVector LStheta0,
    NumericVector sumByCellAll,
    NumericVector sumByCellBio,
    NumericVector sumByGeneAll,
    NumericVector sumByGeneBio,
    int StoreAdapt,
    int EndAdapt,
    int PrintProgress,
    double geneExponent,
    double cellExponent)
{
  using arma::ones;
  using arma::zeros;
  using Rcpp::Rcout;

  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol();
  int nBatch = BatchDesign.ncol();
  int q0 = delta0.size();
  int Naux = N / Thin - Burn / Thin;

  // TRANSFORMATION TO ARMA ELEMENTS
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll);
  arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll);
  arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts);
  arma::mat BatchDesign_arma = as_arma(BatchDesign);

  // OTHER GLOBAL QUANTITIES
  double SumSpikeInput = sum(muSpikes);
  arma::vec BatchSizes = sum(BatchDesign_arma, 0).t();

  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = zeros(q0, Naux);
  arma::mat delta = zeros(q0, Naux);
  arma::mat phi = ones(n, Naux);
  arma::mat s = zeros(n, Naux);
  arma::mat nu = zeros(n, Naux);
  arma::mat theta = zeros(nBatch, Naux);
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::vec LSphi;
  arma::mat LSnu;
  arma::mat LStheta;

  // LOG-PROPOSAL VARIANCES
  if (StoreAdapt == 1) {
    LSmu = zeros(q0, Naux);
    LSdelta = zeros(q0, Naux);
    LSphi = ones(Naux);
    LSnu = zeros(n, Naux);
    LStheta = zeros(nBatch, Naux);
  }

  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = zeros(q0);
  arma::vec PmuAux = zeros(q0);
  arma::vec deltaAccept = zeros(q0);
  arma::vec PdeltaAux = zeros(q0);
  double phiAccept = 0;
  double PphiAux = 0;
  arma::vec nuAccept = zeros(n);
  arma::vec PnuAux = zeros(n);
  arma::vec thetaAccept = zeros(nBatch);
  arma::vec PthetaAux = zeros(nBatch);

  // INITIALIZATION OF PARAMETER VALUES FOR MCMC RUN
  arma::mat muAux = zeros(q0, 2);
  muAux.col(0) = as_arma(mu0);
  arma::mat deltaAux = zeros(q0, 2);
  deltaAux.col(0) = as_arma(delta0);
  arma::vec phiAux = as_arma(phi0);
  Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0);
  arma::mat nuAux = zeros(n, 2);
  nuAux.col(0) = as_arma(nu0);
  arma::mat thetaAux = zeros(nBatch, 2);
  thetaAux.col(0) = as_arma(theta0);
  arma::vec thetaBatch = BatchDesign_arma * as_arma(theta0);

  // INITIALIZATION OF ADAPTIVE VARIANCES
  arma::vec LSmuAux = as_arma(LSmu0);
  arma::vec LSdeltaAux = as_arma(LSdelta0);
  double LSphiAux = LSphi0;
  arma::vec LSnuAux = as_arma(LSnu0);
  arma::vec LSthetaAux = as_arma(LStheta0);

  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = zeros(q0);
  arma::vec PdeltaAux0 = zeros(q0);
  double PphiAux0 = 0;
  arma::vec PnuAux0 = zeros(n);
  arma::vec PthetaAux0 = zeros(nBatch);

  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES
  // (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0;

  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  // To avoid repeated initialisation
  arma::vec y_q0 = ones(q0);
  arma::vec y_n = ones(n);
  arma::vec ind_q0 = zeros(q0);
  arma::vec ind_n = zeros(n);
  arma::vec u_q0 = zeros(q0);
  arma::vec u_n = zeros(n);

  StartSampler(N);

  // START OF MCMC LOOP
  for (int i=0; i<N; i++) {

    Rcpp::checkUserInterrupt();

    if (i == Burn) EndBurn();

    Ibatch++;

    // UPDATE OF PHI:
    // 1st ELEMENT IS THE UPDATE,
    // 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    phiAuxList = phiUpdate(phiAux,
                           exp(LSphiAux),
                           Counts_arma,
                           muAux.col(0),
                           1 / deltaAux.col(0),
                           nuAux.col(0),
                           aphi,
                           sumByGeneBio_arma,
                           q0,
                           n,
                           y_n,
                           cellExponent);

    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]);
    if (i >= Burn) {
      phiAccept += Rcpp::as<double>(phiAuxList["ind"]);
    }

    // theta and regression parameters are global; thus, they require scaling
    // for both cases
    double exponent = 1;
    if (geneExponent != 1) {
      exponent = geneExponent;
    } else if (cellExponent != 1) {
      exponent = cellExponent;
    }

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
                                exponent);
    PthetaAux += thetaAux.col(1);
    if (i >= Burn) {
      thetaAccept += thetaAux.col(1);
    }
    thetaBatch = BatchDesign_arma * thetaAux.col(0);

    // UPDATE OF MU:
    // 1st COLUMN IS THE UPDATE,
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    muAux = Hidden_muUpdate(muAux.col(0),
                            exp(LSmuAux),
                            Counts_arma,
                            1 / deltaAux.col(0),
                            phiAux % nuAux.col(0),
                            sumByCellBio_arma,
                            s2mu,
                            q0,
                            n,
                            y_q0,
                            u_q0,
                            ind_q0,
                            geneExponent);
    PmuAux += muAux.col(1);
    if (i >= Burn) {
      muAccept += muAux.col(1);
    }

    // UPDATE OF S
    sAux = sUpdateBatch(sAux,
                        nuAux.col(0),
                        thetaBatch,
                        as,
                        bs,
                        BatchDesign_arma,
                        n,
                        y_n,
                        cellExponent);

    // UPDATE OF DELTA:
    // 1st COLUMN IS THE UPDATE,
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0),
                           exp(LSdeltaAux),
                           Counts_arma,
                           muAux.col(0),
                           phiAux % nuAux.col(0),
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
    if (i >= Burn) {
      deltaAccept += deltaAux.col(1);
    }

    // UPDATE OF NU:
    // 1st COLUMN IS THE UPDATE,
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateBatch(nuAux.col(0),
                          exp(LSnuAux),
                          Counts_arma,
                          SumSpikeInput,
                          BatchDesign_arma,
                          muAux.col(0),
                          1 / deltaAux.col(0),
                          phiAux,
                          sAux,
                          thetaBatch,
                          sumByGeneAll_arma,
                          q0,
                          n,
                          y_n,
                          u_n,
                          ind_n,
                          cellExponent);
    PnuAux += nuAux.col(1);
    if (i >= Burn) {
      nuAccept += nuAux.col(1);
    }

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if (i < EndAdapt) {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if (Ibatch == 50) {
        PmuAux = PmuAux / 50;
        PmuAux = -1 + 2 * arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux = LSmuAux + PmuAux * 0.1;
        PdeltaAux = PdeltaAux / 50;
        PdeltaAux = -1 + 2 * arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux = LSdeltaAux + PdeltaAux * 0.1;
        PphiAux = PphiAux / 50;
        PphiAux = -1 + 2 * (PphiAux>ar);
        LSphiAux = LSphiAux - PphiAux*0.1;
        PnuAux = PnuAux/50;
        PnuAux = -1 + 2 * arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux = LSnuAux + PnuAux * 0.1;
        PthetaAux = PthetaAux/50;
        PthetaAux = -1 + 2 * arma::conv_to<arma::mat>::from(PthetaAux>ar);
        LSthetaAux = LSthetaAux + PthetaAux * 0.1;

        Ibatch = 0;
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0;
        PphiAux = PphiAux0;
        PnuAux = PnuAux0; PthetaAux = PthetaAux0;
      }
    }

    // STORAGE OF DRAWS
    if ((i % Thin == 0) & (i >= Burn)) {
      int ind = i / Thin - Burn / Thin;
      mu.col(ind) = muAux.col(0);
      delta.col(ind) = deltaAux.col(0);
      phi.col(ind) = phiAux;
      s.col(ind) = sAux;
      nu.col(ind) = nuAux.col(0);
      theta.col(ind) = thetaAux.col(0);

      if (StoreAdapt == 1) {
        LSmu.col(ind) = LSmuAux;
        LSdelta.col(ind) = LSdeltaAux;
        LSphi(ind) = LSphiAux;
        LSnu.col(ind) = LSnuAux;
        LStheta.col(ind) = LSthetaAux;
      }
    }

    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if ((i % (2 * Thin) == 0) & (PrintProgress == 1)) {
      CurrentIter(i, N);
      Rcout << "mu (gene 1): " << muAux(0,0) << std::endl;
      Rcout << "delta (gene 1): " << deltaAux(0,0) << std::endl;
      Rcout << "phi (cell 1): " << phiAux(0) << std::endl;
      Rcout << "s (cell 1): " << sAux(0) << std::endl;
      Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
      Rcout << "theta (batch 1): " << thetaAux(0,0) << std::endl;
      Rcout << "-----------------------------------------------------" << std::endl;
      Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
      Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
      Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl;
      Rcout << "LSphi: " << LSphiAux << std::endl;
      Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
      Rcout << "LStheta (batch 1): " << LSthetaAux(0) << std::endl;
    }
  }

  // END OF MCMC SAMPLER AND ACCEPTANCE RATE CONSOLE OUTPUT
  EndSampler(N);
  ReportAR(muAccept / (N - Burn), "mu[i]'s");
  ReportAR(deltaAccept / (N - Burn), "delta[i]'s");
  Rcout << " " << std::endl;
  Rcout << "Acceptance rate for phi (joint): " <<
    phiAccept / (N - Burn) << std::endl;
  Rcout << " " << std::endl;
  ReportAR(nuAccept / (N - Burn), "nu[j]'s");
  ReportAR(thetaAccept / (N - Burn), "theta[k]'s");

  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;

  if (StoreAdapt == 1) {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu") = mu.t(),
        Rcpp::Named("delta") = delta.t(),
        Rcpp::Named("phi") = phi.t(),
        Rcpp::Named("s") = s.t(),
        Rcpp::Named("nu") = nu.t(),
        Rcpp::Named("theta") = theta.t(),
        Rcpp::Named("ls.mu") = LSmu.t(),
        Rcpp::Named("ls.delta") = LSdelta.t(),
        Rcpp::Named("ls.phi") = LSphi,
        Rcpp::Named("ls.nu") = LSnu.t(),
        Rcpp::Named("ls.theta") = LStheta.t()));
  }

  else {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu") = mu.t(),
        Rcpp::Named("delta") = delta.t(),
        Rcpp::Named("phi") = phi.t(),
        Rcpp::Named("s") = s.t(),
        Rcpp::Named("nu") = nu.t(),
        Rcpp::Named("theta") = theta.t()));
  }
}


