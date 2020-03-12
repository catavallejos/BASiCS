#include "utils.h"

arma::mat BASiCS_DenoisedRates(
    NumericMatrix CountsBio, 
    NumericMatrix Mu,
    NumericMatrix TransInvDelta,
    NumericMatrix PhiNu, 
    int N,
    int q0,
    int n) {

  // Transformations to arma objects
  arma::mat CountsBio_arma = as_arma(CountsBio);
  arma::mat Mu_arma = as_arma(Mu);
  arma::mat TransInvDelta_arma = as_arma(TransInvDelta);
  arma::mat PhiNu_arma = as_arma(PhiNu);
  
  // Where to store the results
  arma::mat Rho = arma::zeros(q0, n);
  // Auxiliary matrices
  arma::mat m1; arma::mat m2;
  
  for (int i = 0; i<N; i++) {
    
    Rcpp::checkUserInterrupt();
    
    m1 = CountsBio_arma; 
    m1.each_col() += TransInvDelta_arma.col(i); 
    m2 = Mu_arma.row(i).t() * PhiNu_arma.row(i);
    m2.each_col() += TransInvDelta_arma.col(i); 
    Rho += m1 / m2; 
  } 
  return(Rho / N);
}


  
  
  

