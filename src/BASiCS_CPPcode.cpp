/* C++ implementation of BASiCS 
 * Author: Catalina A. Vallejos (cnvallej@uc.cl)
 */

#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
//File: matprod_arma.cpp
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

#define ZTOL sqrt(DOUBLE_EPS)

/* 
 * Rgig is an adaptation of the rgig.c function implemented by 
 * Ester Pantaleo and Robert B. Gramacy, 2010
 * (which in turn, was adapted from the C code in the ghyp 
 * v_1.5.2 package for R, rgig function)
 * Our adaptation makes the function compatible with the Rcpp classes
 */

/* Auxiliary function taken from rgig.c */
double gig_y_gfn(double y, double m, double beta, double lambda)
{  
/* 
 * gig_y_gfn: 
 *
 * evaluate the function that we need to find the root 
 * of in order to construct the optimal rejection
 * sampler to obtain GIG samples
 */
  double y2, g;
  y2 = y * y;
  g = 0.5 * beta * y2 * y;
  g -= y2 * (0.5 * beta * m + lambda + 1.0);
  g += y * ((lambda - 1.0) * m - 0.5 * beta) + 0.5 * beta * m;
  return(g);
}

/* Auxiliary function taken from rgig.c */
double rinvgauss(const double mu, const double lambda)
{
/*
 * rinvgauss:
 *
 * Michael/Schucany/Haas Method for generating Inverse Gaussian
 * random variable with mean mu and shape parameter lambda, 
 * as given in Gentle's book on page 193
 */
  double u, y, x1, mu2, l2;

  y = norm_rand(); /// ******
  y *= y;
  mu2 = mu * mu;
  l2 = 2.0*lambda;
  x1 = mu + mu2*y/l2 - (mu/l2)* sqrt(4.0*mu*lambda*y + mu2*y*y);

  u = unif_rand(); // ******
  if(u <= mu/(mu + x1)) return x1;
  else return mu2/x1;
}

/* Auxiliary function taken from rgig.c */
double zeroin_gig(double ax,double bx,
                  double f(double x, double m, double beta, double lambda),
                  double tol,double m,double beta,double lambda)  
/* An estimate to the root  */
//double ax;				/* Left border | of the range	*/
//double bx;  				/* Right border| the root is seeked*/
/* Function under investigation	*/
//double (*f)(double x, double m, double beta, double lambda);	
//double tol;				/* Acceptable tolerance	*/
//double m;                               /* specific to gig_y_gfn */
//double beta;                            /* specific to gig_y_gfn */
//double lambda;                          /* specific to gig_y_gfn */
{
 
  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/

  a = ax;  b = bx;  fa = (f)(a, m, beta, lambda);  fb = (f)(b, m, beta, lambda);
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b-a;		/* Distance from the last but one*/
					/* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  					/* sion operations is delayed   */
 					/* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
   
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
	    a = b;  b = c;  c = a;      /* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2.0*DOUBLE_EPS*fabs(b) + tol/2.0;
    new_step = (c-b)/2.0;

    if( fabs(new_step) <= tol_act || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/

    /* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	        && fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	    register double t1,cb,t2;
	    cb = c-b;
	    if( a==c )			/* If we have only two distinct	*/
	    {				/* points linear interpolation 	*/
	      t1 = fb/fa;			/* can only be applied		*/
	      p = cb*t1;
	      q = 1.0 - t1;
 	    }
	    else				/* Quadric inverse interpolation*/
	    {
	      q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
	      p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	      q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	  
	    if( p>(double)0 )		/* p was calculated with the op-*/
	      q = -q;			/* posite sign; make p positive	*/
	    else				/* and assign possible minus to	*/
	      p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2.0)	/* If b+p/q falls in [b,c]*/
	          && p < fabs(prev_step*q/2.0) )	/* and isn't too large	*/
	      new_step = p/q;			/* it is accepted	*/
					/* If p/q is too large then the	*/
					/* bissection procedure can 	*/
					/* reduce [b,c] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act ) 
    {	/* Adjust the step to be not less*/
      if( new_step > (double)0 )	/* than tolerance		*/
	      new_step = tol_act;
      else
	      new_step = -tol_act;
    }

    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;  fb = (*f)(b, m, beta, lambda);  /* Do step to a new approxim. */
    if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) )
    {                 			/* Adjust c for it to have a sign*/
      c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }
}

/* Random sample generator from a Generalized Inverse Gaussian distribution 
 * (modified version of rgig.c using Rcpp classes)
 */
NumericVector Rgig(const int n, 
                   const double lambda, 
                   const double chi, 
                   const double psi) 
{
  NumericVector samps(n);
  /* special case which is basically a gamma distribution */
  if((chi < ZTOL) & (lambda > 0.0)) {
    int i;
    for(i=0; i<n; i++) samps(i) = R::rgamma(lambda, 2.0/psi);
    return samps;
  }

  /* special cases which is basically an inverse gamma distribution */
  if((psi < ZTOL) & (lambda < 0.0)) { 
    int i;
    for(i=0; i<n; i++) samps(i) = 1.0/R::rgamma(0.0-lambda, 2.0/chi);
    return samps;
  }

  /* special case which is basically an inverse gaussian distribution */
  if(lambda == -0.5) {
    double alpha;
    int i;
    alpha = sqrt(chi/psi);
    for(i=0; i<n; i++) samps(i) = rinvgauss(alpha, chi);
    return samps;
  }  

  /*
   * begin general purpose rgig code, which was basically 
   * translated from the R function rgig in the ghyp package v_1.5.2
   */

  double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2, Y;
  int i, need;

  alpha = sqrt(chi/psi);
  beta2 = psi*chi;
  beta = sqrt(psi*chi);
  lm1 = lambda - 1.0;
  lm12 = lm1*lm1;
  m = (lm1 + sqrt(lm12 + beta2))/beta;
  m1 = m + 1.0/m;
    
  upper = m;
  while (gig_y_gfn(upper, m, beta, lambda) <= 0) { upper *= 2.0; }

  yM = zeroin_gig(0.0, m, gig_y_gfn, ZTOL, m, beta, lambda);
  yP = zeroin_gig(m, upper, gig_y_gfn, ZTOL, m, beta, lambda);
  
  a = (yP - m) * pow(yP/m, 0.5 * lm1);
  a *=  exp(-0.25 * beta * (yP + 1.0/yP - m1));
  b = (yM - m) * pow(yM/m, 0.5 * lm1);
  b *= exp(-0.25 * beta * (yM + 1/yM - m1));
  c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);

  for (i=0; i<n; i++) 
  {
    need = 1;
    while (need) 
    {
      R1 = unif_rand();
      R2 = unif_rand();

      Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
      if (Y > 0.0) 
      {
	      if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) 
	      {
	      need = 0;
	      }
      }
    }
    samps[i] = Y*alpha;
  }
  return(samps);
}

/* Random sample generator from a Generalized Inverse Gaussian distribution 
 * (modified version of rgig.c using Rcpp classes)
 */
double RgigDouble(const double lambda, const double chi, const double psi)
{
  double samps;
  /* special case which is basically a gamma distribution */
  if((chi < ZTOL) & (lambda > 0.0)) {
    samps = R::rgamma(lambda, 2.0/psi);
    return samps;
  }

  /* special cases which is basically an inverse gamma distribution */
  if((psi < ZTOL) & (lambda < 0.0)) { 
    samps = 1.0/R::rgamma(0.0-lambda, 2.0/chi);
    return samps;
  }

  /* special case which is basically an inverse gaussian distribution */
  if(lambda == -0.5) {
    double alpha;
    alpha = sqrt(chi/psi);
    samps = rinvgauss(alpha, chi);
    return samps;
  }  

  /*
   * begin general purpose rgig code, which was basically 
   * translated from the R function rgig in the ghyp package v_1.5.2
   */

  double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2, Y;
  int need;

  alpha = sqrt(chi/psi);
  beta2 = psi*chi;
  beta = sqrt(psi*chi);
  lm1 = lambda - 1.0;
  lm12 = lm1*lm1;
  m = (lm1 + sqrt(lm12 + beta2))/beta;
  m1 = m + 1.0/m;
    
  upper = m;
  while (gig_y_gfn(upper, m, beta, lambda) <= 0) { upper *= 2.0; }

  yM = zeroin_gig(0.0, m, gig_y_gfn, ZTOL, m, beta, lambda);
  yP = zeroin_gig(m, upper, gig_y_gfn, ZTOL, m, beta, lambda);
  
  a = (yP - m) * pow(yP/m, 0.5 * lm1);
  a *=  exp(-0.25 * beta * (yP + 1.0/yP - m1));
  b = (yM - m) * pow(yM/m, 0.5 * lm1);
  b *= exp(-0.25 * beta * (yM + 1/yM - m1));
  c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);

  need = 1;
  while (need) 
  {
    R1 = unif_rand();
    R2 = unif_rand();

    Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
    if (Y > 0.0) 
    {
      if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) 
      {
        need = 0;
	    }
    }
  }
  samps = Y*alpha;
  return(samps);
}

/* END OF ADAPTED CODE*/

/* Auxiliary function that converts Rcpp::NumericVector 
 * objects into arma::vec objects
 */ 
arma::vec as_arma(NumericVector& x) 
{
  return arma::vec(x.begin(), x.length(), false);
}

/* Auxiliary function that converts Rcpp::NumericMatrix 
 * objects into arma::mac objects
 */ 
arma::mat as_arma(NumericMatrix& x) 
{
  return arma::mat(x.begin(), x.nrow(), x.ncol(), false);
}

/* Auxiliary function: logarithm of a Gamma function applied  
 * element-wise to an arma::mat object
 */
arma::mat lgamma_cpp(arma::mat const& x) 
{
  arma::mat output = x;
  for (unsigned int i=0; i<arma::size(x,0); i++)
  {
    for(unsigned int j=0; j<arma::size(x,1); j++)
    {
      output(i,j) = R::lgammafn(x(i,j));
    }
  }
  return output;
}

/* Auxiliary function: logarithm of a Gamma function applied 
 * element-wise to an arma::mat object
 */ 
arma::vec lgamma_cpp_vec(arma::vec const& x) 
{
  arma::vec output = x;
  for (unsigned int i=0; i < output.size(); i++)
  {
      output(i) = R::lgammafn(x(i));
  }
  return output;
}

/* Auxiliary function: to avoid numerical overflows/underflows
 * when computing a sum of exponentiated values
 */  
double log_sum_exp_cpp(arma::vec const& x_arma) 
{
  double offset;
  if ( max(abs(x_arma)) > max(x_arma) ) { offset = min(x_arma);}
  else { offset = max(x_arma);}
  return log(sum(exp(x_arma - offset))) + offset; 
}

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes
 * mu0: current value of mu
 * prop_var: current value of the proposal variances for mu
 * Counts: matrix of expression counts
 * delta: current value of delta
 * phi: current value of phi
 * nu: current value of nu
 * sum_bycell_bio: sum of expression counts by cell (biological genes only)
 * s2_mu: prior hyper-variance for mu
 * q0: number of biological genes
 * n: number of cells
 * mu: auxiliary vector for storage
 * ind: auxiliary vector for storage
 */
arma::mat muUpdate(
  arma::vec const& mu0, 
  arma::vec const& prop_var, 
  arma::mat const& Counts, 
  arma::vec const& invdelta, 
  arma::vec const& phinu, 
  arma::vec const& sum_bycell_bio, 
  double const& s2_mu,
  int const& q0,
  int const& n,
  arma::vec & mu1,
  arma::vec & u, 
  arma::vec & ind)
{
  
  /* PROPOSAL STEP */
  mu1 = exp( arma::randn(q0) % sqrt(prop_var) + log(mu0) );
  u = arma::randu(q0);
    
  /* ACCEPT/REJECT STEP 
   * Note: there is a -1 factor coming from the log-normal prior. 
   * However, it cancels out as using log-normal proposals.
   */
  arma::vec log_aux = (log(mu1) - log(mu0)) % sum_bycell_bio; 
  for (int i=0; i < q0; i++)
  {
    for (int j=0; j < n; j++) 
    {
      log_aux(i) -= ( Counts(i,j) + invdelta(i) ) *  
                      log( ( phinu(j)*mu1(i) + invdelta(i) ) / 
                           ( phinu(j)*mu0(i) + invdelta(i) ));
    }
  }
  log_aux -= (0.5/s2_mu) * (pow(log(mu1),2) - pow(log(mu0),2));

  /* CREATING OUTPUT VARIABLE & DEBUG 
   * Proposed values are automatically rejected in the following cases:
   * - If smaller than 1e-3
   * - If the proposed value is not finite
   * - When the acceptance rate cannot be numerally computed
   */
  for (int i=0; i < q0; i++)
  {
    if( arma::is_finite(log_aux(i)) & arma::is_finite(mu1(i)) )
    {
      if((log(u(i)) < log_aux(i)) & (mu1(i) > 1e-3)) { ind(i) = 1; }
      else{ind(i) = 0; mu1(i) = mu0(i);}            
    }
    else
    {
      ind(i) = 0; mu1(i) = mu0(i);
      Rcpp::Rcout << "Error when updating mu " << i+1 << std::endl;
      Rcpp::warning("Consider additional data filter if error is frequent.");        
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
 */
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
  arma::vec & ind)
{

  /* PROPOSAL STEP */
  delta1 = exp(arma::randn(q0) % sqrt(prop_var) + log(delta0));
  u = arma::randu(q0);
        
  /* ACCEPT/REJECT STEP 
   * Note: there is a -1 factor coming from the log-normal prior. 
   * However, it cancels out as using log-normal proposals.
   */
  arma::vec log_aux = - n * (lgamma_cpp(1/delta1) - lgamma_cpp(1/delta0));
  
  for (int i=0; i < q0; i++)
  {
    for (int j=0; j < n; j++) 
    {
      log_aux(i) += R::lgammafn(Counts(i,j) + (1/delta1(i)));
      log_aux(i) -= R::lgammafn(Counts(i,j) + (1/delta0(i)));
      log_aux(i) -= ( Counts(i,j)+(1/delta1(i)) ) * log( phinu(j)*mu(i)+(1/delta1(i)) );
      log_aux(i) += ( Counts(i,j)+(1/delta0(i)) ) * log( phinu(j)*mu(i)+(1/delta0(i)) );
    }
  }
  log_aux -= n * ( (log(delta1)/delta1) - (log(delta0)/delta0) );
  // Component related to the prior
  if(prior_delta == 1) 
  {
    log_aux += (log(delta1)-log(delta0))*a_delta - b_delta * (delta1 - delta0);
  }
  else 
  { 
    log_aux -= (0.5/s2delta) * (pow(log(delta1),2) - pow(log(delta0),2)); 
  }

  /* CREATING OUTPUT VARIABLE & DEBUG 
   * Proposed values are automatically rejected in the following cases:
   * - If smaller than 1e-3
   * - If the proposed value is not finite
   * - When the acceptance rate cannot be numerally computed
   */    
  for (int i=0; i < q0; i++)
  {
    if( arma::is_finite(log_aux(i)) & arma::is_finite(delta1(i)) )
    {
      if((log(u(i)) < log_aux(i)) & (delta1(i) > 1e-3)) { ind(i) = 1; }
      else {ind(i) = 0; delta1(i) = delta0(i); }
    }      
    else
    {
      ind(i) = 0; delta1(i) = delta0(i);
      Rcpp::Rcout << "Error when updating delta " << i+1 << std::endl;
      Rcpp::warning("Consider additional data filter if error is frequent.");        
    }
  }
  
  // OUTPUT
  return join_rows(delta1, ind);
}

/* Draws for cell-specific normalising constants s[j]
 * Metropolis-Hastings updates are not required as full conditionals 
 * have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells.
 */
arma::vec sUpdate(
  arma::vec const& s0, 
  arma::vec const& nuovertheta, 
  double const& theta, 
  double const& as, 
  double const& bs, 
  int const& n,
  arma::vec & s1)
{   
  // Calculating parameters to the passed as input to the Rgig function 
  // (common for all cells)
  double p = as - 1/theta; 
  double b = 2*bs;
   
  /* GIG draws
   * DEBUG: return original value of s0 if input values are not 
   * within the appropriate range 
   */ 
  if(!R_IsNA(p))  
  {     
    for (int j=0; j<n; j++) 
    {
      if(nuovertheta(j) > 0) 
      {
        s1(j) = RgigDouble(p, 2 * nuovertheta(j), b);
      }
      else
      {
        s1(j) = s0(j);
        Rcpp::Rcout << "Error when updating s " << j+1 << std::endl;
        Rcpp::warning("Consider additional data filter if error is frequent.");          
      }               
    }
    return s1;     
  }
  else return s0;
}

/* Metropolis-Hastings updates of nu 
 * Updates are implemented simulateaneously for all cells.
 */
arma::mat nuUpdate(
  arma::vec const& nu0, 
  arma::vec const& prop_var, 
  arma::mat const& Counts, 
  double const& SumSpikeInput,
  arma::vec const& mu, 
  arma::vec const& invdelta, 
  arma::vec const& phi, 
  arma::vec const& invthetas, 
  double const& theta, 
  arma::vec const& sum_bygene_all, 
  int const& q0,
  int const& n,
  arma::vec & nu1,
  arma::vec & u, 
  arma::vec & ind)
{
  // PROPOSAL STEP    
  nu1 = exp(arma::randn(n) % sqrt(prop_var) + log(nu0));
  u = arma::randu(n);
    
  /* ACCEPT/REJECT STEP 
   * Note: there is a -1 factor coming from the log-normal prior. 
   * However, it cancels out as using log-normal proposals.
   */
  arma::vec log_aux = arma::zeros(n);

  // Loop to replace matrix operations, through genes and cells
  for (int j=0; j < n; j++) 
  {
    for (int i=0; i < q0; i++) 
    {
      log_aux(j) -= ( Counts(i,j) + invdelta(i) ) *  
                      log( ( phi(j)*nu1(j)*mu(i) + invdelta(i) ) / 
                           ( phi(j)*nu0(j)*mu(i) + invdelta(i) ));
    } 
  }
  log_aux += ( log(nu1) - log(nu0) ) % (sum_bygene_all + 1/theta);
  log_aux -= ( nu1 - nu0 )  % ( SumSpikeInput + invthetas );

  /* CREATING OUTPUT VARIABLE & DEBUG 
   * Proposed values are automatically rejected in the following cases:
   * - If smaller than 1e-5
   * - If the proposed value is not finite
   * - When the acceptance rate cannot be numerally computed
   */  
  for (int j=0; j < n; j++)
  {
    if(arma::is_finite(log_aux(j)) & arma::is_finite(nu1(j)))
    {
      if( (log(u(j)) < log_aux(j)) & (nu1(j) > 1e-5) ) { ind(j) = 1; }
      else {ind(j) = 0; nu1(j) = nu0(j); }
    }      
    else
    {
      ind(j) = 0; nu1(j) = nu0(j);
      Rcpp::Rcout << "Error when updating nu " << j+1 << std::endl;
      Rcpp::warning("Consider additional data filter if error is frequent.");    
    }
  }
   
  // OUTPUT
  return join_rows(nu1, ind);
}

/* Metropolis-Hastings updates of theta 
 */
arma::vec thetaUpdate(
  double const& theta0, 
  double const& prop_var, 
  arma::vec const& s, 
  arma::vec const& nu, 
  double const& a_theta, 
  double const& b_theta, 
  int const& n)
{
  // CREATING VARIABLES WHERE TO STORE DRAWS
  double logtheta = log(theta0);

  // PROPOSAL STEP
  double y = as<double>(wrap(arma::randn(1))) * sqrt(prop_var) + logtheta;
  double u = as<double>(wrap(arma::randu(1)));
    
  // ACCEPT/REJECT STEP
  double log_aux = (y-logtheta)*a_theta;
  log_aux -= n * (logtheta/theta0) * ((y/logtheta)*exp(-y+logtheta)-1);  
  log_aux -= n * (R::lgammafn(exp(-y)) - R::lgammafn(1/theta0));
  log_aux += ((exp(-y+logtheta)-1)/theta0) * sum(log(nu/s)-(nu/s));
  log_aux -= b_theta*theta0*(exp(y-logtheta)-1);
  double aux = exp(y);
  bool ind = log(u) < log_aux;
  // DEBUG: Reject proposed values below 0.0001 (to avoid numerical innacuracies)
  if(aux<0.0001) {ind=0;}
    
  // CREATING OUTPUT VARIABLE
  double theta = ind * aux + (1 - ind) * theta0;
    
  // OUTPUT
  arma::vec output(2); output(0)=theta; output(1)=ind;
  return output;
}

/*Dirichlet sampler*/
arma::vec rDirichlet(
  arma::vec alpha)
{
  arma::vec aux = arma::ones(alpha.size());
  unsigned int i;
  for(i=0; i<alpha.size(); i++)
  {
    aux(i) = Rcpp::as<double>(rgamma(1,alpha(i),1));   
  }
  aux = aux / sum(aux);    
  return aux;  
}

/* Metropolis-Hastings updates of phi 
 * Joint updates using Dirichlet proposals
 */
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
  arma::vec & phi1) 
{
  int ind;
  
  // PROPOSAL STEP
  phi1 = n * rDirichlet(prop_var * phi0); 
  double u = R::runif(0,1);
  
  // ACCEPT/REJECT STEP (REJECT VALUES OUTSIDE VALID RANGE)  
  if(all(prop_var * phi1 < 2.5327372760800758e+305)  & 
     all(prop_var * phi0 < 2.5327372760800758e+305) &
     all(phi1 > 0) & all(phi0 > 0)) 
  {
    // There is an extra -1 but it cancels out with the proposal component
    double log_aux = sum( (sum_bygene_bio + aphi) % (log(phi1) - log(phi0)));
        
    // Loop to replace matrix operations, through genes and cells
    // There is an extra factor in the prior n^(-n); it cancels out in the ratio
    // There is an extra factor n^(-(sum(aphi) - 1));it cancels out in the ratio
    for (int j=0; j < n; j++) 
    {
      for (int i=0; i < q0; i++) 
      {
        log_aux -= ( Counts(i,j) + invdelta(i) ) *  
                      log( (phi1(j)*nu(j)*mu(i) + invdelta(i) ) / 
                           (phi0(j)*nu(j)*mu(i) + invdelta(i) ));
      } 
    }
    // There is an extra factor 
    // n^(-(sum(prop_var * phi1))) / n^(-(sum(prop_var * phi0)));
    // it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
    // There is an extra factor 
    // gamma(sum(prop_var * phi0)) / gamma(sum(prop_var * phi1));
    // it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
    log_aux += prop_var * sum(phi1 % log(phi0) - phi0 % log(phi1));
    log_aux -= sum(lgamma_cpp_vec(prop_var * phi1) - lgamma_cpp_vec(prop_var * phi0));    
    
    if(!R_IsNA(log_aux))
    {
        if(log(u) < log_aux) { ind = 1; }
        else {ind = 0; phi1 = phi0;}
    }
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
    else
    {
      Rcpp::Rcout << "Error when updating phi" << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset."); 
      ind = 0; phi1 = phi0;
    }     
  }
  else
  {
    ind = 0; phi1 = phi0;     
  }      
  return(Rcpp::List::create(
         Rcpp::Named("phi") = phi1,
         Rcpp::Named("ind") = ind)); 
}




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
  double theta0, 
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
  double LStheta0, 
  NumericVector sumByCellAll, 
  NumericVector sumByCellBio,
  NumericVector sumByGeneAll, 
  NumericVector sumByGeneBio,
  int StoreAdapt, 
  int EndAdapt,
  int PrintProgress) 
{
  using arma::ones;
  using arma::zeros;
  using Rcpp::Rcout; 
  
  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); 
  int nBatch = BatchDesign.ncol();
  int q0 = delta0.size();
  int Naux = N/Thin - Burn/Thin;
  
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); 
  arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); 
  arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts);
  arma::mat BatchDesign_arma = as_arma(BatchDesign);
  
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = zeros(Naux, q0); 
  arma::mat delta = zeros(Naux, q0); 
  arma::mat phi = ones(Naux, n);
  arma::mat s = zeros(Naux, n);  
  arma::mat nu = zeros(Naux, n); 
  arma::vec theta = zeros(Naux); 
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::vec LSphi;
  arma::mat LSnu;
  arma::vec LStheta; 
  
  // LOG-PROPOSAL VARIANCES 
  if(StoreAdapt == 1)
  {
    LSmu = zeros(Naux,q0); 
    LSdelta = zeros(Naux,q0); 
    LSphi = ones(Naux);   
    LSnu = zeros(Naux,n); 
    LStheta = zeros(Naux);   
  }
  
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = zeros(q0); arma::vec PmuAux = zeros(q0);
  arma::vec deltaAccept = zeros(q0); arma::vec PdeltaAux = zeros(q0);
  double phiAccept = 0; double PphiAux = 0;
  arma::vec nuAccept = zeros(n); arma::vec PnuAux = zeros(n);
  double thetaAccept = 0; double PthetaAux = 0;
  
  // INITIALIZATION OF PARAMETER VALUES FOR MCMC RUN
  arma::mat muAux = zeros(q0,2); muAux.col(0)=as_arma(mu0); 
  arma::mat deltaAux = zeros(q0,2); deltaAux.col(0)=as_arma(delta0); 
  arma::vec phiAux = as_arma(phi0); Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0); 
  arma::mat nuAux = zeros(n,2); nuAux.col(0) = as_arma(nu0);
  arma::vec thetaAux = zeros(2); thetaAux(0) = theta0; 
  
  // INITIALIZATION OF ADAPTIVE VARIANCES
  arma::vec LSmuAux = as_arma(LSmu0);
  arma::vec LSdeltaAux = as_arma(LSdelta0); 
  double LSphiAux = LSphi0; 
  arma::vec LSnuAux = as_arma(LSnu0);
  double LSthetaAux = LStheta0;  
  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = zeros(q0); arma::vec PdeltaAux0 = zeros(q0);
  double PphiAux0 = 0; 
  arma::vec PnuAux0 = zeros(n); double PthetaAux0 = 0;
  
  double SumSpikeInput = sum(muSpikes);
  
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES 
  // (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; 
  
  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  // To avoid repeated initialisation
  arma::vec y_q0 = ones(q0); arma::vec y_n = ones(n); 
  arma::vec ind_q0 = zeros(q0); arma::vec ind_n = zeros(n);
  arma::vec u_q0 = zeros(q0); arma::vec u_n = zeros(n);
  
  Rcout << "-------------------------------------------------------------" << std::endl;  
  Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  
  // START OF MCMC LOOP
  for (int i=0; i<N; i++) {
    
    Rcpp::checkUserInterrupt();
    
    if(i==Burn)
    {
      Rcout << "-------------------------------------------------------------" << std::endl; 
      Rcout << "End of Burn-in period."<< std::endl;
      Rcout << "-------------------------------------------------------------" << std::endl; 
    }
    
    Ibatch++; 
    
    // UPDATE OF PHI: 
    // 1st ELEMENT IS THE UPDATE, 
    // 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, 
                           muAux.col(0), 1/deltaAux.col(0),
                           nuAux.col(0), aphi, sumByGeneBio_arma, q0, n,
                           y_n); 
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]);
    if(i>=Burn) { phiAccept += Rcpp::as<double>(phiAuxList["ind"]); }

    // UPDATE OF THETA: 
    // 1st ELEMENT IS THE UPDATE, 
    // 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdate(thetaAux(0), exp(LSthetaAux), 
                           sAux, nuAux.col(0), atheta, btheta, n);
    PthetaAux += thetaAux(1); if(i>=Burn) {thetaAccept += thetaAux(1);}
    
    // UPDATE OF MU: 
    // 1st COLUMN IS THE UPDATE, 
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR       
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, 
                     1/deltaAux.col(0), phiAux % nuAux.col(0), 
                     sumByCellBio_arma, s2mu, q0, n,
                     y_q0, u_q0, ind_q0);     
    PmuAux += muAux.col(1); if(i>=Burn) {muAccept += muAux.col(1);}

    // UPDATE OF S
    sAux = sUpdate(sAux, nuAux.col(0)/thetaAux(0), thetaAux(0), 
                   as, bs, n, y_n); 

    // UPDATE OF DELTA: 
    // 1st COLUMN IS THE UPDATE, 
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, 
                           muAux.col(0), phiAux % nuAux.col(0), 
                           adelta, bdelta, s2delta, prior_delta, 
                           q0, n, y_q0, u_q0, ind_q0);  
    PdeltaAux += deltaAux.col(1); if(i>=Burn) {deltaAccept += deltaAux.col(1);}    
    
    // UPDATE OF NU: 
    // 1st COLUMN IS THE UPDATE, 
    // 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdate(nuAux.col(0), exp(LSnuAux), Counts_arma, SumSpikeInput,
                     muAux.col(0), 1/deltaAux.col(0), phiAux, 
                     1/(thetaAux(0) * sAux), thetaAux(0), 
                     sumByGeneAll_arma, q0, n, y_n, u_n, ind_n); 
    PnuAux += nuAux.col(1); if(i>=Burn) {nuAccept += nuAux.col(1);}

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux = PmuAux/50;
        PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux = LSmuAux + PmuAux*0.1;
        PdeltaAux=PdeltaAux/50; 
        PdeltaAux = -1 + 2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux = LSdeltaAux + PdeltaAux*0.1;
        PphiAux = PphiAux/50; PphiAux = -1 + 2*(PphiAux>ar); 
        LSphiAux = LSphiAux - PphiAux*0.1;  
        PnuAux = PnuAux/50; 
        PnuAux = -1 + 2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux = LSnuAux + PnuAux*0.1;
        PthetaAux = PthetaAux/50; PthetaAux = -1+2*(PthetaAux>ar); 
        LSthetaAux = LSthetaAux + PthetaAux*0.1;
      
        Ibatch = 0; 
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0; 
        PphiAux = PphiAux0; 
        PnuAux = PnuAux0; PthetaAux = PthetaAux0; 
      }
      
    }
        
    // STORAGE OF DRAWS
    if((i%Thin==0) & (i>=Burn))
    {      
      mu.row(i/Thin - Burn/Thin) = muAux.col(0).t(); 
      delta.row(i/Thin - Burn/Thin) = deltaAux.col(0).t(); 
      phi.row(i/Thin - Burn/Thin) = phiAux.t();
      s.row(i/Thin - Burn/Thin) = sAux.t();
      nu.row(i/Thin - Burn/Thin) = nuAux.col(0).t();       
      theta(i/Thin - Burn/Thin) = thetaAux(0);       
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/Thin - Burn/Thin) = LSmuAux.t();
        LSdelta.row(i/Thin - Burn/Thin) = LSdeltaAux.t();
        LSphi(i/Thin - Burn/Thin) = LSphiAux;
        LSnu.row(i/Thin - Burn/Thin) = LSnuAux.t();
        LStheta(i/Thin - Burn/Thin) = LSthetaAux; 
      }
    }
    
    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if((i%(2*Thin) == 0) & (PrintProgress == 1))
    {
        Rcout << "-------------------------------------------------------------" << std::endl;
        Rcout << "MCMC iteration " << i << " out of " << N << " has been completed." << std::endl;
        Rcout << "-------------------------------------------------------------" << std::endl;
        Rcout << "Current draws of some selected parameters are displayed below." << std::endl;
        Rcout << "mu (gene 1): " << muAux(0,0) << std::endl; 
        Rcout << "delta (gene 1): " << deltaAux(0,0) << std::endl; 
        Rcout << "phi (cell 1): " << phiAux(0) << std::endl;
        Rcout << "s (cell 1): " << sAux(0) << std::endl;
        Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
        Rcout << "theta: " << thetaAux(0) << std::endl;
        Rcout << "-------------------------------------------------------------" << std::endl;
        Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
        Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
        Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl; 
        Rcout << "LSphi: " << LSphiAux << std::endl;
        Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
        Rcout << "LStheta: " << LSthetaAux << std::endl;
    }    
  }
  
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcout << " " << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << "All " << N << " MCMC iterations have been completed." << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << "Please see below a summary of the overall acceptance rates." << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
  
  Rcout << "Minimum acceptance rate among mu[i]'s: " << 
    min(muAccept/(N-Burn)) << std::endl;
  Rcout << "Average acceptance rate among mu[i]'s: " << 
    mean(muAccept/(N-Burn)) << std::endl;
  Rcout << "Maximum acceptance rate among mu[i]'s: " << 
    max(muAccept/(N-Burn)) << std::endl;
  Rcout << " " << std::endl;
  Rcout << "Minimum acceptance rate among delta[i]'s: " << 
    min(deltaAccept/(N-Burn)) << std::endl;
  Rcout << "Average acceptance rate among delta[i]'s: " << 
    mean(deltaAccept/(N-Burn)) << std::endl;
  Rcout << "Maximum acceptance rate among delta[i]'s: " << 
    max(deltaAccept/(N-Burn)) << std::endl;
  Rcout << " " << std::endl;
  Rcout << "Acceptance rate for phi (joint): " << 
    phiAccept/(N-Burn) << std::endl;
  Rcout << " " << std::endl;
  Rcout << "Minimum acceptance rate among nu[j]'s: " << 
    min(nuAccept/(N-Burn)) << std::endl;
  Rcout << "Average acceptance rate among nu[j]'s: " << 
    mean(nuAccept/(N-Burn)) << std::endl;
  Rcout << "Maximum acceptance rate among nu[j]'s: " << 
    max(nuAccept/(N-Burn)) << std::endl;
  Rcout << " " << std::endl;
  Rcout << "Acceptance rate for theta: " << thetaAccept/(N-Burn) << std::endl;
  Rcout << "-------------------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;

      
  if(StoreAdapt == 1)
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu")=mu,
           Rcpp::Named("delta")=delta,
           Rcpp::Named("phi")=phi,
           Rcpp::Named("s")=s,
           Rcpp::Named("nu")=nu,
           Rcpp::Named("theta")=theta,
           Rcpp::Named("ls.mu")=LSmu,
           Rcpp::Named("ls.delta")=LSdelta,
           Rcpp::Named("ls.phi")=LSphi,
           Rcpp::Named("ls.nu")=LSnu,
           Rcpp::Named("ls.theta")=LStheta)); 
  }
  
  else
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu")=mu,
           Rcpp::Named("delta")=delta,
           Rcpp::Named("phi")=phi,
           Rcpp::Named("s")=s,
           Rcpp::Named("nu")=nu,
           Rcpp::Named("theta")=theta)); 
  }

}

/* Draws for cell-specific normalising constants s[j]
 * Metropolis-Hastings updates are not required as full conditionals have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::vec sUpdateBatch(
    arma::vec const& s0_arma, /* Current value of $s=(s_1,...,s_n)$' */
arma::vec const& nu_arma, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
arma::vec const& theta, /* Current value of $\theta$ */
double const& as, /* Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
double const& bs, /* Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
arma::mat const& BatchDesign,
int const& n)
{
  arma::vec thetaBatch = BatchDesign * theta; 
  
  // Creating a vector where draws will be stored
  arma::vec aux = arma::zeros(n);
  
  // Calculating parameters to the passed as input to the Rgig function (common for all cells)
  arma::vec p = as - 1 / thetaBatch; 
  double b = 2 * bs;
  
  // GIG draws
  /* DEBUG: return original value of s0 if input values are not wiThin the appropriate range (to avoid problems related to numerical innacuracy) */
  
  // Calculating parameter to the passed as input to the Rgig function (specific to each cell)
  //   NumericVector a = 2 * nu / thetaBatch;     
  arma::vec a = 2 * nu_arma / thetaBatch;
  int j;
  for (j=0; j<n; j++) 
  {
    if(!R_IsNA(p(j))) 
    {
      if(!R_IsNA(a(j)) & (a(j)>0)) 
      {
        aux(j) = Rcpp::as<double>(Rgig(1,p(j),a(j),b));
        /* DEBUG: break in case of undefined values */
        if(R_IsNA(aux(j))) 
        {
          Rcpp::Rcout << "SomeThing went wrong when updating s" << j << std::endl;
          Rcpp::stop("Please consider additional filter of the input dataset.");
        }
      }
      else if(!(a(j)<0) & (p(j)>0)) {aux(j) = Rcpp::as<double>(Rgig(1,p(j),a(j),b));}
      else{aux(j)=s0_arma(j);}         
    }
    else{aux(j)=s0_arma(j);}
  }
  return aux;     
}


/* Auxiliary function required for some of the Metropolis-Hastings updates of nu */
arma::mat UpdateAux_nuTrick(   
    arma::vec const& nu, /* Auxiliary variable (function of the current and proposed value of nu) */
arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
arma::vec const& phi) /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
{
  arma::mat x = ((phi % nu) * mu.t()).t();
  x.each_col() += 1/delta;
  return x;
}

/* Metropolis-Hastings updates of nu 
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::mat nuUpdateBatch(
  arma::vec const& nu0, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\nu=(\nu_1,...,\nu_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  double const& SumSpikeInput,
  arma::mat const& BatchDesign, 
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */  
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& theta, /* Current value of $\theta$' */
  arma::vec const& sum_bygene_all, /* Sum of expression counts by gene (all genes) */
  int const& q0,
  int const& n)
{
    using arma::span;
    
    // CREATING VARIABLES WHERE TO STORE DRAWS
    arma::vec lognu = log(nu0);
    arma::vec thetaBatch = BatchDesign * theta; 

    // PROPOSAL STEP    
    arma::vec y = arma::randn(n) % sqrt(prop_var) + lognu;
    arma::vec u = arma::randu(n);
       
    // ACCEPT/REJECT STEP
    
// CODE CRUSHES WHEN USING THIS LOOP, CHECK WHY
//    arma::vec log_aux = arma::zeros(n);

    // Loop to replace matrix operations, through genes and cells
//    for (int j=0; j < n; j++) 
//    {
//      for (int i=0; i < q0; i++) 
//      {
//        log_aux(j) -= ( Counts(i,j) + (1/delta(i)) ) *  
//                      log( ( phi(j)*exp(y(j))*mu(i)+(1/delta(i)) ) / ( phi(j)*nu0(j)*mu(i)+(1/delta(i)) ));
//      } 
//      log_aux(j) += (y(j) - lognu(j)) * (sum_bygene_all(j) + 1/thetaBatch(j)); // - ( y(j)-nu0(j) )  * (SumSpikeInput + (1/(thetaBatch(j)*s(j))));
//      log_aux(j) -= nu0(j) * (exp(y(j)-lognu(j))-1)  * (SumSpikeInput + (1/(thetaBatch(j) * s(j))));
//    }

    
    arma::vec log_aux = (y - lognu) % (sum_bygene_all + 1/thetaBatch);
    log_aux -= nu0 % (exp(y-lognu)-1)  % (SumSpikeInput + (1/(thetaBatch % s)));   
    arma::mat m = Counts.rows(0, q0 - 1);
    m.each_col() += 1 / delta;   
    arma::mat num = UpdateAux_nuTrick(exp(y), mu(span(0, q0 - 1)), delta,phi);
    num /= UpdateAux_nuTrick(nu0,    mu(span(0, q0 - 1)), delta,phi);
    m %= log(num);
    log_aux -= sum(m, 0).t();
    arma::umat ind = log(u) < log_aux;
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message    
    ind.elem(find_nonfinite(log_aux)).zeros();
    if(size(find_nonfinite(log_aux),0)>0)
    {
      Rcpp::Rcout << "SomeThing went wrong when updating nu" << size(find_nonfinite(log_aux),0) << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset.");
    }

    // CREATING OUTPUT VARIABLE        
    arma::vec nu = ind % exp(y) + (1 - ind) % nu0;
   
    // OUTPUT
    return join_rows(nu, arma::conv_to<arma::mat>::from(ind));
}

/* Metropolis-Hastings updates of theta 
 */
arma::mat thetaUpdateBatch(
  arma::vec const& theta0, /* Current value of $\theta$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\theta$ */
  arma::mat const& BatchDesign, 
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& a_theta, /* Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  double const& b_theta, /* Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  int const& n,
  int const& nBatch)
  {
    using arma::span;

    // CREATING VARIABLES WHERE TO STORE DRAWS
    arma::vec logtheta = log(theta0);

    // PROPOSAL STEP
    arma::vec y = arma::randn(nBatch) % sqrt(prop_var) + logtheta;
    arma::vec u = arma::randu(nBatch);
    
    arma::mat BatchDesignAux = BatchDesign.cols(0, nBatch - 1);
    arma::vec nBatches = sum(BatchDesignAux,0).t();
    
    BatchDesignAux.each_col() %= log(nu / s) - (nu / s);
    
    // ACCEPT/REJECT STEP
    arma::vec log_aux = (y-logtheta) * a_theta;
    log_aux -= nBatches % (logtheta/theta0) % ((y/logtheta) % exp(-y+logtheta)-1) + nBatches % (lgamma_cpp(exp(-y))-lgamma_cpp(1/theta0));
    log_aux += ((exp(-y+logtheta)-1)/theta0) % sum(BatchDesignAux,0).t() - b_theta * theta0 % (exp(y-logtheta)-1);
    arma::umat ind = log(u) < log_aux;
    // DEBUG: Reject proposed values below 0.0001 (to avoid numerical innacuracies)
    ind %= 0.0001 < exp(y);

    // CREATING OUTPUT VARIABLE
    arma::vec theta = ind % exp(y) + (1 - ind) % theta0;
    
    // OUTPUT
    return join_rows(theta, arma::conv_to<arma::mat>::from(ind));
}

/* MCMC sampler 
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcppBatch(
  int N, // Total number of MCMC draws 
  int Thin, // Thinning period for MCMC chain 
  int Burn, // Burning period for MCMC chain 
  NumericMatrix Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom) 
  NumericMatrix BatchDesign, // Design matrix representing batch information (number of columns must be equal to number of batches)
  NumericVector muSpikes, 
  NumericVector mu0, // Starting value of $\mu=(\mu_1,...,\mu_q)'$ (true mRNA content for technical genes)  
  NumericVector delta0, // Starting value of $\delta=(\delta_1,...,\delta_{q_0})'$  
  NumericVector phi0, // Starting value of $\phi=(\phi_1,...,\phi_n)$'$ 
  NumericVector s0, // Starting value of $s=(s_1,...,s_n)$'$ 
  NumericVector nu0, // Starting value of $\nu=(\nu_1,...,\nu_n)$'$   
  double theta0, // Starting value of $\theta$ 
  double s2mu, 
  double adelta, // Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  double bdelta, // Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  NumericVector aphi, // Dirichlet hyper-parameter for $\phi / n$ 
  double as, // Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double bs, // Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$  
  double atheta, // Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
  double btheta, // Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
  double ar, // Optimal acceptance rate for adaptive Metropolis-Hastings updates
  NumericVector LSmu0, // Starting value of adaptive proposal variance of $\mu=(\mu_1,...,\mu_q)'$ (log-scale)
  NumericVector LSdelta0, // Starting value of adaptive proposal variance of $\delta=(\delta_1,...,\delta_{q_0})'$ (log-scale)
  double LSphi0, // Starting value of adaptive proposal precision of $\phi=(\phi_1,...,\phi_n)'$ (log-scale)
  NumericVector LSnu0, // Starting value of adaptive proposal variance of $\nu=(\nu_1,...,\nu_n)'$ (log-scale)
  double LStheta0, // Starting value of adaptive proposal variance of $\theta$ (log-scale)  
  NumericVector sumByCellAll, // Sum of expression counts by cell (all genes)
  NumericVector sumByCellBio, // Sum of expression counts by cell (biological genes only)
  NumericVector sumByGeneAll, // Sum of expression counts by gene (all genes)
  NumericVector sumByGeneBio,
  int StoreAdapt, 
  int EndAdapt,
  int PrintProgress,
  double s2delta, 
  double prior_delta) 
{
  using arma::ones;
  using arma::zeros;
  using Rcpp::Rcout; 

  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); 
  int q0 = delta0.size(); 
  int Naux = N/Thin - Burn/Thin;
  int nBatch = BatchDesign.ncol();
  
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); 
  arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); 
  arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts); 
  arma::mat BatchDesign_arma = as_arma(BatchDesign);
  
  double SumSpikeInput = sum(muSpikes);
  
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = arma::zeros(Naux, q0); 
  arma::mat delta = arma::zeros(Naux, q0); 
  arma::mat phi = arma::ones(Naux, n);
  arma::mat s = arma::zeros(Naux, n);  
  arma::mat nu = arma::zeros(Naux, n); 
  arma::mat theta = arma::zeros(Naux, nBatch); 
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::vec LSphi;
  arma::mat LSnu;
  arma::mat LStheta; 

  // LOG-PROPOSAL VARIANCES 
  if(StoreAdapt == 1)
  {
    LSmu = arma::zeros(Naux, q0); 
    LSdelta = arma::zeros(Naux, q0); 
    LSphi = arma::ones(Naux);   
    LSnu = arma::zeros(Naux, n); 
    LStheta = arma::zeros(Naux, nBatch);   
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(q0); arma::vec PmuAux = arma::zeros(q0);
  arma::vec deltaAccept = arma::zeros(q0); arma::vec PdeltaAux = arma::zeros(q0);
  double phiAccept = 0; double PphiAux = 0;
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  arma::vec thetaAccept = arma::zeros(nBatch); arma::vec PthetaAux = arma::zeros(nBatch);
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(q0,2); muAux.col(0)=as_arma(mu0); arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(q0,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::vec phiAux = as_arma(phi0); double LSphiAux = LSphi0; Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0); 
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::mat thetaAux = arma::zeros(nBatch, 2); thetaAux.col(0) = theta0 * arma::ones(nBatch); 
  arma::vec LSthetaAux = LStheta0 * arma::ones(nBatch);  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q0); arma::vec PdeltaAux0 = arma::zeros(q0);
  double PphiAux0 = 0; 
  arma::vec PnuAux0 = arma::zeros(n); arma::vec PthetaAux0 = arma::ones(nBatch);
  
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;
  
  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  // To avoid repeated initialisation
  arma::vec y_q0 = ones(q0); arma::vec y_n = ones(n); 
  arma::vec ind_q0 = zeros(q0); arma::vec ind_n = zeros(n);
  arma::vec u_q0 = zeros(q0); arma::vec u_n = zeros(n);
  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;  
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  
  // START OF MCMC LOOP
  for (i=0; i<N; i++) {
    
    Rcpp::checkUserInterrupt();
    
    if(i==Burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
      Rcpp::Rcout << "End of Burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
    }
    
    Ibatch++; 
    
//    struct timespec time0_1 = orwl_gettime();
    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, muAux.col(0), deltaAux.col(0),
                           nuAux.col(0), aphi, sumByGeneBio_arma, q0,n, y_n); 
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]); if(i>=Burn) {phiAccept += Rcpp::as<double>(phiAuxList["ind"]);}
//    struct timespec time1_1 = orwl_gettime();    
//    Rcpp::Rcout << "Time phi: " << (time1_1.tv_nsec - time0_1.tv_nsec) / ((float)(n)) << std::endl;

//    struct timespec time0_2 = orwl_gettime();
    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdateBatch(thetaAux.col(0), exp(LSthetaAux), BatchDesign_arma,
                                sAux, nuAux.col(0), atheta, btheta, n, nBatch);
    PthetaAux += thetaAux.col(1); if(i>=Burn) {thetaAccept += thetaAux.col(1);}
//    struct timespec time1_2 = orwl_gettime();    
//    Rcpp::Rcout << "Time theta: "  <<  (time1_2.tv_nsec - time0_2.tv_nsec) / ((float)(nBatch)) << std::endl;
    
//    struct timespec time0_3 = orwl_gettime();
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR 
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, 
                     1/deltaAux.col(0), phiAux % nuAux.col(0), 
                     sumByCellBio_arma, s2mu, q0, n,
                     y_q0, u_q0, ind_q0);     
    PmuAux += muAux.col(1); if(i>=Burn) {muAccept += muAux.col(1);}

    
//    struct timespec time0_4 = orwl_gettime();
    // UPDATE OF S
    sAux = sUpdateBatch(sAux, nuAux.col(0), thetaAux.col(0), as, bs, BatchDesign_arma, n);  
//    struct timespec time1_4 = orwl_gettime();    
//    Rcpp::Rcout << "Time s: "  <<  (time1_4.tv_nsec - time0_4.tv_nsec) / ((float)(n)) << std::endl;
    

    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, 
                           muAux.col(0), phiAux % nuAux.col(0), 
                           adelta, bdelta, s2delta, prior_delta, 
                           q0, n, y_q0, u_q0, ind_q0);  
    PdeltaAux += deltaAux.col(1); if(i>=Burn) {deltaAccept += deltaAux.col(1);}   

    
//    struct timespec time0_6 = orwl_gettime();
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateBatch(nuAux.col(0), exp(LSnuAux), Counts_arma, SumSpikeInput,
                     BatchDesign_arma,
                     muAux.col(0), deltaAux.col(0),
                     phiAux, sAux, thetaAux.col(0), sumByGeneAll_arma, q0, n); 
    PnuAux += nuAux.col(1); if(i>=Burn) {nuAccept += nuAux.col(1);}
//    struct timespec time1_6 = orwl_gettime();    
//    Rcpp::Rcout << "Time nu: " <<  (time1_6.tv_nsec - time0_6.tv_nsec) / ((float)(n)) << std::endl;

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/50; PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux=LSmuAux+PmuAux*0.1;
        //LSmuAux=LSmuAux+PmuAux*std::min(0.01,1/sqrt(i));
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux=LSdeltaAux+PdeltaAux*0.1;
        //LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.01,1/sqrt(i));
        PphiAux=PphiAux/50; PphiAux = -1+2*(PphiAux>ar);
        LSphiAux=LSphiAux - PphiAux*0.1;
        //LSphiAux=LSphiAux - PphiAux*std::min(0.01,1/sqrt(i));
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux=LSnuAux+PnuAux*0.1;
        //LSnuAux=LSnuAux+PnuAux*std::min(0.01,1/sqrt(i));
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*arma::conv_to<arma::mat>::from(PthetaAux>ar); 
        LSthetaAux=LSthetaAux+PthetaAux*0.1;
        //LSthetaAux= LSthetaAux + PthetaAux*std::min(0.01,1/sqrt(i));
                
        Ibatch = 0; 
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0; 
        PphiAux = PphiAux0; 
        PnuAux = PnuAux0; PthetaAux = PthetaAux0; 
      }
      
    }
        
    // STORAGE OF DRAWS
    if(i%Thin==0 & i>=Burn)
    {      
      mu.row(i/Thin - Burn/Thin) = muAux.col(0).t(); 
      delta.row(i/Thin - Burn/Thin) = deltaAux.col(0).t(); 
      phi.row(i/Thin - Burn/Thin) = phiAux.t();
      s.row(i/Thin - Burn/Thin) = sAux.t();
      nu.row(i/Thin - Burn/Thin) = nuAux.col(0).t();       
      theta.row(i/Thin - Burn/Thin) = thetaAux.col(0).t();   
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/Thin - Burn/Thin) = LSmuAux.t();
        LSdelta.row(i/Thin - Burn/Thin) = LSdeltaAux.t();
        LSphi(i/Thin - Burn/Thin) = LSphiAux;
        LSnu.row(i/Thin - Burn/Thin) = LSnuAux.t();
        LStheta.row(i/Thin - Burn/Thin) = LSthetaAux.t(); 
      }
    }
    
    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if(i%(2*Thin) == 0 & PrintProgress == 1)
    {
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "MCMC iteration " << i << " out of " << N << " has been completed." << std::endl;
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "Current draws of some selected parameters are displayed below." << std::endl;
        Rcpp::Rcout << "mu (gene 1): " << muAux(0,0) << std::endl; 
        Rcpp::Rcout << "delta (gene 1): " << deltaAux(0,0) << std::endl; 
        Rcpp::Rcout << "phi (cell 1): " << phiAux(0) << std::endl;
        Rcpp::Rcout << "s (cell 1): " << sAux(0) << std::endl;
        Rcpp::Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
        Rcpp::Rcout << "theta (batch 1): " << thetaAux(0,0) << std::endl;
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
        Rcpp::Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
        Rcpp::Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl; 
        Rcpp::Rcout << "LSphi: " << LSphiAux << std::endl;
        Rcpp::Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
        Rcpp::Rcout << "LStheta (batch 1): " << LSthetaAux(0) << std::endl;
    }    
  }
  
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "All " << N << " MCMC iterations have been completed." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Please see below a summary of the overall acceptance rates." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  
  Rcpp::Rcout << "Minimum acceptance rate among mu[i]'s: " << min(muAccept(arma::span(0, q0 - 1))/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << mean(muAccept(arma::span(0, q0 - 1))/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among mu[i]'s: " << max(muAccept(arma::span(0, q0 - 1))/(N-Burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among delta[i]'s: " << min(deltaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << mean(deltaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << max(deltaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for phi (joint): " << phiAccept/(N-Burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among nu[j]'s: " << min(nuAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among nu[j]'s: " << mean(nuAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among nu[j]'s: " << max(nuAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among theta's: " << min(thetaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among theta's: " << mean(thetaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among theta's: " << max(thetaAccept/(N-Burn)) << std::endl;  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
      
  if(StoreAdapt == 1)
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu") = mu,
           Rcpp::Named("delta") = delta,
           Rcpp::Named("phi") = phi,
           Rcpp::Named("s") = s,
           Rcpp::Named("nu") = nu,
           Rcpp::Named("theta") = theta,
           Rcpp::Named("ls.mu") = LSmu,
           Rcpp::Named("ls.delta") = LSdelta,
           Rcpp::Named("ls.phi") = LSphi,
           Rcpp::Named("ls.nu") = LSnu,
           Rcpp::Named("ls.theta") = LStheta)); 
  }
  
  else
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu") = mu,
           Rcpp::Named("delta") = delta,
           Rcpp::Named("phi") = phi,
           Rcpp::Named("s") = s,
           Rcpp::Named("nu") = nu,
           Rcpp::Named("theta") = theta)); 
  }

}

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdateNoSpikes(
    arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q_0)'$ */
    arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
    double const& Constrain,
    arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */  
    arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
    arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
    arma::vec const& sum_bycell_all, /* Sum of expression counts by cell (biological genes only) */
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu,
    arma::vec & ind,
    int const& RefGene,
    arma::uvec const& ConstrainGene,
    arma::uvec const& NotConstrainGene,
    int const& ConstrainType)
{
  using arma::span;
  
  // PROPOSAL STEP    
  arma::vec y = exp(arma::randn(q0) % sqrt(prop_var) + log(mu0));
  arma::vec u = arma::randu(q0);
  // INITIALIZE MU
  mu = mu0 + 1 - 1;
  double aux; double iAux;
  double sumAux = sum(log(mu0.elem(ConstrainGene))) - log(mu0(RefGene));
  
  // ACCEPT/REJECT STEP
  
  // Step 1: Computing the likelihood contribution of the acceptance rate 
  // Calculated in the same way for all genes, but the reference one (no need to be sequential)
  arma::vec log_aux = (log(y) - log(mu0)) % sum_bycell_all;
  for (int i=0; i < q0; i++)
  {
    if(i != RefGene)
    {
      for (int j=0; j < n; j++) 
      {
        log_aux(i) -= ( Counts(i,j) + (1/delta(i)) ) * 
          log( ( nu(j)*y(i)+(1/delta(i)) ) / ( nu(j)*mu(i)+(1/delta(i)) ));
      }
    }
  }
  
  // Step 2: Computing prior component of the acceptance rate 
  
  // Step 2.1: For genes that are under the constrain (excluding the reference one)
  for (int i=0; i < ConstrainGene.size(); i++)
  {
    iAux = ConstrainGene(i);
    if(iAux != RefGene)
    {
      aux = 0.5 * (ConstrainGene.size() * Constrain - (sumAux - log(mu(iAux))));
      log_aux(iAux) -= (0.5 * 2 /s2_mu) * (pow(log(y(iAux)) - aux,2)); 
      log_aux(iAux) += (0.5 * 2 /s2_mu) * (pow(log(mu0(iAux)) - aux,2));
      // ACCEPT REJECT
      if((log(u(iAux)) < log_aux(iAux)) & (y(iAux) > 1e-3)) 
      {
        ind(iAux) = 1; mu(iAux) = y(iAux);
        sumAux += log(mu(iAux)) - log(mu0(iAux)); 
      }
      else{ind(iAux) = 0; mu(iAux) = mu0(iAux); }      
    }
  }

  // Step 2.2: For the reference gene 
  ind(RefGene) = 1;
  mu(RefGene) = exp(ConstrainGene.size() * Constrain - sumAux);

  // Step 2.3: For genes that are *not* under the constrain
  if(ConstrainType != 1)
  {
    for (int i=0; i < NotConstrainGene.size(); i++)
    {
      iAux = NotConstrainGene(i);
      log_aux(iAux) -= (0.5/s2_mu) * (pow(log(y(iAux)),2) - pow(log(mu0(iAux)),2));
      // ACCEPT REJECT
      if((log(u(iAux)) < log_aux(iAux)) & (y(iAux) > 1e-3)) 
      {
        ind(iAux) = 1; mu(iAux) = y(iAux);
      }
      else{ind(iAux) = 0; mu(iAux) = mu0(iAux);}
    }
  }
  // OUTPUT
  return join_rows(mu, ind);
}

/* Metropolis-Hastings updates of delta
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat deltaUpdateNoSpikes(
  arma::vec const& delta0, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& prop_var,  /* Current value of the proposal variances for $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
  double const& a_delta, /* Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  double const& b_delta, /* Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  int const& q0,
  int const& n,
  double const& s2delta,
  double const& prior_delta)
{
  using arma::span;
  
  // CREATING VARIABLES WHERE TO STORE DRAWS
  arma::vec delta = arma::zeros(q0); 
  arma::vec ind = arma::zeros(q0);
  
  // PROPOSAL STEP
  arma::vec y = exp(arma::randn(q0) % sqrt(prop_var) + log(delta0));
  arma::vec u = arma::randu(q0);
  
  // ACCEPT/REJECT STEP 
  arma::vec log_aux = - n * (lgamma_cpp(1/y)-lgamma_cpp(1/delta0));
  // +1 should appear because we update log(delta) not delta. However, it cancels out with the prior. 
  log_aux -= n * ( (log(y)/y) - (log(delta0)/delta0) );

  // Loop to replace matrix operations, through genes and cells
  for (int i=0; i < q0; i++)
  {
      for (int j=0; j < n; j++) 
      {
        log_aux(i) += R::lgammafn(Counts(i,j) + (1/y(i))) - R::lgammafn(Counts(i,j) + (1/delta0(i)));
        log_aux(i) -= ( Counts(i,j) + (1/y(i)) ) *  log( nu(j)*mu(i)+(1/y(i)) );
        log_aux(i) += ( Counts(i,j) + (1/delta0(i)) ) *  log( nu(j)*mu(i)+(1/delta0(i)) );
      } 
  }
  
  // Component related to the prior
  if(prior_delta == 1) {log_aux += (log(y) - log(delta0)) * a_delta - b_delta * (y - delta0);}
  else {log_aux -= (0.5/s2delta) * (pow(log(y),2) - pow(log(delta0),2));}

  // CREATING OUTPUT VARIABLE & DEBUG
  for (int i=0; i < q0; i++)
  {
    if((arma::is_finite(log_aux(i))))
    {
      if((y(i) > 1e-3) & (log(u(i)) < log_aux(i)))
      {
        // DEBUG: Reject very small values to avoid numerical issues
        if(arma::is_finite(y(i))) {ind(i) = 1; delta(i) = y(i);}
        else{ind(i) = 0; delta(i) = delta0(i);}            
      }
      else{ind(i) = 0; delta(i) = delta0(i);}
    }      
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
    else
    {
      ind(i) = 0; delta(i) = delta0(i);
      Rcpp::Rcout << "delta0(i): " << delta0(i) << std::endl;
      Rcpp::Rcout << "y(i): " << y(i) << std::endl;
      Rcpp::Rcout << "mu(i): " << mu(i) << std::endl;
      Rcpp::Rcout << "SomeThing went wrong when updating delta " << i+1 << std::endl;
      Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");        
    }
  }
  
  // OUTPUT
  return join_rows(delta, ind);
}

/* Metropolis-Hastings updates of nu 
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::mat nuUpdateNoSpikes(
  arma::vec const& nu0, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\nu=(\nu_1,...,\nu_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::mat const& BatchDesign, 
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */  
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& theta, /* Current value of $\theta$' */
  arma::vec const& sum_bygene_all, /* Sum of expression counts by gene (all genes) */
  int const& q0,
  int const& n,
  arma::vec const& BatchInfo,
  arma::vec const& BatchIds,
  int const& nBatch,
  arma::vec const& BatchSizes,
  arma::vec const& BatchOffSet)
{
  using arma::span;
  
  // CREATING VARIABLES WHERE TO STORE DRAWS
  arma::vec lognu = log(nu0);
  arma::vec thetaBatch = BatchDesign * theta; 
  
  // PROPOSAL STEP    
  arma::vec y = arma::randn(n) % sqrt(prop_var) + lognu;
  arma::vec u = arma::randu(n);
  
  // ACCEPT/REJECT STEP
  arma::vec log_aux = (y - lognu) % (sum_bygene_all + 1/thetaBatch);
  log_aux -= nu0 % (exp(y-lognu)-1)  % (1/(thetaBatch % phi));   
  arma::mat m = Counts.rows(0, q0 - 1);
  m.each_col() += 1 / delta;   
  arma::mat num = UpdateAux_nuTrick(exp(y), mu, delta, arma::ones(n));
  num /= UpdateAux_nuTrick(nu0,    mu, delta, arma::ones(n));
  m %= log(num);
  log_aux -= sum(m, 0).t();
  arma::umat ind = log(u) < log_aux;
  // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
  // DEBUG: Print warning message    
  ind.elem(find_nonfinite(log_aux)).zeros();
  if(size(find_nonfinite(log_aux),0)>0)
  {
    Rcpp::Rcout << "SomeThing went wrong when updating nu" << size(find_nonfinite(log_aux),0) << std::endl;
    Rcpp::stop("Please consider additional filter of the input dataset.");
  }
  
  // CREATING OUTPUT VARIABLE        
  arma::vec nu = ind % exp(y) + (1 - ind) % nu0;
  
  // OUTPUT
  return join_rows(nu, arma::conv_to<arma::mat>::from(ind));
}

/* MCMC sampler 
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcppNoSpikes(
    int N, // Total number of MCMC draws 
    int Thin, // Thinning period for MCMC chain 
    int Burn, // Burning period for MCMC chain 
    NumericMatrix Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom) 
    NumericMatrix BatchDesign, // Design matrix representing batch information (number of columns must be equal to number of batches)
    NumericVector mu0, // Starting value of $\mu=(\mu_1,...,\mu_q_0)'$ (true mRNA content for technical genes)  
    NumericVector delta0, // Starting value of $\delta=(\delta_1,...,\delta_{q_0})'$  
    NumericVector phi0, // Starting value of $\phi=(\phi_1,...,\phi_n)$'$ 
    NumericVector nu0, // Starting value of $\nu=(\nu_1,...,\nu_n)$'$   
    double theta0, // Starting value of $\theta$ 
    double s2mu, 
    double adelta, // Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
    double bdelta, // Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
    double aphi,
    double bphi,
    double atheta, // Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
    double btheta, // Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
    double ar, // Optimal acceptance rate for adaptive Metropolis-Hastings updates
    NumericVector LSmu0, // Starting value of adaptive proposal variance of $\mu=(\mu_1,...,\mu_q_0)'$ (log-scale)
    NumericVector LSdelta0, // Starting value of adaptive proposal variance of $\delta=(\delta_1,...,\delta_{q_0})'$ (log-scale)
    NumericVector LSnu0, // Starting value of adaptive proposal variance of $\nu=(\nu_1,...,\nu_n)'$ (log-scale)
    double LStheta0, // Starting value of adaptive proposal variance of $\theta$ (log-scale)  
    NumericVector sumByCellAll, // Sum of expression counts by cell (all genes)
    NumericVector sumByGeneAll, // Sum of expression counts by gene (all genes)
    int StoreAdapt, 
    int EndAdapt,
    int PrintProgress,
    double s2delta, 
    double prior_delta,
    NumericVector BatchInfo,
    NumericVector BatchIds,
    NumericVector BatchSizes,
    NumericVector BatchOffSet,
    double Constrain,
    NumericVector Index,
    int RefGene,
    NumericVector RefGenes,
    NumericVector ConstrainGene,
    NumericVector NotConstrainGene,
    int ConstrainType)
{
  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int q0 = Counts.nrow(); int Naux = N/Thin - Burn/Thin;
  int nBatch = BatchDesign.ncol();
  
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll);
  arma::mat Counts_arma = as_arma(Counts); arma::mat BatchDesign_arma = as_arma(BatchDesign);
  arma::vec BatchInfo_arma = as_arma(BatchInfo);
  arma::vec BatchIds_arma = as_arma(BatchIds);
  arma::vec BatchSizes_arma = as_arma(BatchSizes);
  arma::vec BatchOffSet_arma = as_arma(BatchOffSet);
  arma::vec mu0_arma = as_arma(mu0);
  arma::vec phi0_arma = as_arma(phi0);
  arma::vec Index_arma = as_arma(Index);
  arma::uvec ConstrainGene_arma = Rcpp::as<arma::uvec>(ConstrainGene);
  arma::uvec NotConstrainGene_arma = Rcpp::as<arma::uvec>(NotConstrainGene);
  arma::vec RefGenes_arma = as_arma(RefGenes);
  
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = arma::zeros(Naux, q0); 
  arma::mat delta = arma::zeros(Naux, q0); 
  arma::mat phi = arma::ones(Naux, n);
  arma::mat nu = arma::zeros(Naux, n); 
  arma::mat theta = arma::zeros(Naux, nBatch); 
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::mat LSnu;
  arma::mat LStheta; 
  
  // LOG-PROPOSAL VARIANCES 
  if(StoreAdapt == 1)
  {
    LSmu = arma::zeros(Naux, q0); 
    LSdelta = arma::zeros(Naux, q0); 
    LSnu = arma::zeros(Naux, n); 
    LStheta = arma::zeros(Naux, nBatch);   
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(q0); arma::vec PmuAux = arma::zeros(q0);
  arma::vec deltaAccept = arma::zeros(q0); arma::vec PdeltaAux = arma::zeros(q0);
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  arma::vec thetaAccept = arma::zeros(nBatch); arma::vec PthetaAux = arma::zeros(nBatch);
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(q0,2); muAux.col(0)=as_arma(mu0); arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(q0,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::vec phiAux = phi0_arma; 
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::mat thetaAux = arma::zeros(nBatch, 2); thetaAux.col(0) = theta0 * arma::ones(nBatch); 
  arma::vec LSthetaAux = LStheta0 * arma::ones(nBatch);  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q0); arma::vec PdeltaAux0 = arma::zeros(q0);
  arma::vec PnuAux0 = arma::zeros(n); arma::vec PthetaAux0 = arma::ones(nBatch);
  
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;
  
  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  arma::vec muUpdateAux = arma::ones(q0);
  arma::vec indQ = arma::zeros(q0);
  double LSmuAuxExtra = 0;
  arma::vec RefFreq = arma::zeros(q0); 
  int RefAux;
  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;  
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  
  // START OF MCMC LOOP
  for (i=0; i<N; i++) {
    
    Rcpp::checkUserInterrupt();
    
    if(i==Burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
      Rcpp::Rcout << "End of Burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
    }
    
    Ibatch++; 
    
    // UPDATE OF PHI
    // WE CAN RECYCLE THE SAME FULL CONDITIONAL AS IMPLEMENTED FOR S (BATCH CASE)
    phiAux = sUpdateBatch(phiAux, nuAux.col(0), thetaAux.col(0),
                          aphi, bphi, BatchDesign_arma, n);
//    Rcpp::Rcout << "phi" << phiAux.t() << std::endl;
    
    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdateBatch(thetaAux.col(0), exp(LSthetaAux), BatchDesign_arma,
                                phiAux, nuAux.col(0), atheta, btheta, n, nBatch);
    PthetaAux += thetaAux.col(1); if(i>=Burn) {thetaAccept += thetaAux.col(1);}
//    Rcpp::Rcout << "theta" << thetaAux.col(0).t() << std::endl;
    
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR 
    // Additional steps required for ConstrainType = 3 (stochastic reference)
    if((ConstrainType == 3) | (ConstrainType == 5))
    {
      RefAux = as_scalar(arma::randi( 1, arma::distr_param(0, RefGenes_arma.size()-1) ));
      RefGene = RefGenes(RefAux); 
      if(i >= Burn) {RefFreq(RefGene) += 1;}
    }
    muAux = muUpdateNoSpikes(muAux.col(0), exp(LSmuAux), Constrain, Counts_arma, deltaAux.col(0), 
                             nuAux.col(0), sumByCellAll_arma, s2mu, q0, n,
                             muUpdateAux, indQ, RefGene, ConstrainGene_arma, NotConstrainGene_arma,
                             ConstrainType);
    PmuAux += muAux.col(1); if(i>=Burn) {muAccept += muAux.col(1);}  
//    Rcpp::Rcout << "mu" << muAux.col(0).t() << std::endl;
    
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdateNoSpikes(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, 
                           muAux.col(0), nuAux.col(0), adelta, bdelta, q0, n, s2delta, prior_delta);  
    PdeltaAux += deltaAux.col(1); if(i>=Burn) {deltaAccept += deltaAux.col(1);}
//    Rcpp::Rcout << "delta" << deltaAux.col(0).t() << std::endl;
   
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateNoSpikes(nuAux.col(0), exp(LSnuAux), Counts_arma, 
                            BatchDesign_arma,
                            muAux.col(0), deltaAux.col(0),
                            phiAux, thetaAux.col(0), sumByGeneAll_arma, q0, n,
                            BatchInfo_arma, BatchIds_arma, nBatch,
                            BatchSizes_arma, BatchOffSet_arma); 
    PnuAux += nuAux.col(1); if(i>=Burn) {nuAccept += nuAux.col(1);}
//    Rcpp::Rcout << "nu" << nuAux.col(0).t() << std::endl;

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/(50-RefFreq); PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);//
        LSmuAux.elem(find(Index_arma != RefGene)) = LSmuAux.elem(find(Index_arma != RefGene)) + PmuAux.elem(find(Index_arma != RefGene))*0.1; 
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux=LSdeltaAux+PdeltaAux*0.1;                
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux=LSnuAux+PnuAux*0.1; 
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*arma::conv_to<arma::mat>::from(PthetaAux>ar); 
        LSthetaAux= LSthetaAux + PthetaAux*0.1;
        
        Ibatch = 0; 
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0; 
        PnuAux = PnuAux0; PthetaAux = PthetaAux0;
      }
      
    }
    
    // STORAGE OF DRAWS
    if(i%Thin==0 & i>=Burn)
    {      
      mu.row(i/Thin - Burn/Thin) = muAux.col(0).t(); 
      delta.row(i/Thin - Burn/Thin) = deltaAux.col(0).t(); 
      phi.row(i/Thin - Burn/Thin) = phiAux.t();
      nu.row(i/Thin - Burn/Thin) = nuAux.col(0).t();       
      theta.row(i/Thin - Burn/Thin) = thetaAux.col(0).t();   
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/Thin - Burn/Thin) = LSmuAux.t();
        LSdelta.row(i/Thin - Burn/Thin) = LSdeltaAux.t();
        LSnu.row(i/Thin - Burn/Thin) = LSnuAux.t();
        LStheta.row(i/Thin - Burn/Thin) = LSthetaAux.t(); 
      }
    }
    
    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if(i%(2*Thin) == 0 & PrintProgress == 1)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "MCMC iteration " << i << " out of " << N << " has been completed." << std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "Current draws of some selected parameters are displayed below." << std::endl;
      Rcpp::Rcout << "mu (gene 1): " << muAux(0,0) << std::endl; 
      Rcpp::Rcout << "delta (gene 1): " << deltaAux(0,0) << std::endl; 
      Rcpp::Rcout << "phi (cell 1): " << phiAux(0) << std::endl;
      Rcpp::Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
      Rcpp::Rcout << "theta (batch 1): " << thetaAux(0,0) << std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
      Rcpp::Rcout << "LSmu (gene 1): " << LSmuAuxExtra + LSmuAux(0) << std::endl;
      Rcpp::Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl; 
      Rcpp::Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
      Rcpp::Rcout << "LStheta (batch 1): " << LSthetaAux(0) << std::endl;
    }    
  }
  
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "All " << N << " MCMC iterations have been completed." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Please see below a summary of the overall acceptance rates." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;

  Rcpp::Rcout << "Minimum acceptance rate among mu[i]'s: " << min(muAccept.elem(find(Index_arma != RefGene))/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << mean(muAccept.elem(find(Index_arma != RefGene))/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among mu[i]'s: " << max(muAccept.elem(find(Index_arma != RefGene))/(N-Burn)) << std::endl;
  
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among delta[i]'s: " << min(deltaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << mean(deltaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << max(deltaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among nu[jk]'s: " << min(nuAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among nu[jk]'s: " << mean(nuAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among nu[jk]'s: " << max(nuAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among theta[k]'s: " << min(thetaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among theta[k]'s: " << mean(thetaAccept/(N-Burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among theta[k]'s: " << max(thetaAccept/(N-Burn)) << std::endl;  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  
  if(StoreAdapt == 1)
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu") = mu,
        Rcpp::Named("delta") = delta,
        Rcpp::Named("phi") = phi,
        Rcpp::Named("nu") = nu,
        Rcpp::Named("theta") = theta,
        Rcpp::Named("ls.mu") = LSmu,
        Rcpp::Named("ls.delta") = LSdelta,
        Rcpp::Named("ls.nu") = LSnu,
        Rcpp::Named("ls.theta") = LStheta,
        Rcpp::Named("RefFreq") = RefFreq/(N-Burn))); 
  }
  
  else
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu") = mu,
        Rcpp::Named("delta") = delta,
        Rcpp::Named("phi") = phi,
        Rcpp::Named("nu") = nu,
        Rcpp::Named("theta") = theta,
        Rcpp::Named("RefFreq") = RefFreq/(N-Burn))); 
  }
  
}

// [[Rcpp::export]]
arma::mat HiddenBASiCS_DenoisedRates(
  NumericMatrix CountsBio, 
  NumericMatrix Mu,
  NumericMatrix TransInvDelta,
  NumericMatrix PhiNu, 
  int N,
  int q0,
  int n)
{
  // Transformations to arma objects
  arma::mat CountsBio_arma = as_arma(CountsBio);
  arma::mat Mu_arma = as_arma(Mu);
  arma::mat TransInvDelta_arma = as_arma(TransInvDelta);
  arma::mat PhiNu_arma = as_arma(PhiNu);
  
  // Where to store the results
  arma::mat Rho = arma::zeros(q0, n);
  // Auxiliary matrices
  arma::mat m1; arma::mat m2;
  
  for (int i = 0; i<N; i++) 
  {
    Rcpp::checkUserInterrupt();
    
    m1 = CountsBio_arma; 
    m1.each_col() += TransInvDelta_arma.col(i); 
    m2 = Mu_arma.row(i).t() * PhiNu_arma.row(i);
    m2.each_col() += TransInvDelta_arma.col(i); 
    Rho += m1 / m2; 
  } 
  return(Rho / N);
}
  
  
  

