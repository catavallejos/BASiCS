/* C++ implementation of BASiCS */

#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
//File: matprod_arma.cpp
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <gperftools/profiler.h>

using namespace Rcpp;

//#include <mach/mach_time.h>
//#define ORWL_NANO (+1.0E-9)
//#define ORWL_GIGA UINT64_C(1000000000)

//static double orwl_timebase = 0.0;
//static uint64_t orwl_timestart = 0;

//struct timespec orwl_gettime(void) {
  // be more careful in a multithreaded environement
//  if (!orwl_timestart) {
//    mach_timebase_info_data_t tb = { 0 };
//    mach_timebase_info(&tb);
//    orwl_timebase = tb.numer;
//    orwl_timebase /= tb.denom;
//    orwl_timestart = mach_absolute_time();
//  }
//  struct timespec t;
//  double diff = (mach_absolute_time() - orwl_timestart) * orwl_timebase;
//  t.tv_sec = diff * ORWL_NANO;
//  t.tv_nsec = diff - (t.tv_sec * ORWL_GIGA);
//  return t;
//}

#define ZTOL sqrt(DOUBLE_EPS)

/* 
 * Rgig is an adaptation of the rgig.c function implemented by 
 * Ester Pantaleo and Robert B. Gramacy, 2010
 * (which in turn, was adapted from the C code in the ghyp v_1.5.2 package for R, rgig function)
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
double zeroin_gig(double ax,double bx,double f(double x, double m, double beta, double lambda),double tol,double m,double beta,double lambda)  /* An estimate to the root  */
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
	a = b;  b = c;  c = a;          /* best approximation		*/
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

    if( fabs(new_step) < tol_act ) {	/* Adjust the step to be not less*/
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

/* Random sample generator from a Generalized Inverse Gaussian distribution (modified version of rgig.c using Rcpp classes)*/
NumericVector Rgig(const int n, const double lambda, const double chi, const double psi)
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

  for (i=0; i<n; i++) {
    need = 1;
    while (need) {
      R1 = unif_rand();
      R2 = unif_rand();

      Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
      if (Y > 0.0) {
	  if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) {
	    need = 0;
	  }
      }
    }
    samps[i] = Y*alpha;
  }
  return(samps);
}

/* Random sample generator from a Generalized Inverse Gaussian distribution (modified version of rgig.c using Rcpp classes)*/
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
    while (need) {
      R1 = unif_rand();
      R2 = unif_rand();

      Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
      if (Y > 0.0) {
    if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) {
      need = 0;
	  }
      }
    }
    samps = Y*alpha;
  return(samps);
}

/* Auxiliary function that converts Rcpp::NumericVector elements into arma::vec elements */
arma::vec as_arma(
  NumericVector& x) /* Vector to be converted into an arma::vec class */
{
  return arma::vec(x.begin(), x.length(), false);
}

/* Auxiliary function that converts Rcpp::NumericMatrix elements into arma::mac elements */
arma::mat as_arma(NumericMatrix& x) /* Matrix to be converted into an arma::mat class */
{
  return arma::mat(x.begin(), x.nrow(), x.ncol(), false);
}

/* Auxiliary function: logarithm of a Gamma function applied element-wise to an arma::mat element */
arma::mat lgamma_cpp(arma::mat const& x) 
{
/*
 * NumericMatrix x: matrix to which the logarithm of the Gamma function will be applied element-wise
 */
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

/* Auxiliary function: logarithm of a Gamma function applied element-wise to an arma::mat element */
arma::vec lgamma_cpp_vec(arma::vec const& x) 
{
/*
 * NumericMatrix x: matrix to which the logarithm of the Gamma function will be applied element-wise
 */
  arma::vec output = x;
  for (unsigned int i=0; i < output.size(); i++)
  {
      output(i) = R::lgammafn(x(i));
  }
  return output;
}

/* Auxiliary function: to avoid numerical overflows/underflows when computing a sum of exponentiated values */
double log_sum_exp_cpp(arma::vec const& x_arma) 
{
  double offset;
  if ( max(abs(x_arma)) > max(x_arma) ) { offset = min(x_arma);}
  else { offset = max(x_arma);}
  return log(sum(exp(x_arma - offset))) + offset; 
}

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdate(
  arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */  
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
  arma::vec const& sum_bycell_bio, /* Sum of expression counts by cell (biological genes only) */
  double const& s2_mu,
  int const& q,
  int const& q_bio,
  int const& n,
  arma::vec & mu,
  arma::vec & ind)
{
    using arma::span;
    
    // PROPOSAL STEP    
    arma::vec y = exp(arma::randn(q_bio) % sqrt(prop_var(span(0,q_bio - 1))) + log(mu0(span(0, q_bio - 1))));
    arma::vec u = arma::randu(q_bio);
    
    // ACCEPT/REJECT STEP
    arma::vec log_aux = (log(y) - log(mu0(span(0, q_bio - 1)))) % sum_bycell_bio; 

    // Loop to replace matrix operations, through genes and cells
    for (int i=0; i < q_bio; i++)
    {
      for (int j=0; j < n; j++) 
      {
        log_aux(i) -= ( Counts(i,j) + (1/delta(i)) ) *  
                      log( ( phi(j)*nu(j)*y(i)+(1/delta(i)) ) / ( phi(j)*nu(j)*mu0(i)+(1/delta(i)) ));
      }
    }
    
    // THERE IS A -1 COMMING FROM THE L0G-NORMAL PRIOR. IT CANCELS OUT AS PROPOSING IN THE LOG-SCALE. 
    log_aux -= (0.5/s2_mu) * (pow(log(y),2) - pow(log(mu0(span(0, q_bio - 1))),2));

    // CREATING OUTPUT VARIABLE & DEBUG
    for (int i=0; i < q_bio; i++)
    {
      if(arma::is_finite(log_aux(i)))
      {
        if(log(u(i)) < log_aux(i))
        {
          // DEBUG: Reject very small values to avoid numerical issues
          if(arma::is_finite(y(i)) & (y(i) > 1e-3)) {ind(i) = 1; mu(i) = y(i);}
          else{ind(i) = 0; mu(i) = mu0(i);}            
        }
        else{ind(i) = 0; mu(i) = mu0(i);}
      }      
      // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
      // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
      else
      {
        ind(i) = 0; mu(i) = mu0(i);
        Rcpp::Rcout << "Something went wrong when updating mu " << i+1 << std::endl;
        Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");        
      }
    }
    
    // OUTPUT
    return join_rows(mu, ind);
}

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat deltaUpdate(
  arma::vec const& delta0, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& prop_var,  /* Current value of the proposal variances for $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
  double const& a_delta, /* Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  double const& b_delta, /* Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  int const& q_bio,
  int const& n,
  double const& s2_delta,
  double const& prior_delta)
{
    using arma::span;
    
    // CREATING VARIABLES WHERE TO STORE DRAWS
    arma::vec delta = arma::zeros(q_bio); 
    arma::vec ind = arma::zeros(q_bio);

    // PROPOSAL STEP
    arma::vec y = exp(arma::randn(q_bio) % sqrt(prop_var) + log(delta0));
    arma::vec u = arma::randu(q_bio);
        
    // ACCEPT/REJECT STEP 
    arma::vec log_aux = - n * (lgamma_cpp(1/y)-lgamma_cpp(1/delta0));
    
    // Loop to replace matrix operations, through genes and cells
    for (int i=0; i < q_bio; i++)
    {
      for (int j=0; j < n; j++) 
      {
          log_aux(i) += R::lgammafn(Counts(i,j) + (1/y(i))) - R::lgammafn(Counts(i,j) + (1/delta0(i)));
          log_aux(i) -= ( Counts(i,j) + (1/y(i)) ) *  log( phi(j)*nu(j)*mu(i)+(1/y(i)) );
          log_aux(i) += ( Counts(i,j) + (1/delta0(i)) ) *  log( phi(j)*nu(j)*mu(i)+(1/delta0(i)) );
      }
    }
    
    // +1 should appear because we update log(delta) not delta. However, it cancels out with the prior. 
    log_aux -= n * ( (log(y)/y) - (log(delta0)/delta0) );
    // Component related to the prior
    if(prior_delta == 1) {log_aux += (log(y) - log(delta0)) * a_delta - b_delta * (y - delta0);}
    else {log_aux -= (0.5/s2_delta) * (pow(log(y),2) - pow(log(delta0),2));}
    
    // CREATING OUTPUT VARIABLE & DEBUG
    for (int i=0; i < q_bio; i++)
    {
      if(arma::is_finite(log_aux(i)))
      {
        if(log(u(i)) < log_aux(i))
        {
          // DEBUG: Reject very small values to avoid numerical issues
          if(arma::is_finite(y(i)) & (y(i) > 1e-3)) {ind(i) = 1; delta(i) = y(i);}
          else{ind(i) = 0; delta(i) = delta0(i);}            
        }
        else{ind(i) = 0; delta(i) = delta0(i);}
      }      
      // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
      // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
      else
      {
        ind(i) = 0; delta(i) = delta0(i);
        Rcpp::Rcout << "Something went wrong when updating delta " << i+1 << std::endl;
        Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");        
      }
    }
    
    // OUTPUT
    return join_rows(delta, ind);
}

/* Draws for cell-specific normalising constants s[j]
 * Metropolis-Hastings updates are not required as full conditionals have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::vec sUpdate(
  arma::vec const& s0, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& theta, /* Current value of $\theta$ */
  double const& as, /* Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double const& bs, /* Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  int const& n)
{   
   // Creating a vector where draws will be stored
   arma::vec aux = arma::zeros(n);
   
   // Calculating parameters to the passed as input to the Rgig function (common for all cells)
   double p = as-1/theta; 
   double b = 2*bs;
   
   // GIG draws
   /* DEBUG: return original value of s0 if input values are not within the appropriate range (to avoid problems related to numerical innacuracy) */
   if(!R_IsNA(p))  
   {     
      for (int j=0; j<n; j++) 
      {
        if(nu(j) > 0) 
        {
          aux(j) = RgigDouble(p,2*nu(j)/theta,b);
        }
        /* DEBUG: break in case of undefined values */
        else
        {
          aux(j) = s0(j);
          Rcpp::Rcout << "Something went wrong when updating s" << j << std::endl;
          Rcpp::stop("Please consider additional filter of the input dataset.");         
        }               
      }
      return aux;     
   }
  else return s0;
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
arma::mat nuUpdate(
  arma::vec const& nu0, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\nu=(\nu_1,...,\nu_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  double const& SumSpikeInput,
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */  
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  double const& theta, /* Current value of $\theta$' */
  arma::vec const& sum_bygene_all, /* Sum of expression counts by gene (all genes) */
  int const& q,
  int const& q_bio,
  int const& n)
{
    using arma::span;
    
    // CREATING VARIABLES WHERE TO STORE DRAWS
    arma::vec nu = arma::zeros(n); 
    arma::vec ind = arma::zeros(n);

    // PROPOSAL STEP    
    arma::vec y = exp(arma::randn(n) % sqrt(prop_var) + log(nu0));
    arma::vec u = arma::randu(n);
    
    // ACCEPT/REJECT STEP
    arma::vec log_aux = arma::zeros(n);

    // Loop to replace matrix operations, through genes and cells
    for (int j=0; j < n; j++) 
    {
      for (int i=0; i < q_bio; i++) 
      {
        log_aux(j) -= ( Counts(i,j) + (1/delta(i)) ) *  
                      log( ( phi(j)*y(j)*mu(i)+(1/delta(i)) ) / ( phi(j)*nu0(j)*mu(i)+(1/delta(i)) ));
      } 
      log_aux(j) += (log(y(j)) - log(nu0(j))) * (sum_bygene_all(j) + 1/theta) - ( y(j)-nu0(j) )  * (SumSpikeInput + (1/(theta*s(j))));
    }

    // CREATING OUTPUT VARIABLE & DEBUG
    for (int j=0; j < n; j++)
    {
      if(arma::is_finite(log_aux(j)))
      {
        if(log(u(j)) < log_aux(j))
        {
          // DEBUG: Reject very small values to avoid numerical issues
          if(arma::is_finite(y(j)) & (y(j) > 1e-5)) {ind(j) = 1; nu(j) = y(j);}
          else{ind(j) = 0; nu(j) = nu0(j);}            
        }
        else{ind(j) = 0; nu(j) = nu0(j);}
      }      
      // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
      // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
      else
      {
        ind(j) = 0; nu(j) = nu0(j);
        Rcpp::Rcout << "Something went wrong when updating nu " << j+1 << std::endl;
        Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");        
      }
    }
   
    // OUTPUT
    return join_rows(nu, ind);
}

/* Metropolis-Hastings updates of theta 
 */
arma::vec thetaUpdate(
  double const& theta0, /* Current value of $\theta$ */
  double const& prop_var, /* Current value of the proposal variances for $\theta$ */
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& a_theta, /* Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  double const& b_theta, /* Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  int const& n)
  {
    using arma::span;

    // CREATING VARIABLES WHERE TO STORE DRAWS
    double logtheta = log(theta0);

    // PROPOSAL STEP
    double y = as<double>(wrap(arma::randn(1))) * sqrt(prop_var) + logtheta;
    double u = as<double>(wrap(arma::randu(1)));
    
    // ACCEPT/REJECT STEP
    double log_aux = (y-logtheta)*a_theta;
    log_aux -= n * (logtheta/theta0) * ((y/logtheta)*exp(-y+logtheta)-1) + n*(R::lgammafn(exp(-y))-R::lgammafn(1/theta0));
    log_aux += ((exp(-y+logtheta)-1)/theta0) * sum(log(nu/s)-(nu/s)) - b_theta*theta0*(exp(y-logtheta)-1);
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

/*Uniform sample from a discrete range*/
double rUnifDisc(arma::vec const& range, double const& rangeSize)
{
  double u = as<double>(wrap(arma::randu(1)));
  double out;
  
  arma::vec probs = arma::ones(rangeSize);
  probs = cumsum(probs) / rangeSize;
  if(u <= probs(0)) {out = range(0);}
  else
  {
    for(int i=1; i < rangeSize; i++)
    {
      if(u <= probs(i) & u > probs(i-1)) {out = range(i);}
    } 
  } 
  return out;
}

/*Uniform sample from a discrete range, multiple draws*/
arma::vec rUnifDiscMult(int const& ndraws, arma::vec const& range, double const& rangeSize)
{
  arma::vec out = arma::zeros(ndraws);

  for(int i=0; i < ndraws; i++)
  {
    out(i) = rUnifDisc(range, rangeSize);
  } 
  return out;
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

/* Auxiliary function required Metropolis-Hastings updates of phi */
arma::mat UpdateAux_phiTrick(
  arma::vec const& mu, /* Current value of mu or mu*exp(tau), depending on the group) */
  arma::vec const& delta, /* Current value of $\delta$ or $\delta \exp(omega_{bio})$, depending on the group)$*/
  arma::vec const& phinu) /* Current value of $\phi \nu$*/
{  
  arma::mat x = (phinu * mu.t()).t();
  x.each_col() += 1 / delta;
  
  return(log(x));
}

/* Metropolis-Hastings updates of phi 
 * Joint updates using Dirichlet proposals
 */
Rcpp::List phiUpdate(
  arma::vec const& phi0, // Current value of $\phi=(\phi_1,...,\phi_n)'$
  double const& prop_var, // Current value of the proposal precision
  arma::mat const& Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom)
  arma::vec const& mu, // Current value of $\mu=(\mu_1,...,\mu_q)'$
  arma::vec const& delta, // Current value of $\delta=(\delta_1,...,\delta_{q_0})'$
  arma::vec const& nu, // Current value of $\nu=(\nu_1,...,\nu_n)'$
  arma::vec const& p_phi, // Dirichlet hyper-parameter of the prior for $\phi / n$
  arma::vec const& sum_bygene_bio, // Sum of expression counts by gene (biological genes only)
  int const& q_bio, // Number of biological genes
  int const& n) // Total number of cells 
{
  using arma::span;
  using arma::pow;
    
  // PROPOSAL STEP
  arma::vec y = n * rDirichlet(prop_var * phi0); 
  double u = R::runif(0,1);
  
  arma::vec phi; 
  int ind;
    
  // ACCEPT/REJECT STEP (REJECT VALUES OUTSIDE VALID RANGE)  
  if(all(prop_var * y < 2.5327372760800758e+305)  & 
     all(prop_var * phi0 < 2.5327372760800758e+305) &
     all(y > 0) &
     all(phi0 > 0)) 
  {
    // There is an extra -1 but it cancels out with the proposal component
    double log_aux = sum( (sum_bygene_bio + p_phi) % (log(y) - log(phi0)));
        
    // Loop to replace matrix operations, through genes and cells
    // There is an extra factor in the prior n^(-n) but it cancels out in the ratio
    // There is an extra factor n^(-(sum(p_phi) - 1)) but it cancels out in the ratio
    for (int j=0; j < n; j++) 
    {
      for (int i=0; i < q_bio; i++) 
      {
        log_aux -= ( Counts(i,j) + (1/delta(i)) ) *  
                    log( (y(j)*nu(j)*mu(i)+(1/delta(i)) ) / ( phi0(j)*nu(j)*mu(i)+(1/delta(i)) ));
      } 
    }
    // There is an extra factor n^(-(sum(prop_var * y))) / n^(-(sum(prop_var * phi0))), it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
    // There is an extra factor gamma(sum(prop_var * phi0)) / gamma(sum(prop_var * y)), it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
    log_aux += prop_var * sum(y % log(phi0) - phi0 % log(y));
    log_aux -= sum(lgamma_cpp_vec(prop_var * y) - lgamma_cpp_vec(prop_var * phi0));    
    
    if(!R_IsNA(log_aux))
    {
        if(log(u) < log_aux) {ind = 1;}
        else {ind = 0;}
        phi = ind * y + (1-ind) * phi0;
    }
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
    else
    {
      Rcpp::Rcout << "Something went wrong when updating phi" << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset."); 
      ind = 0;
      phi = phi0;
    }     
  }
  else
  {
     ind = 0;
     phi = phi0;     
  }      
  return(Rcpp::List::create(
         Rcpp::Named("phi")=phi,
         Rcpp::Named("ind")=ind)); 
}




/* MCMC sampler 
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcpp(
  int N, // Total number of MCMC draws 
  int thin, // Thinning period for MCMC chain 
  int burn, // Burning period for MCMC chain 
  NumericMatrix Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom) 
  NumericVector mu0, // Starting value of $\mu=(\mu_1,...,\mu_q)'$ (true mRNA content for technical genes)  
  NumericVector delta0, // Starting value of $\delta=(\delta_1,...,\delta_{q_0})'$  
  NumericVector phi0, // Starting value of $\phi=(\phi_1,...,\phi_n)$'$ 
  NumericVector s0, // Starting value of $s=(s_1,...,s_n)$'$ 
  NumericVector nu0, // Starting value of $\nu=(\nu_1,...,\nu_n)$'$   
  double theta0, // Starting value of $\theta$ 
  double s2mu,
  double adelta, // Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  double bdelta, // Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  NumericVector p_Phi, // Dirichlet hyper-parameter for $\phi / n$ 
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
  double s2_delta, 
  double prior_delta) 
{

  
// NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int q = Counts.nrow(); int qbio = delta0.size(); int Naux = N/thin - burn/thin;
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts);
  arma::vec mu0_arma = as_arma(mu0);
  arma::vec phi0_arma = as_arma(phi0);
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = arma::zeros(Naux,q); 
  arma::mat delta = arma::zeros(Naux,qbio); 
  arma::mat phi = arma::ones(Naux,n);
  arma::mat s = arma::zeros(Naux,n);  
  arma::mat nu = arma::zeros(Naux,n); 
  arma::vec theta = arma::zeros(Naux); 
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::vec LSphi;
  arma::mat LSnu;
  arma::vec LStheta; 

  // LOG-PROPOSAL VARIANCES 
  if(StoreAdapt == 1)
  {
    LSmu = arma::zeros(Naux,q); 
    LSdelta = arma::zeros(Naux,qbio); 
    LSphi = arma::ones(Naux);   
    LSnu = arma::zeros(Naux,n); 
    LStheta = arma::zeros(Naux);   
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(q); arma::vec PmuAux = arma::zeros(q);
  arma::vec deltaAccept = arma::zeros(qbio); arma::vec PdeltaAux = arma::zeros(qbio);
  double phiAccept = 0; double PphiAux = 0;
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  double thetaAccept=0; double PthetaAux=0;
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(q,2); muAux.col(0)=mu0_arma; arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(qbio,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::vec phiAux = phi0_arma; double LSphiAux = LSphi0; Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0); 
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::vec thetaAux = arma::zeros(2); thetaAux(0) = theta0; double LSthetaAux = LStheta0;  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q); arma::vec PdeltaAux0 = arma::zeros(qbio);
  double PphiAux0 = 0; 
  arma::vec PnuAux0 = arma::zeros(n); double PthetaAux0=0;
  
  double SumSpikeInput = sum(mu0_arma(arma::span(qbio,q -1)));
  

  
  NumericVector LSphiValuesNV = NumericVector::create(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 50000000);
  arma::vec LSphiValues = as_arma(LSphiValuesNV);
  NumericVector LSthetaValuesNV = NumericVector::create(-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2);
  arma::vec LSthetaValues = as_arma(LSthetaValuesNV);
  
  
  
                                                          
    
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;
  
  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  arma::vec muUpdateAux = arma::ones(q); muUpdateAux(arma::span(qbio, q-1)) = mu0_arma(arma::span(qbio, q-1));
  arma::vec indQ = arma::zeros(q);
  
  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;  
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  
  // START OF MCMC LOOP
  for (i=0; i<N; i++) {
    
//    Rcpp::checkUserInterrupt();
    
    if(i==burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
      Rcpp::Rcout << "End of burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
    }
    
    Ibatch++; 
    
    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
//    LSphiAux = rUnifDisc(LSphiValues, LSphiValuesNV.size());    
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, muAux.col(0), deltaAux.col(0),
                           nuAux.col(0), p_Phi, sumByGeneBio_arma, qbio, n); 
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]); if(i>=burn) {phiAccept += Rcpp::as<double>(phiAuxList["ind"]);}

    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
//    LSthetaAux = rUnifDisc(LSthetaValues, LSthetaValuesNV.size());  
    thetaAux = thetaUpdate(thetaAux(0), exp(LSthetaAux), sAux, nuAux.col(0), atheta, btheta, n);
    PthetaAux += thetaAux(1); if(i>=burn) {thetaAccept += thetaAux(1);}
    
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR       
//    LSmuAux = rUnifDiscMult(q, LSthetaValues, LSthetaValuesNV.size());
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, deltaAux.col(0), 
                     phiAux, nuAux.col(0), sumByCellBio_arma, s2mu, q, qbio, n,
                     muUpdateAux, indQ);     
    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}
    
    

    
    // UPDATE OF S
    sAux = sUpdate(sAux, nuAux.col(0), thetaAux(0), as, bs, n); 
    
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
//    LSdeltaAux = rUnifDiscMult(qbio, LSthetaValues, LSthetaValuesNV.size());
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, muAux.col(0), 
                           nuAux.col(0), phiAux, adelta, bdelta, qbio, n, s2_delta, prior_delta);  
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);}    
    
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdate(nuAux.col(0), exp(LSnuAux), Counts_arma, SumSpikeInput,
                     muAux.col(0), deltaAux.col(0), phiAux, sAux, thetaAux(0), sumByGeneAll_arma, q, qbio, n); 
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/50; PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
//        LSmuAux=LSmuAux+PmuAux*std::min(0.03,1/sqrt(i)); 
        LSmuAux=LSmuAux+PmuAux*0.1;
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
//        LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.03,1/sqrt(i));
        LSdeltaAux=LSdeltaAux+PdeltaAux*0.1;
//        if(i < 1000) {LSdeltaAux=LSdeltaAux+PdeltaAux*0.5;}     
//        else {LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.01,1/sqrt(i));}  
        PphiAux=PphiAux/50; PphiAux = -1+2*(PphiAux>ar); 
        LSphiAux=LSphiAux - PphiAux*0.1;  
//        LSphiAux=LSphiAux - PphiAux*std::min(0.03,1/sqrt(i));  
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
//        LSnuAux=LSnuAux+PnuAux*std::min(0.03,1/sqrt(i)); 
        LSnuAux=LSnuAux+PnuAux*0.1;
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*(PthetaAux>ar); 
//        LSthetaAux=LSthetaAux+PthetaAux*std::min(0.03,1/sqrt(i));
        LSthetaAux=LSthetaAux+PthetaAux*0.1;
      
        Ibatch = 0; 
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0; 
        PphiAux = PphiAux0; 
        PnuAux = PnuAux0; PthetaAux = PthetaAux0; 
      }
      
    }
        
    // STORAGE OF DRAWS
    if((i%thin==0) & (i>=burn))
    {      
      mu.row(i/thin - burn/thin) = muAux.col(0).t(); 
      delta.row(i/thin - burn/thin) = deltaAux.col(0).t(); 
      phi.row(i/thin - burn/thin) = phiAux.t();
      s.row(i/thin - burn/thin) = sAux.t();
      nu.row(i/thin - burn/thin) = nuAux.col(0).t();       
      theta(i/thin - burn/thin) = thetaAux(0);       
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/thin - burn/thin) = LSmuAux.t();
        LSdelta.row(i/thin - burn/thin) = LSdeltaAux.t();
        LSphi(i/thin - burn/thin) = LSphiAux;
        LSnu.row(i/thin - burn/thin) = LSnuAux.t();
        LStheta(i/thin - burn/thin) = LSthetaAux; 
      }
    }
    
    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if((i%(2*thin) == 0) & (PrintProgress == 1))
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
        Rcpp::Rcout << "theta: " << thetaAux(0) << std::endl;
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
        Rcpp::Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
        Rcpp::Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl; 
        Rcpp::Rcout << "LSphi: " << LSphiAux << std::endl;
        Rcpp::Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
        Rcpp::Rcout << "LStheta: " << LSthetaAux << std::endl;
    }    
  }
  
//  ProfilerStop();
  
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
  
  Rcpp::Rcout << "Minimum acceptance rate among mu[i]'s: " << min(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << mean(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among mu[i]'s: " << max(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among delta[i]'s: " << min(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << mean(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among delta[i]'s: " << max(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for phi (joint): " << phiAccept/(N-burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among nu[j]'s: " << min(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among nu[j]'s: " << mean(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among nu[j]'s: " << max(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for theta: " << thetaAccept/(N-burn) << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;

      
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
   // Allowing the use of rgig within C++ code (transform arma::vec elements into NumericVector elements)
//   NumericVector s0 = wrap(s0_arma);
//   NumericVector nu = wrap(nu_arma);
   arma::vec thetaBatch = BatchDesign * theta; 
   
   // Creating a vector where draws will be stored
   arma::vec aux = arma::zeros(n);
   
   // Calculating parameters to the passed as input to the Rgig function (common for all cells)
   arma::vec p = as - 1 / thetaBatch; 
   double b = 2 * bs;
   
   // GIG draws
   /* DEBUG: return original value of s0 if input values are not within the appropriate range (to avoid problems related to numerical innacuracy) */

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
          Rcpp::Rcout << "Something went wrong when updating s" << j << std::endl;
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
  int const& q,
  int const& q_bio,
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
//      for (int i=0; i < q_bio; i++) 
//      {
//        log_aux(j) -= ( Counts(i,j) + (1/delta(i)) ) *  
//                      log( ( phi(j)*exp(y(j))*mu(i)+(1/delta(i)) ) / ( phi(j)*nu0(j)*mu(i)+(1/delta(i)) ));
//      } 
//      log_aux(j) += (y(j) - lognu(j)) * (sum_bygene_all(j) + 1/thetaBatch(j)); // - ( y(j)-nu0(j) )  * (SumSpikeInput + (1/(thetaBatch(j)*s(j))));
//      log_aux(j) -= nu0(j) * (exp(y(j)-lognu(j))-1)  * (SumSpikeInput + (1/(thetaBatch(j) * s(j))));
//    }

    
    arma::vec log_aux = (y - lognu) % (sum_bygene_all + 1/thetaBatch);
    log_aux -= nu0 % (exp(y-lognu)-1)  % (SumSpikeInput + (1/(thetaBatch % s)));   
    arma::mat m = Counts.rows(0, q_bio - 1);
    m.each_col() += 1 / delta;   
    arma::mat num = UpdateAux_nuTrick(exp(y), mu(span(0, q_bio - 1)), delta,phi);
    num /= UpdateAux_nuTrick(nu0,    mu(span(0, q_bio - 1)), delta,phi);
    m %= log(num);
    log_aux -= sum(m, 0).t();
    arma::umat ind = log(u) < log_aux;
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message    
    ind.elem(find_nonfinite(log_aux)).zeros();
    if(size(find_nonfinite(log_aux),0)>0)
    {
      Rcpp::Rcout << "Something went wrong when updating nu" << size(find_nonfinite(log_aux),0) << std::endl;
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
  int thin, // Thinning period for MCMC chain 
  int burn, // Burning period for MCMC chain 
  NumericMatrix Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom) 
  NumericMatrix BatchDesign, // Design matrix representing batch information (number of columns must be equal to number of batches)
  NumericVector mu0, // Starting value of $\mu=(\mu_1,...,\mu_q)'$ (true mRNA content for technical genes)  
  NumericVector delta0, // Starting value of $\delta=(\delta_1,...,\delta_{q_0})'$  
  NumericVector phi0, // Starting value of $\phi=(\phi_1,...,\phi_n)$'$ 
  NumericVector s0, // Starting value of $s=(s_1,...,s_n)$'$ 
  NumericVector nu0, // Starting value of $\nu=(\nu_1,...,\nu_n)$'$   
  double theta0, // Starting value of $\theta$ 
  double s2mu, 
  double adelta, // Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  double bdelta, // Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  NumericVector p_Phi, // Dirichlet hyper-parameter for $\phi / n$ 
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
  double s2_delta, 
  double prior_delta) 
{

  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int q = Counts.nrow(); int qbio = delta0.size(); int Naux = N/thin - burn/thin;
  int nBatch = BatchDesign.ncol();
  
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts); arma::mat BatchDesign_arma = as_arma(BatchDesign);
  arma::vec mu0_arma = as_arma(mu0);
  arma::vec  phi0_arma = as_arma(phi0);
  
  double SumSpikeInput = sum(mu0_arma(arma::span(qbio,q -1)));
  
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = arma::zeros(Naux, q); 
  arma::mat delta = arma::zeros(Naux, qbio); 
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
    LSmu = arma::zeros(Naux, q); 
    LSdelta = arma::zeros(Naux, qbio); 
    LSphi = arma::ones(Naux);   
    LSnu = arma::zeros(Naux, n); 
    LStheta = arma::zeros(Naux, nBatch);   
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(q); arma::vec PmuAux = arma::zeros(q);
  arma::vec deltaAccept = arma::zeros(qbio); arma::vec PdeltaAux = arma::zeros(qbio);
  double phiAccept = 0; double PphiAux = 0;
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  arma::vec thetaAccept = arma::zeros(nBatch); arma::vec PthetaAux = arma::zeros(nBatch);
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(q,2); muAux.col(0)=as_arma(mu0); arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(qbio,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::vec phiAux = phi0_arma; double LSphiAux = LSphi0; Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0); 
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::mat thetaAux = arma::zeros(nBatch, 2); thetaAux.col(0) = theta0 * arma::ones(nBatch); 
  arma::vec LSthetaAux = LStheta0 * arma::ones(nBatch);  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q); arma::vec PdeltaAux0 = arma::zeros(qbio);
  double PphiAux0 = 0; 
  arma::vec PnuAux0 = arma::zeros(n); arma::vec PthetaAux0 = arma::ones(nBatch);
  
//  arma::vec ThetaExample = thetaAux.col(0);
//  arma::vec Example = BatchDesign_arma * ThetaExample; 
//  Rcpp::Rcout << Example << std::endl;
  
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;
  
  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  arma::vec muUpdateAux = arma::ones(q); muUpdateAux(arma::span(qbio, q-1)) = mu0_arma(arma::span(qbio, q-1));
  arma::vec indQ = arma::zeros(q);
  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;  
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  
  // START OF MCMC LOOP
  for (i=0; i<N; i++) {
    
    Rcpp::checkUserInterrupt();
    
    if(i==burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
      Rcpp::Rcout << "End of burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
    }
    
    Ibatch++; 
    
//    struct timespec time0_1 = orwl_gettime();
    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, muAux.col(0), deltaAux.col(0),
                           nuAux.col(0), p_Phi, sumByGeneBio_arma, qbio,n); 
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]); if(i>=burn) {phiAccept += Rcpp::as<double>(phiAuxList["ind"]);}
//    struct timespec time1_1 = orwl_gettime();    
//    Rcpp::Rcout << "Time phi: " << (time1_1.tv_nsec - time0_1.tv_nsec) / ((float)(n)) << std::endl;

//    struct timespec time0_2 = orwl_gettime();
    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdateBatch(thetaAux.col(0), exp(LSthetaAux), BatchDesign_arma,
                                sAux, nuAux.col(0), atheta, btheta, n, nBatch);
    PthetaAux += thetaAux.col(1); if(i>=burn) {thetaAccept += thetaAux.col(1);}
//    struct timespec time1_2 = orwl_gettime();    
//    Rcpp::Rcout << "Time theta: "  <<  (time1_2.tv_nsec - time0_2.tv_nsec) / ((float)(nBatch)) << std::endl;
    
//    struct timespec time0_3 = orwl_gettime();
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR 
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, deltaAux.col(0), 
                     phiAux, nuAux.col(0), sumByCellBio_arma, s2mu, q, qbio, n,
                     muUpdateAux, indQ);     
    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}
//    struct timespec time1_3 = orwl_gettime();    
//    Rcpp::Rcout << "Time mu: "  <<  (time1_3.tv_nsec - time0_3.tv_nsec) / ((float)(qbio)) << std::endl;
    
//    struct timespec time0_4 = orwl_gettime();
    // UPDATE OF S
    sAux = sUpdateBatch(sAux, nuAux.col(0), thetaAux.col(0), as, bs, BatchDesign_arma, n);  
//    struct timespec time1_4 = orwl_gettime();    
//    Rcpp::Rcout << "Time s: "  <<  (time1_4.tv_nsec - time0_4.tv_nsec) / ((float)(n)) << std::endl;
    
//    struct timespec time0_5 = orwl_gettime();
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, 
                           muAux.col(0), nuAux.col(0), phiAux, adelta, bdelta, qbio, n, s2_delta, prior_delta);  
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);} 
//    struct timespec time1_5 = orwl_gettime();    
//    Rcpp::Rcout << "Time delta: "  << (time1_5.tv_nsec - time0_5.tv_nsec) / ((float)(qbio)) << std::endl;
    
//    struct timespec time0_6 = orwl_gettime();
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateBatch(nuAux.col(0), exp(LSnuAux), Counts_arma, SumSpikeInput,
                     BatchDesign_arma,
                     muAux.col(0), deltaAux.col(0),
                     phiAux, sAux, thetaAux.col(0), sumByGeneAll_arma, q, qbio, n); 
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}
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
    if(i%thin==0 & i>=burn)
    {      
      mu.row(i/thin - burn/thin) = muAux.col(0).t(); 
      delta.row(i/thin - burn/thin) = deltaAux.col(0).t(); 
      phi.row(i/thin - burn/thin) = phiAux.t();
      s.row(i/thin - burn/thin) = sAux.t();
      nu.row(i/thin - burn/thin) = nuAux.col(0).t();       
      theta.row(i/thin - burn/thin) = thetaAux.col(0).t();   
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/thin - burn/thin) = LSmuAux.t();
        LSdelta.row(i/thin - burn/thin) = LSdeltaAux.t();
        LSphi(i/thin - burn/thin) = LSphiAux;
        LSnu.row(i/thin - burn/thin) = LSnuAux.t();
        LStheta.row(i/thin - burn/thin) = LSthetaAux.t(); 
      }
    }
    
    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if(i%(2*thin) == 0 & PrintProgress == 1)
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
  
  Rcpp::Rcout << "Minimum acceptance rate among mu[i]'s: " << min(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << mean(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among mu[i]'s: " << max(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among delta[i]'s: " << min(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << mean(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << max(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for phi (joint): " << phiAccept/(N-burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among nu[j]'s: " << min(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among nu[j]'s: " << mean(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among nu[j]'s: " << max(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among theta's: " << min(thetaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among theta's: " << mean(thetaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among theta's: " << max(thetaAccept/(N-burn)) << std::endl;  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
      
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

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdateNoSpikes(
  arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q_0)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */  
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
  arma::vec const& sum_bycell_all, /* Sum of expression counts by cell (biological genes only) */
  double const& s2_mu,
  int const& q_bio,
  int const& n,
  arma::vec & mu,
  arma::vec & ind)
{
  using arma::span;
  
  // PROPOSAL STEP    
  arma::vec y = exp(arma::randn(q_bio) % sqrt(prop_var) + log(mu0));
  arma::vec u = arma::randu(q_bio);
  
  // ACCEPT/REJECT STEP
  arma::vec log_aux = (log(y) - log(mu0)) % sum_bycell_all; 
  
  // Loop to replace matrix operations, through genes and cells
  for (int i=0; i < q_bio; i++)
  {
    for (int j=0; j < n; j++) 
    {
      log_aux(i) -= ( Counts(i,j) + (1/delta(i)) ) *  
        log( ( nu(j)*y(i)+(1/delta(i)) ) / ( nu(j)*mu0(i)+(1/delta(i)) ));
    }
  }
  
  log_aux -= (0.5/s2_mu) * (pow(log(y),2) - pow(log(mu0),2));
  
  // CREATING OUTPUT VARIABLE & DEBUG
  for (int i=0; i < q_bio; i++)
  {
    if(arma::is_finite(log_aux(i)))
    {
      if(log(u(i)) < log_aux(i))
      {
        // DEBUG: Reject very small values to avoid numerical issues
        if(arma::is_finite(y(i)) & (y(i) > 1e-3)) {ind(i) = 1; mu(i) = y(i);}
        else{ind(i) = 0; mu(i) = mu0(i);}            
      }
      else{ind(i) = 0; mu(i) = mu0(i);}
    }      
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
    else
    {
      ind(i) = 0; mu(i) = mu0(i);
      Rcpp::Rcout << "Something went wrong when updating mu " << i+1 << std::endl;
      Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");        
    }
  }
  
  // OUTPUT
  return join_rows(mu, ind);
}

/* Multivariate normal random variable generator.
Code has been taken from Rcpp gallery: http://gallery.rcpp.org/articles/simulate-multivariate-normal/
Author: Ahmadou Dicko. Date: March 12, 2013*/
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdateNoSpikesConstrain(
    arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q_0)'$ */
    arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
    double const& Constrain,
    arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */  
    arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
    arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
    arma::vec const& sum_bycell_all, /* Sum of expression counts by cell (biological genes only) */
    double const& s2_mu,
    int const& q_bio,
    int const& n,
    arma::vec & mu,
    arma::vec & ind,
    arma::mat & InvCovMu,
    arma::vec const& mu0_const,
    int const& ref,
    arma::vec const& Index)
{
  //     arma::mat & CholCovMu,
  // 
  using arma::span;
  
  arma::uvec IndexAux = find(Index != ref); 
  
  // PROPOSAL STEP  
  arma::vec y = arma::ones(q_bio);
  y.elem(IndexAux) = exp(log(mu0.elem(IndexAux)) + sqrt(prop_var.elem(IndexAux)) % arma::randn(q_bio - 1));
  y(ref) = exp(q_bio * Constrain - sum(log(y.elem(IndexAux))));
  double u = R::runif(0,1);
  
  // There is no problem when calculating log(mu0), starting values are all above 1. 
//  arma::vec y = exp(q_bio * Constrain * rDirichlet(prop_var * log(mu0)));
//  arma::vec y = mu0_const % exp(rDirichlet(prop_var * log(mu0)));
//  arma::vec y = arma::ones(q_bio);

////  y(span(0, q_bio-2)) = exp(log(mu0(span(0, q_bio-2))) + prop_var * CholCovMu * arma::randn(q_bio - 1));
//  y(span(0, q_bio-2)) = exp(log(mu0(span(0, q_bio-2))) + sqrt(prop_var(span(0, q_bio-2))) % arma::randn(q_bio - 1));
//  y(q_bio-1) = exp(q_bio * Constrain - sum(log(y(span(0, q_bio-2)))));
  
  
  // ACCEPT/REJECT STEP
  
  // Likelihood contribution
  arma::vec log_aux = (log(y) - log(mu0)) % sum_bycell_all; 
  // Loop to replace matrix operations, through genes and cells
  for (int i=0; i < q_bio; i++)
  {
    for (int j=0; j < n; j++) 
    {
      log_aux(i) -= ( Counts(i,j) + (1/delta(i)) ) *  
        log( ( nu(j)*y(i)+(1/delta(i)) ) / ( nu(j)*mu0(i)+(1/delta(i)) ));
    }
  }
  // Prior contribution
  arma::vec aux = Constrain * arma::ones(q_bio-1);
//  double log_aux_sum = sum(log_aux(span(0, q_bio-2)));
//  log_aux_sum -= as_scalar(0.5 * (log(y(span(0, q_bio-2))) - aux).t() * InvCovMu * (log(y(span(0, q_bio-2))) - aux));
//  log_aux_sum += as_scalar(0.5 * (log(mu0(span(0, q_bio-2))) - aux).t() * InvCovMu * (log(mu0(span(0, q_bio-2))) - aux));
  double log_aux_sum = sum(log_aux.elem(IndexAux));
  log_aux_sum -= as_scalar(0.5 * (log(y.elem(IndexAux)) - aux).t() * InvCovMu * (log(y.elem(IndexAux)) - aux));
  log_aux_sum += as_scalar(0.5 * (log(mu0.elem(IndexAux)) - aux).t() * InvCovMu * (log(mu0.elem(IndexAux)) - aux));
  
//  // There is a -1 comming from the log-normal prior. It cancels out as proposing in the log-scale.
//  // There is a -1 comming from the log-normal prior. It cancels out with the -1 in the Dirichlet density
//  log_aux -= (0.5/s2_mu) * (pow(log(y),2) - pow(log(mu0),2));
//  double log_aux_sum = sum(log_aux);
  // Proposal factor
  // There is an extra factor (q_bio * Constrain)^(-(sum(prop_var * log(y)))) / (q_bio * Constrain)^(-(sum(prop_var * log(mu0)))), it cancels out as sum(prop_var*y) = sum(prop_var*phi0)
  // There is an extra factor gamma(sum(prop_var * log(mu0))) / gamma(sum(prop_var * log(y))), it cancels out as sum(prop_var*log(y)) = sum(prop_var*log(mu0))
//  log_aux_sum += prop_var * sum(log(y) % log(log(mu0)) - log(mu0) % log(log(y)));
//  log_aux_sum -= sum(lgamma_cpp_vec(prop_var * log(y)) - lgamma_cpp_vec(prop_var * log(mu0)));  
//  log_aux_sum += prop_var * sum((y / mu0_const) % log(mu0 / mu0_const) - (mu0 / mu0_const) % log(y / mu0_const));
//  log_aux_sum -= sum(lgamma_cpp_vec(prop_var * y / mu0_const) - lgamma_cpp_vec(prop_var * mu0 / mu0_const));
//  log_aux_sum += lgamma(prop_var * sum(y / mu0_const)) - lgamma(prop_var * sum(mu0 / mu0_const));
 
//  arma::vec aux = Constrain * arma::ones(q_bio-1);
//  double log_aux_sum = sum(log_aux(span(0, q_bio-2)));
//  log_aux_sum -= as_scalar(0.5 * (log(y(span(0, q_bio-2))) - aux).t() * InvCovMu * (log(y(span(0, q_bio-2))) - aux));
//  log_aux_sum += as_scalar(0.5 * (log(mu0(span(0, q_bio-2))) - aux).t() * InvCovMu * (log(mu0(span(0, q_bio-2))) - aux));
//  Rcpp::Rcout << "log_aux_sum: " << log_aux_sum << std::endl;
//  Rcpp::Rcout << "log(u): " << log(u(0)) << std::endl;
  
  // DEBUG: Reject very small values to avoid numerical issues
  if(log(u) < log_aux_sum & min(y) > 1e-3)
  {
    ind(0) = 1; mu = y;
  }
  else{ind(0) = 0; mu = mu0;}
  
  // CREATING OUTPUT VARIABLE & DEBUG
//  for (int i=0; i < q_bio; i++)
//  {
//    if(arma::is_finite(log_aux(i)))
//    {
//      if(log(u(i)) < log_aux(i))
//      {
//        // DEBUG: Reject very small values to avoid numerical issues
//        if(arma::is_finite(y(i)) & (y(i) > 1e-3)) {ind(i) = 1; mu(i) = y(i);}
//        else{ind(i) = 0; mu(i) = mu0(i);}            
//      }
//      else{ind(i) = 0; mu(i) = mu0(i);}
//    }      
//    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
//    // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
//    else
//    {
//      ind(i) = 0; mu(i) = mu0(i);
//      Rcpp::Rcout << "Something went wrong when updating mu " << i+1 << std::endl;
//      Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");        
//    }
//  }
  
  // OUTPUT
  return join_rows(mu, ind);
}

/* Metropolis-Hastings updates of mu 
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdateNoSpikesConstrainSequential(
    arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q_0)'$ */
    arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
    double const& Constrain,
    arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */  
    arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
    arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
    arma::vec const& sum_bycell_all, /* Sum of expression counts by cell (biological genes only) */
    double const& s2_mu,
    int const& q_bio,
    int const& n,
    arma::vec & mu,
    arma::vec & ind,
    int const& ref,
    arma::vec const& Index)
{
  using arma::span;
  
  arma::uvec IndexAux = find(Index != ref);
  
  // PROPOSAL STEP    
//  arma::vec y = exp(SigmaAuxChol*arma::randn(q_bio-1) + log(mu0(span(0, q_bio-2))));
//  arma::vec y = exp(arma::randn(q_bio-1) % sqrt(prop_var(span(0, q_bio-2))) + log(mu0(span(0, q_bio-2))));
  arma::vec y = exp(arma::randn(q_bio) % sqrt(prop_var) + log(mu0));
  arma::vec u = arma::randu(q_bio);
  // INITIALIZE MU
//  Rcpp::Rcout << "mu0(0) before mu" << mu0(0) << std::endl; 
  mu = mu0 + 1 - 1;
//  Rcpp::Rcout << "mu0(0) after mu" << mu0(0) << std::endl;
  double aux;
  double aux2; 
  double sumAux = sum(log(mu0(IndexAux)));

  // ACCEPT/REJECT STEP
//  arma::vec log_aux = (log(y) - log(mu0(span(0, q_bio-2)))) % sum_bycell_all(span(0, q_bio-2));
  arma::vec log_aux = (log(y) - log(mu0)) % sum_bycell_all;
  
  // Independent prior 
//  log_aux -= (0.5/s2_mu) * (pow(log(y),2) - pow(log(mu0(span(0, q_bio - 2))),2));
  // COMPUTING THE LIKELIHOOD CONTRIBUTION OF THE ACCEPTANCE RATE (NO NEED TO BE SEQUENTIAL)
//  for (int i=0; i < q_bio - 1; i++)
  for (int i=0; i < q_bio; i++)
  {
    if(i != ref)
    {
      for (int j=0; j < n; j++) 
      {
        log_aux(i) -= ( Counts(i,j) + (1/delta(i)) ) * 
          log( ( nu(j)*y(i)+(1/delta(i)) ) / ( nu(j)*mu(i)+(1/delta(i)) ));
      }
//      aux = 0.5 * (q_bio * Constrain - (sum(log(mu.elem(IndexAux))) - log(mu(i))));
      aux = 0.5 * (q_bio * Constrain - (sumAux - log(mu(i))));
      //    aux = 0.5 * (log(mu(i)) + log(mu(ref)));
//      Rcpp::Rcout << "aux" << aux << std::endl;
//      Rcpp::Rcout << "aux2" << aux2 << std::endl;
      //    log_aux(i) -= (0.5 * q_bio /s2_mu) * (pow(log(y(i)) - aux,2)); 
      //    log_aux(i) += (0.5 * q_bio /s2_mu) * (pow(log(mu0(i)) - aux,2));   
      log_aux(i) -= (0.5 * 2 /s2_mu) * (pow(log(y(i)) - aux,2)); 
      log_aux(i) += (0.5 * 2 /s2_mu) * (pow(log(mu0(i)) - aux,2));   
      if(log(u(i)) < log_aux(i) & y(i) > 1e-3) 
      {
        ind(i) = 1; mu(i) = y(i);
        sumAux += log(mu(i)) - log(mu0(i)); 
      }
      else{ind(i) = 0; mu(i) = mu0(i); }
      // FINAL GENE
      //    mu(q_bio-1) = exp(q_bio * Constrain - sum(log(mu(span(0, q_bio-2)))));
      //    mu(ref) = exp(q_bio * Constrain - sum(log(mu.elem(IndexAux))));
    }
  }
  // FINAL GENE
  ind(ref) = 0;
  //    mu(q_bio-1) = exp(q_bio * Constrain - sum(log(mu(span(0, q_bio-2)))));  
//  mu(ref) = exp(q_bio * Constrain - sum(log(mu.elem(IndexAux))));
  mu(ref) = exp(q_bio * Constrain - sumAux);
  Rcpp::Rcout << "sumAux: " << sumAux << std::endl;
//  Rcpp::Rcout << "mu(ref): " << mu(ref) << std::endl;
//  Rcpp::Rcout << "mu(ref) alt:" << exp(q_bio * Constrain - sumAux) << std::endl;
  
//  Rcpp::Rcout << "mu0(0) after update" << mu0(0) << std::endl;
    
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
  int const& q_bio,
  int const& n,
  double const& s2_delta,
  double const& prior_delta)
{
  using arma::span;
  
  // CREATING VARIABLES WHERE TO STORE DRAWS
  arma::vec delta = arma::zeros(q_bio); 
  arma::vec ind = arma::zeros(q_bio);
  
  // PROPOSAL STEP
  arma::vec y = exp(arma::randn(q_bio) % sqrt(prop_var) + log(delta0));
  arma::vec u = arma::randu(q_bio);
  
  // ACCEPT/REJECT STEP 
  arma::vec log_aux = - n * (lgamma_cpp(1/y)-lgamma_cpp(1/delta0));
  
  // Loop to replace matrix operations, through genes and cells
  for (int i=0; i < q_bio; i++)
  {
    for (int j=0; j < n; j++) 
    {
      log_aux(i) += R::lgammafn(Counts(i,j) + (1/y(i))) - R::lgammafn(Counts(i,j) + (1/delta0(i)));
      log_aux(i) -= ( Counts(i,j) + (1/y(i)) ) *  log( nu(j)*mu(i)+(1/y(i)) );
      log_aux(i) += ( Counts(i,j) + (1/delta0(i)) ) *  log( nu(j)*mu(i)+(1/delta0(i)) );
    }
  }
  
  // +1 should appear because we update log(delta) not delta. However, it cancels out with the prior. 
  log_aux -= n * ( (log(y)/y) - (log(delta0)/delta0) );
  // Component related to the prior
  if(prior_delta == 1) {log_aux += (log(y) - log(delta0)) * a_delta - b_delta * (y - delta0);}
  else {log_aux -= (0.5/s2_delta) * (pow(log(y),2) - pow(log(delta0),2));}
  
  // CREATING OUTPUT VARIABLE & DEBUG
  for (int i=0; i < q_bio; i++)
  {
    if(arma::is_finite(log_aux(i)))
    {
      if(log(u(i)) < log_aux(i))
      {
        // DEBUG: Reject very small values to avoid numerical issues
        if(arma::is_finite(y(i)) & (y(i) > 1e-3)) {ind(i) = 1; delta(i) = y(i);}
        else{ind(i) = 0; delta(i) = delta0(i);}            
      }
      else{ind(i) = 0; delta(i) = delta0(i);}
    }      
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
    else
    {
      ind(i) = 0; delta(i) = delta0(i);
      Rcpp::Rcout << "Something went wrong when updating delta " << i+1 << std::endl;
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
  int const& q_bio,
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
  
  // CODE CRUSHES WHEN USING THIS LOOP, CHECK WHY
  //    arma::vec log_aux = arma::zeros(n);
  
  // Loop to replace matrix operations, through genes and cells
  //    for (int j=0; j < n; j++) 
  //    {
  //      for (int i=0; i < q_bio; i++) 
  //      {
  //        log_aux(j) -= ( Counts(i,j) + (1/delta(i)) ) *  
  //                      log( ( phi(j)*exp(y(j))*mu(i)+(1/delta(i)) ) / ( phi(j)*nu0(j)*mu(i)+(1/delta(i)) ));
  //      } 
  //      log_aux(j) += (y(j) - lognu(j)) * (sum_bygene_all(j) + 1/thetaBatch(j)); // - ( y(j)-nu0(j) )  * (SumSpikeInput + (1/(thetaBatch(j)*s(j))));
  //      log_aux(j) -= nu0(j) * (exp(y(j)-lognu(j))-1)  * (SumSpikeInput + (1/(thetaBatch(j) * s(j))));
  //    }
  
  
  arma::vec log_aux = (y - lognu) % (sum_bygene_all + 1/thetaBatch);
  log_aux -= nu0 % (exp(y-lognu)-1)  % (1/(thetaBatch % phi));   
  arma::mat m = Counts.rows(0, q_bio - 1);
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
    Rcpp::Rcout << "Something went wrong when updating nu" << size(find_nonfinite(log_aux),0) << std::endl;
    Rcpp::stop("Please consider additional filter of the input dataset.");
  }
  
  // CREATING OUTPUT VARIABLE        
  arma::vec nu = ind % exp(y) + (1 - ind) % nu0;
  
  // QUICK FIX FOR OFFSET ISSUE (NO FORMAL JUSTIFICATION)
//  for (int k=0; k < nBatch; k++)
//  {
//    nu.elem(find(BatchInfo == BatchIds(k))) = BatchSizes(k) * BatchOffSet(k) * nu.elem(find(BatchInfo == BatchIds(k))) / sum(nu.elem(find(BatchInfo == BatchIds(k))));
//  }
  
  // OUTPUT
  return join_rows(nu, arma::conv_to<arma::mat>::from(ind));
}

/* Metropolis-Hastings updates of phi 
 * Joint updates using Dirichlet proposals
 */
Rcpp::List phiUpdateNoSpikes(
    arma::vec const& phi0, // Current value of $\phi=(\phi_1,...,\phi_n)'$
    arma::vec const& prop_var, // Current value of the proposal precision
    arma::vec const& nu, // Current value of $\nu=(\nu_1,...,\nu_n)'$
    arma::vec const& theta, 
    arma::vec const& p_phi, // Dirichlet hyper-parameter of the prior for $\phi / n$
    int const& n, // Total number of cells 
    arma::vec const& BatchInfo,
    arma::vec const& BatchIds,
    int const& nBatch,
    arma::vec const& BatchSizes,
    arma::vec const& BatchOffSet)
{
  using arma::span;
  using arma::pow;
 
  arma::vec phi = arma::ones(n);
  arma::vec y; 
  double u;
  arma::vec log_aux = arma::ones(nBatch);
  arma::vec ind = arma::zeros(nBatch);
  
  // LOOP OVER BATCHES
  for (int k=0; k < nBatch; k++)
  {
    // PROPOSAL STEP
    y = BatchOffSet(k) * BatchSizes(k) * rDirichlet(prop_var(k) * phi0.elem(find(BatchInfo == BatchIds(k))));
    u = R::runif(0,1);
    
//    Rcpp::Rcout << "Proposal"  << std::endl;
//    Rcpp::Rcout << y << std::endl;
//    Rcpp::Rcout << "Old value"  << std::endl;
//    Rcpp::Rcout << phi0.elem(find(BatchInfo == BatchIds(k))) << std::endl;
   if(all(prop_var(k) * y < 2.5327372760800758e+305)  & 
      all(prop_var(k) *  phi0.elem(find(BatchInfo == BatchIds(k))) < 2.5327372760800758e+305) &
      all(y > 0) &
      all( phi0.elem(find(BatchInfo == BatchIds(k))) > 0 )) 
   {
     // ACCEPTANCE STEP 
     log_aux(k) = prop_var(k) * sum(y % log(phi0.elem(find(BatchInfo == BatchIds(k)))) - phi0.elem(find(BatchInfo == BatchIds(k))) % log(y));
     log_aux(k) += sum((p_phi.elem(find(BatchInfo == BatchIds(k))) - 1/theta(k)) % (log(y) - log(phi0.elem(find(BatchInfo == BatchIds(k))))));
     log_aux(k) -= (1/theta(k)) * sum(nu.elem(find(BatchInfo == BatchIds(k))) % ( (1/y) - (1/phi0.elem(find(BatchInfo == BatchIds(k)))) ) );
     log_aux(k) -= sum(lgamma_cpp_vec(prop_var(k) * y) - lgamma_cpp_vec(prop_var(k) * phi0.elem(find(BatchInfo == BatchIds(k)))));

     if(!R_IsNA(log_aux(k)))
     {
       if(log(u) < log_aux(k)) {ind(k) = 1;}
       else {ind(k) = 0;}
     }
     // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
     // DEBUG: Print warning message
     else
     {
        Rcpp::Rcout << "Something went wrong when updating phi" << std::endl;
        Rcpp::stop("Please consider additional filter of the input dataset."); 
        ind(k) = 0;
     }      
   }
   else
   {
     ind(k) = 0;
   }
   
   phi.elem(find(BatchInfo == BatchIds(k))) = ind(k) * y + (1-ind(k)) * phi0.elem(find(BatchInfo == BatchIds(k)));
  }
  
//  Rcpp::Rcout << "Sampling phi"  << std::endl;
//  Rcpp::Rcout << phi << std::endl;
  // ACCEPT/REJECT STEP (REJECT VALUES OUTSIDE VALID RANGE)  
//  if(all(prop_var * y < 2.5327372760800758e+305)  & 
//     all(prop_var * phi0 < 2.5327372760800758e+305) &
//     all(y > 0) &
//     all(phi0 > 0)) 
//  {
//    double log_aux = sum( (sum_bygene_all + p_phi) % (log(y) - log(phi0)));
    
    // Loop to replace matrix operations, through genes and cells
//    for (int j=0; j < n; j++) 
//    {
//      for (int i=0; i < q_bio; i++) 
//      {
//        log_aux -= ( Counts(i,j) + (1/delta(i)) ) *  
//          log( (y(j)*nu(j)*mu(i)+(1/delta(i)) ) / ( phi0(j)*nu(j)*mu(i)+(1/delta(i)) ));
//      } 
//    }
//    log_aux += prop_var * sum(y % log(phi0) - phi0 % log(y));
//    log_aux -= sum(lgamma_cpp_vec(prop_var * y) - lgamma_cpp_vec(prop_var * phi0));    
    
//    if(!R_IsNA(log_aux))
//    {
//      if(log(u) < log_aux) {ind = 1;}
//      else {ind = 0;}
//      phi = ind * y + (1-ind) * phi0;
//    }
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
//    else
//    {
//      Rcpp::Rcout << "Something went wrong when updating phi" << std::endl;
//      Rcpp::stop("Please consider additional filter of the input dataset."); 
//      ind = 0;
//      phi = phi0;
//    }     
//  }
//  else
//  {
//    ind = 0;
//    phi = phi0;     
//  }      
  return(Rcpp::List::create(
      Rcpp::Named("phi")=phi,
      Rcpp::Named("ind")=ind)); 
}

/* MCMC sampler 
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcppNoSpikes(
    int N, // Total number of MCMC draws 
    int thin, // Thinning period for MCMC chain 
    int burn, // Burning period for MCMC chain 
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
    NumericVector p_Phi, // Dirichlet hyper-parameter for $\phi / n$
    double aphi,
    double bphi,
    double atheta, // Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
    double btheta, // Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
    double ar, // Optimal acceptance rate for adaptive Metropolis-Hastings updates
    NumericVector LSmu0, // Starting value of adaptive proposal variance of $\mu=(\mu_1,...,\mu_q_0)'$ (log-scale)
    NumericVector LSdelta0, // Starting value of adaptive proposal variance of $\delta=(\delta_1,...,\delta_{q_0})'$ (log-scale)
    NumericVector LSphi0, // Starting value of adaptive proposal precision of $\phi=(\phi_1,...,\phi_n)'$ (log-scale)
    NumericVector LSnu0, // Starting value of adaptive proposal variance of $\nu=(\nu_1,...,\nu_n)'$ (log-scale)
    double LStheta0, // Starting value of adaptive proposal variance of $\theta$ (log-scale)  
    NumericVector sumByCellAll, // Sum of expression counts by cell (all genes)
    NumericVector sumByGeneAll, // Sum of expression counts by gene (all genes)
    int StoreAdapt, 
    int EndAdapt,
    int PrintProgress,
    double s2_delta, 
    double prior_delta,
    NumericVector BatchInfo,
    NumericVector BatchIds,
    NumericVector BatchSizes,
    NumericVector BatchOffSet,
    double Constrain,
    NumericMatrix InvCovMu,
    NumericVector Index,
    int ref)
{

   //   NumericMatrix CholCovMu,
    //  ) 
  
  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int qbio = Counts.nrow(); int Naux = N/thin - burn/thin;
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
//  arma::mat CholCovMu_arma = as_arma(CholCovMu);
  arma::mat InvCovMu_arma = as_arma(InvCovMu);
  arma::vec Index_arma = as_arma(Index);
  
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = arma::zeros(Naux, qbio); 
  arma::mat delta = arma::zeros(Naux, qbio); 
  arma::mat phi = arma::ones(Naux, n);
  arma::mat nu = arma::zeros(Naux, n); 
  arma::mat theta = arma::zeros(Naux, nBatch); 
  arma::mat LSmu;
  arma::mat LSdelta;
//  arma::mat LSphi;
  arma::mat LSnu;
  arma::mat LStheta; 
  
  // LOG-PROPOSAL VARIANCES 
  if(StoreAdapt == 1)
  {
    LSmu = arma::zeros(Naux, qbio); 
    LSdelta = arma::zeros(Naux, qbio); 
//    LSphi = arma::ones(Naux, nBatch);   
    LSnu = arma::zeros(Naux, n); 
    LStheta = arma::zeros(Naux, nBatch);   
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(qbio); arma::vec PmuAux = arma::zeros(qbio);
  arma::vec deltaAccept = arma::zeros(qbio); arma::vec PdeltaAux = arma::zeros(qbio);
//  arma::vec phiAccept = arma::zeros(nBatch);  arma::vec PphiAux = arma::zeros(nBatch);
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  arma::vec thetaAccept = arma::zeros(nBatch); arma::vec PthetaAux = arma::zeros(nBatch);
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(qbio,2); muAux.col(0)=as_arma(mu0); arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(qbio,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::vec phiAux = phi0_arma; //arma::vec LSphiAux = LSphi0; Rcpp::List phiAuxList;
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::mat thetaAux = arma::zeros(nBatch, 2); thetaAux.col(0) = theta0 * arma::ones(nBatch); 
  arma::vec LSthetaAux = LStheta0 * arma::ones(nBatch);  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(qbio); arma::vec PdeltaAux0 = arma::zeros(qbio);
//  arma::vec PphiAux0 = arma::zeros(nBatch); 
  arma::vec PnuAux0 = arma::zeros(n); arma::vec PthetaAux0 = arma::ones(nBatch);
  
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;
  
  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  arma::vec muUpdateAux = arma::ones(qbio);
  arma::vec indQ = arma::zeros(qbio);
  double LSmuAuxExtra = 0;
//  int ref = qbio-1;
  arma::vec RefFreq = arma::zeros(qbio); 
//  arma::mat D; arma::mat SigmaAux; arma::mat SigmaAuxChol;
//  arma::mat OnesMat = arma::ones(qbio-1, qbio-1);
//  D = diagmat(exp(LSmuAux(arma::span(0, qbio-2))));
//  SigmaAux = D - (1/sum(exp(LSmuAux(arma::span(0, qbio-2))))) * D * OnesMat * D;
//  SigmaAuxChol = arma::chol(SigmaAux);
  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;  
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  
  // START OF MCMC LOOP
  for (i=0; i<N; i++) {
    
    Rcpp::checkUserInterrupt();
    
    if(i==burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
      Rcpp::Rcout << "End of burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl; 
    }
    
    Ibatch++; 
    
    //    struct timespec time0_1 = orwl_gettime();
    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    // WE CAN RECYCLE THE SAME FULL CONDITIONAL AS IMPLEMENTED FOR S (BATCH CASE)
    
    phiAux = sUpdateBatch(phiAux, nuAux.col(0), thetaAux.col(0),
                          aphi, bphi, BatchDesign_arma, n);
    //phiAuxList = phiUpdateNoSpikes(phiAux, exp(LSphiAux), 
    //                       nuAux.col(0), thetaAux.col(0), 
    //                       p_Phi, n,
    //                       BatchInfo_arma, BatchIds_arma, nBatch, BatchSizes_arma, BatchOffSet_arma); 
    //phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    //PphiAux += as<arma::vec>(phiAuxList["ind"]); if(i>=burn) {phiAccept += as<arma::vec>(phiAuxList["ind"]);}
    //    struct timespec time1_1 = orwl_gettime();    
    //    Rcpp::Rcout << "Time phi: " << (time1_1.tv_nsec - time0_1.tv_nsec) / ((float)(n)) << std::endl;
    
    //    struct timespec time0_2 = orwl_gettime();
    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdateBatch(thetaAux.col(0), exp(LSthetaAux), BatchDesign_arma,
                                phiAux, nuAux.col(0), atheta, btheta, n, nBatch);
    PthetaAux += thetaAux.col(1); if(i>=burn) {thetaAccept += thetaAux.col(1);}
    //    struct timespec time1_2 = orwl_gettime();    
    //    Rcpp::Rcout << "Time theta: "  <<  (time1_2.tv_nsec - time0_2.tv_nsec) / ((float)(nBatch)) << std::endl;
    
    //    struct timespec time0_3 = orwl_gettime();
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR 

//    if(i > 2000)
//    {
//      ref = as_scalar(arma::randi( 1, arma::distr_param(0,qbio-1) ));
//      Rcpp::Rcout << "ref: " << ref << std::endl; 
//      muAux = muUpdateNoSpikesConstrain(muAux.col(0), exp(LSmuAuxExtra + LSmuAux), Constrain, Counts_arma, deltaAux.col(0), 
//                                        nuAux.col(0), sumByCellAll_arma, s2mu, qbio, n,
//                                        muUpdateAux, indQ, InvCovMu_arma, mu0_arma, ref, Index_arma); // CholCovMu_arma,
//      PmuAux += muAux.col(1); if(i>=burn) {muAccept(0) += muAux(0,1);}
//    }
//    else
//    {
//      ref = as_scalar(arma::randi( 1, arma::distr_param(0,qbio-1) )); 
      Rcpp::Rcout << "ref: " << ref << std::endl;
      RefFreq(ref) += 1;
      muAux = muUpdateNoSpikesConstrainSequential(muAux.col(0), exp(LSmuAux), Constrain, Counts_arma, deltaAux.col(0), 
                                                  nuAux.col(0), sumByCellAll_arma, s2mu, qbio, n,
                                                  muUpdateAux, indQ, ref, Index_arma);
      PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}
//    }
//    muAux = muUpdateNoSpikes(muAux.col(0), exp(LSmuAux), Constrain, Counts_arma, deltaAux.col(0), 
//                             nuAux.col(0), sumByCellAll_arma, s2mu, qbio, n,
//                             muUpdateAux, indQ); 
//    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}  
    //    struct timespec time1_3 = orwl_gettime();    
    //    Rcpp::Rcout << "Time mu: "  <<  (time1_3.tv_nsec - time0_3.tv_nsec) / ((float)(qbio)) << std::endl;
    
    //    struct timespec time0_5 = orwl_gettime();
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdateNoSpikes(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, 
                           muAux.col(0), nuAux.col(0), adelta, bdelta, qbio, n, s2_delta, prior_delta);  
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);} 
    //    struct timespec time1_5 = orwl_gettime();    
    //    Rcpp::Rcout << "Time delta: "  << (time1_5.tv_nsec - time0_5.tv_nsec) / ((float)(qbio)) << std::endl;
    
    //    struct timespec time0_6 = orwl_gettime();
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateNoSpikes(nuAux.col(0), exp(LSnuAux), Counts_arma, 
                            BatchDesign_arma,
                            muAux.col(0), deltaAux.col(0),
                            phiAux, thetaAux.col(0), sumByGeneAll_arma, qbio, n,
                            BatchInfo_arma, BatchIds_arma, nBatch,
                            BatchSizes_arma, BatchOffSet_arma); 
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}
    //    struct timespec time1_6 = orwl_gettime();    
    //    Rcpp::Rcout << "Time nu: " <<  (time1_6.tv_nsec - time0_6.tv_nsec) / ((float)(n)) << std::endl;
    
    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
//        if(i > 2000)
//        {
//          PmuAux=PmuAux/50; PmuAux(0) = -1+2*(PmuAux(0)>ar);
//          LSmuAuxExtra = LSmuAuxExtra + PmuAux(0)*0.1;     
//        }
//        else
//        {
          PmuAux=PmuAux/(50-RefFreq); PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);//
          LSmuAux.elem(find(Index_arma != ref)) = LSmuAux.elem(find(Index_arma != ref)) + PmuAux.elem(find(Index_arma != ref))*0.1; 
//        }
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux=LSdeltaAux+PdeltaAux*0.1;                
//        PphiAux=PphiAux/50; PphiAux = -1+2*arma::conv_to<arma::mat>::from(PphiAux>ar);//-1+2*(PphiAux(0)>ar); 
//        LSphiAux = LSphiAux - PphiAux*0.1;  
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux=LSnuAux+PnuAux*0.1; 
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*arma::conv_to<arma::mat>::from(PthetaAux>ar); 
        LSthetaAux= LSthetaAux + PthetaAux*0.1;
        
        Ibatch = 0; 
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0; 
//        PphiAux = PphiAux0; 
        PnuAux = PnuAux0; PthetaAux = PthetaAux0;
        RefFreq = arma::zeros(qbio);
      }
      
    }
    
    // STORAGE OF DRAWS
    if(i%thin==0 & i>=burn)
    {      
      mu.row(i/thin - burn/thin) = muAux.col(0).t(); 
      delta.row(i/thin - burn/thin) = deltaAux.col(0).t(); 
      phi.row(i/thin - burn/thin) = phiAux.t();
      nu.row(i/thin - burn/thin) = nuAux.col(0).t();       
      theta.row(i/thin - burn/thin) = thetaAux.col(0).t();   
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/thin - burn/thin) = LSmuAux.t();
        LSdelta.row(i/thin - burn/thin) = LSdeltaAux.t();
//        LSphi.row(i/thin - burn/thin) = LSphiAux.t();
        LSnu.row(i/thin - burn/thin) = LSnuAux.t();
        LStheta.row(i/thin - burn/thin) = LSthetaAux.t(); 
      }
    }
    
    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if(i%(2*thin) == 0 & PrintProgress == 1)
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
//      Rcpp::Rcout << "LSphi: " << LSphiAux << std::endl;
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

  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << muAccept(0)/(N-burn) << std::endl;
  
  Rcpp::Rcout << "Minimum acceptance rate among mu[i]'s: " << min(muAccept(arma::span(0, qbio - 2))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << mean(muAccept(arma::span(0, qbio - 2))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among mu[i]'s: " << max(muAccept(arma::span(0, qbio - 2))/(N-burn)) << std::endl;
  
//  Rcpp::Rcout << "LS mu[i]'s: " << LSmuAuxExtra + LSmuAux << std::endl;
  
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among delta[i]'s: " << min(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << mean(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << max(deltaAccept/(N-burn)) << std::endl;
//  Rcpp::Rcout << " " << std::endl;
//  Rcpp::Rcout << "Acceptance rate for phi (joint): " << phiAccept/(N-burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among nu[jk]'s: " << min(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among nu[jk]'s: " << mean(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among nu[jk]'s: " << max(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among theta[k]'s: " << min(thetaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among theta[k]'s: " << mean(thetaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among theta[k]'s: " << max(thetaAccept/(N-burn)) << std::endl;  
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  
  if(StoreAdapt == 1)
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
        Rcpp::Named("mu")=mu,
        Rcpp::Named("delta")=delta,
        Rcpp::Named("phi")=phi,
        Rcpp::Named("nu")=nu,
        Rcpp::Named("theta")=theta,
        Rcpp::Named("ls.mu")=LSmu,
        Rcpp::Named("ls.delta")=LSdelta,
//        Rcpp::Named("ls.phi")=LSphi,
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
        Rcpp::Named("nu")=nu,
        Rcpp::Named("theta")=theta)); 
  }
  
}

