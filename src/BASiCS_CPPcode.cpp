/* C++ implementation of BASiCS */

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

  y = norm_rand();
  y *= y;
  mu2 = mu * mu;
  l2 = 2.0*lambda;
  x1 = mu + mu2*y/l2 - (mu/l2)* sqrt(4.0*mu*lambda*y + mu2*y*y);

  u = unif_rand();
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

/* Auxiliary function: to avoid numerical overflows/underflows when computing a sum of exponentiated values */
double log_sum_exp_cpp(arma::vec const& x_arma) 
{
  double offset;
  if ( max(abs(x_arma)) > max(x_arma) ) { offset = min(x_arma);}
  else { offset = max(x_arma);}
  return log(sum(exp(x_arma - offset))) + offset; 
}

/* Auxiliary function required for some of the Metropolis-Hastings updates of mu */
arma::mat UpdateAux_muTrick(
  arma::vec const& aux1mu, /* Auxiliary variable (function of the current and proposed value of mu) */
  arma::vec const& mu_bio, /* Current value of $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu) /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
{
  arma::mat x = ((phi % nu) * aux1mu.t()).t();
  x.each_col() += 1 / (mu_bio % delta);
  return x;
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
  int const& q_bio)
{
    using arma::span;

    // CREATING VARIABLES WHERE TO STORE DRAWS
    arma::vec logmu = log(mu0(span(0, q_bio - 1)));

    // PROPOSAL STEP    
    arma::vec y = arma::randn(q_bio) % sqrt(prop_var(span(0,q_bio - 1))) + logmu;
    arma::vec u = arma::randu(q_bio);
    
    // ACCEPT/REJECT STEP
    arma::mat m = Counts.rows(0, q_bio - 1);
    m.each_col() += 1 / delta;
    arma::mat num = UpdateAux_muTrick(exp(y-logmu), mu0(span(0, q_bio - 1)), delta, phi, nu );
    num /= UpdateAux_muTrick(arma::ones(q_bio), mu0(span(0, q_bio - 1)), delta, phi, nu );
    m %= log(num);   
    arma::vec log_aux = (y - logmu) % sum_bycell_bio - sum(m, 1); 
    log_aux -= (0.5/s2_mu) * (pow(y,2) - pow(logmu,2));
    arma::umat ind = log(u) < log_aux;
    ind = join_cols(ind, arma::uvec(q - q_bio, arma::fill::zeros));
    arma::vec aux = join_cols(exp(y), arma::vec(q - q_bio, arma::fill::zeros));
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
    ind.elem(find_nonfinite(log_aux)).zeros();
    if(size(find_nonfinite(log_aux),0)>0) 
    {
      Rcpp::Rcout << "Something went wrong when updating mu" << size(find_nonfinite(log_aux),0) << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset.");
    }
    // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
    // DEBUG: Print warning message
    ind.elem(find_nonfinite(aux)).zeros();
    if(size(find_nonfinite(aux),0)>0) 
    {
      Rcpp::Rcout << "Something went wrong when updating mu" << size(find_nonfinite(aux),0) << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset.");
    }

    // CREATING OUTPUT VARIABLE
    arma::vec mu = ind % aux + (1 - ind) % mu0;
    
    // OUTPUT
    return join_rows(mu, arma::conv_to<arma::mat>::from(ind));
}

/* Auxiliary function required for some of the Metropolis-Hastings updates of mu */
arma::mat UpdateAux_deltaTrick(
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu) /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
{
  arma::mat x = ((phi % nu) * mu.t()).t();
  x.each_col() += 1 / delta;
  return x;
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
  int const& n)
{
    using arma::span;
    
    // CREATING VARIABLES WHERE TO STORE DRAWS
    arma::vec logdelta = log(delta0);

    // PROPOSAL STEP
    arma::vec y = arma::randn(q_bio) % sqrt(prop_var) + logdelta;
    arma::vec u = arma::randu(q_bio);
    
    // ACCEPT/REJECT STEP
    arma::vec log_aux = (y - logdelta) * a_delta - n * (logdelta/delta0) % ((y/logdelta) % exp(-y+logdelta) - 1);    
    arma::mat m1 = Counts.rows(0, q_bio - 1);
    m1.each_col() += exp(-y);
    arma::mat m2 = Counts.rows(0, q_bio - 1);
    m2.each_col() += 1/delta0;
    log_aux += sum(lgamma_cpp(m1)-lgamma_cpp(m2),1) - n*(lgamma_cpp(exp(-y))-lgamma_cpp(1/delta0));
    arma::mat num = UpdateAux_deltaTrick(mu(span(0, q_bio - 1)), exp(y), phi, nu);
    arma::mat den = UpdateAux_deltaTrick(mu(span(0, q_bio - 1)), delta0, phi, nu);
    m1 %= log(num);
    m2 %= log(den);
    log_aux -= sum(m1-m2,1) + b_delta* delta0 % (exp(y-logdelta)-1);
    arma::umat ind = log(u) < log_aux;
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
    ind.elem(find_nonfinite(log_aux)).zeros();
    if(size(find_nonfinite(log_aux),0)>0)
    {
      Rcpp::Rcout << "Something went wrong when updating delta" << size(find_nonfinite(log_aux),0) << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset.");
    }

    // CREATING OUTPUT VARIABLE
    arma::vec delta = ind % exp(y) + (1 - ind) % delta0;
    
    // OUTPUT
    return join_rows(delta, arma::conv_to<arma::mat>::from(ind));
}

/* Draws for cell-specific normalising constants s[j]
 * Metropolis-Hastings updates are not required as full conditionals have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::vec sUpdate(
  arma::vec const& s0_arma, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu_arma, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& theta, /* Current value of $\theta$ */
  double const& as, /* Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double const& bs, /* Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  int const& n)
{
   // Allowing the use of rgig within C++ code (transform arma::vec elements into NumericVector elements)
   NumericVector s0 = wrap(s0_arma);
   NumericVector nu = wrap(nu_arma);
   
   // Creating a vector where draws will be stored
   arma::vec aux = arma::zeros(n);
   
   // Calculating parameters to the passed as input to the Rgig function (common for all cells)
   double p = as-1/theta; 
   double b = 2*bs;
   
   // GIG draws
   /* DEBUG: return original value of s0 if input values are not within the appropriate range (to avoid problems related to numerical innacuracy) */
   if(!R_IsNA(p)) 
   {
     // Calculating parameter to the passed as input to the Rgig function (specific to each cell)
      NumericVector a = 2*nu/theta;      
      int j;
      for (j=0; j<n; j++) 
      {
          if(!R_IsNA(a(j)) & (a(j)>0)) 
          {
            aux(j) = Rcpp::as<double>(Rgig(1,p,a(j),b));
            /* DEBUG: break in case of undefined values */
            if(R_IsNA(aux(j))) 
            {
              Rcpp::Rcout << "Something went wrong when updating s" << j << std::endl;
              Rcpp::stop("Please consider additional filter of the input dataset.");
            }
          }
      else if(!(a(j)<0) & (p>0)) {aux(j) = Rcpp::as<double>(Rgig(1,p,a(j),b));}
      else{aux(j)=s0_arma(j);}         
      }
      return aux;     
   }
  else return s0_arma;
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
    arma::vec lognu = log(nu0);

    // PROPOSAL STEP    
    arma::vec y = arma::randn(n) % sqrt(prop_var) + lognu;
    arma::vec u = arma::randu(n);
       
    // ACCEPT/REJECT STEP
    arma::vec log_aux = (y - lognu) % (sum_bygene_all + 1/theta);
    log_aux -= nu0 % (exp(y-lognu)-1)  % (sum(mu(span(q_bio,q -1))) + (1/(theta*s)));    
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
  double u = as<double>(wrap(arma::randu(1)));;
  
  arma::vec phi; 
  bool ind;
    
  // ACCEPT/REJECT STEP (REJECT VALUES OUTSIDE VALID RANGE)  
  if(all(prop_var * y < 2.5327372760800758e+305)  & 
     all(prop_var * phi0 < 2.5327372760800758e+305) &
     all(y > 0) &
     all(phi0 > 0)) 
  {
    arma::mat m = Counts.rows(0, q_bio - 1);
    m.each_col() += 1/delta;
    arma::mat log_num = UpdateAux_phiTrick(arma::ones(q_bio), 
                                           delta % mu(span(0, q_bio - 1)), 
                                           y % nu);
    arma::mat log_den = UpdateAux_phiTrick(arma::ones(q_bio), 
                                           delta % mu(span(0, q_bio - 1)), 
                                           phi0 % nu);  
    m %= (log_num - log_den);
  
    double log_aux1 = sum( (sum_bygene_bio + p_phi) % (log(y) - log(phi0))) - sum(sum(m, 0));
    double log_aux2 = prop_var * sum(y % log(phi0) - phi0 % log(y));
    double log_aux3 = as<double>(wrap(sum(lgamma_cpp(prop_var * y) - lgamma_cpp(prop_var * phi0))));
    double log_aux = log_aux1 + log_aux2 - log_aux3;

    if(!R_IsNA(log_aux))
    { 
      ind = log(u) < log_aux;
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
         Rcpp::Named("ind")=Rcpp::as<double>(wrap(ind)))); 
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
  int PrintProgress) 
{

// NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int q = Counts.nrow(); int qbio = delta0.size(); int Naux = N/thin - burn/thin;
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts);
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
  arma::mat muAux = arma::zeros(q,2); muAux.col(0)=as_arma(mu0); arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(qbio,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::vec phiAux = as_arma(phi0); double LSphiAux = LSphi0; Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0); 
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::vec thetaAux = arma::zeros(2); thetaAux(0) = theta0; double LSthetaAux = LStheta0;  
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q); arma::vec PdeltaAux0 = arma::zeros(qbio);
  double PphiAux0 = 0; 
  arma::vec PnuAux0 = arma::zeros(n); double PthetaAux0=0;
  
  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;
  
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
    
    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, muAux.col(0), deltaAux.col(0),
                           nuAux.col(0), p_Phi, sumByGeneBio_arma, qbio, n); 
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]); if(i>=burn) {phiAccept += Rcpp::as<double>(phiAuxList["ind"]);}

    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdate(thetaAux(0), exp(LSthetaAux), sAux, nuAux.col(0), atheta, btheta, n);
    PthetaAux += thetaAux(1); if(i>=burn) {thetaAccept += thetaAux(1);}
    
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR    
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, deltaAux.col(0), 
                     phiAux, nuAux.col(0), sumByCellBio_arma, s2mu, q, qbio);     
    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}
    
    // UPDATE OF S
    sAux = sUpdate(sAux, nuAux.col(0), thetaAux(0), as, bs, n); 
    
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, muAux.col(0), 
                           nuAux.col(0), phiAux, adelta, bdelta, qbio, n);  
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);}    
    
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdate(nuAux.col(0), exp(LSnuAux), Counts_arma, 
                     muAux.col(0), deltaAux.col(0), phiAux, sAux, thetaAux(0), sumByGeneAll_arma, q, qbio, n); 
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/50; PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux=LSmuAux+PmuAux*std::min(0.01,1/sqrt(i)); 
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.01,1/sqrt(i));                 
        PphiAux=PphiAux/50; PphiAux = -1+2*(PphiAux>ar); 
        LSphiAux=LSphiAux - PphiAux*std::min(0.01,1/sqrt(i));  
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux=LSnuAux+PnuAux*std::min(0.01,1/sqrt(i)); 
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*(PthetaAux>ar); 
        LSthetaAux=LSthetaAux+PthetaAux*std::min(0.01,1/sqrt(i));
      
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
   NumericVector s0 = wrap(s0_arma);
   NumericVector nu = wrap(nu_arma);
   NumericVector thetaBatch = wrap(BatchDesign * theta); 
   
   // Creating a vector where draws will be stored
   arma::vec aux = arma::zeros(n);
   
   // Calculating parameters to the passed as input to the Rgig function (common for all cells)
   NumericVector p = as - 1 / thetaBatch; 
   double b = 2 * bs;
   
   // GIG draws
   /* DEBUG: return original value of s0 if input values are not within the appropriate range (to avoid problems related to numerical innacuracy) */

   // Calculating parameter to the passed as input to the Rgig function (specific to each cell)
   NumericVector a = 2 * nu / thetaBatch;      
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
    arma::vec log_aux = (y - lognu) % (sum_bygene_all + 1/thetaBatch);
    log_aux -= nu0 % (exp(y-lognu)-1)  % (sum(mu(span(q_bio,q -1))) + (1/(thetaBatch % s)));    
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
  int PrintProgress) 
{

  // NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int q = Counts.nrow(); int qbio = delta0.size(); int Naux = N/thin - burn/thin;
  int nBatch = BatchDesign.ncol();
  
  // TRANSFORMATION TO ARMA ELEMENTS 
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts); arma::mat BatchDesign_arma = as_arma(BatchDesign);
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
  arma::vec phiAux = as_arma(phi0); double LSphiAux = LSphi0; Rcpp::List phiAuxList;
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
    
    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, muAux.col(0), deltaAux.col(0),
                           nuAux.col(0), p_Phi, sumByGeneBio_arma, qbio,n); 
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]); if(i>=burn) {phiAccept += Rcpp::as<double>(phiAuxList["ind"]);}

    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdateBatch(thetaAux.col(0), exp(LSthetaAux), BatchDesign_arma,
                                sAux, nuAux.col(0), atheta, btheta, n, nBatch);
    PthetaAux += thetaAux.col(1); if(i>=burn) {thetaAccept += thetaAux.col(1);}
    
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR    
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, deltaAux.col(0), 
                     phiAux, nuAux.col(0), sumByCellBio_arma, s2mu, q, qbio);     
    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}
    
    // UPDATE OF S
    sAux = sUpdateBatch(sAux, nuAux.col(0), thetaAux.col(0), as, bs, BatchDesign_arma, n);  
    
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, 
                           muAux.col(0), nuAux.col(0), phiAux, adelta, bdelta, qbio, n);  
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);} 
    
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdateBatch(nuAux.col(0), exp(LSnuAux), Counts_arma, BatchDesign_arma,
                     muAux.col(0), deltaAux.col(0),
                     phiAux, sAux, thetaAux.col(0), sumByGeneAll_arma, q, qbio, n); 
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/50; PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux=LSmuAux+PmuAux*std::min(0.01,1/sqrt(i)); 
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.01,1/sqrt(i));                 
        PphiAux=PphiAux/50; PphiAux = -1+2*(PphiAux>ar); 
        LSphiAux=LSphiAux - PphiAux*std::min(0.01,1/sqrt(i));  
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux=LSnuAux+PnuAux*std::min(0.01,1/sqrt(i)); 
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*arma::conv_to<arma::mat>::from(PthetaAux>ar); 
        LSthetaAux= LSthetaAux + PthetaAux*std::min(0.01,1/sqrt(i));
                
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