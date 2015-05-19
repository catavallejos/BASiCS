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
  for (int i=0; i<arma::size(x,0); i++)
  {
    for(int j=0; j<arma::size(x,1); j++)
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
  arma::vec const& sum_bycell_bio) /* Sum of expression counts by cell (biological genes only) */
{
    using arma::span;

    // CREATING VARIABLES WHERE TO STORE DRAWS
    int q = arma::size(mu0,0);
    int q_bio = arma::size(delta,0);
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
  arma::mat const& prop_var,  /* Current value of the proposal variances for $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */  
  double const& a_delta, /* Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  double const& b_delta) /* Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
{
    using arma::span;
    
    // CREATING VARIABLES WHERE TO STORE DRAWS
    int n = size(nu,0);
    int q_bio = arma::size(delta0,0);
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

/* Auxiliary function required for some of the Metropolis-Hastings updates of kappa */
arma::mat UpdateAux_kappaTrick(
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)'$ */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& nu) /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
{
  arma::mat x = ((phi % nu) * mu.t()).t();
  x.each_col() += 1/delta;
  return x;
}

/* Metropolis-Hastings updates of kappa 
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::mat kappaUpdate(
  arma::vec const& kappa0, /* Current value of $\kappa=(\kappa_1,...,\kappa_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\kappa=(\kappa_1,...,\kappa_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& s2_kappa, /* Variance hyper-parameter of the Normal($0$,$s^2_{\kappa}$) prior assigned to each $\kappa_j$ */
  arma::vec const& sum_bygene_bio) /* Sum of expression counts by gene (biological genes only) */
{
    using arma::span;
    using arma::pow;
  
    // CREATING VARIABLES WHERE TO STORE DRAWS
    int n = size(kappa0,0);
    int q_bio = arma::size(delta,0);
    arma::vec logphi0 = log(n) + kappa0 - log_sum_exp_cpp(kappa0);
    arma::vec phi0 = exp(logphi0);
    
    // PROPOSAL STEP
    arma::vec y = arma::randn(n) % sqrt(prop_var) + kappa0;
    y(0) = 0;     
    arma::vec logy = log(n) + y - log_sum_exp_cpp(y);
    arma::vec u = arma::randu(n);
       
    // ACCEPT/REJECT STEP   
    arma::vec log_aux = (logy-logphi0) % sum_bygene_bio - 0.5*(pow(y,2)-pow(kappa0,2))/s2_kappa;
    arma::mat m = Counts.rows(0, q_bio - 1);
    m.each_col() += 1 / delta;
    arma::mat num = UpdateAux_kappaTrick(exp(logy), mu(span(0, q_bio - 1)), delta, nu);
    num /= UpdateAux_kappaTrick(phi0, mu(span(0, q_bio - 1)), delta, nu);    
    m %= log(num);
    log_aux -= sum(m, 0).t();
    arma::umat ind = log(u) < log_aux;
    // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
    // DEBUG: Print warning message
    ind.elem(find_nonfinite(log_aux)).zeros();
    if(size(find_nonfinite(log_aux),0)>0)
    {
      Rcpp::Rcout << "Something went wrong when updating kappa" << size(find_nonfinite(log_aux),0) << std::endl;
      Rcpp::stop("Please consider additional filter of the input dataset.");
    }

    // CREATING OUTPUT VARIABLE
    arma::vec kappa = ind % y + (1 - ind) % kappa0;
    kappa(0) = 0; /* To mantain $\kappa_1=0$ */
   
    // TRANFORMING KAPPA IN TERMS OF PHI
    arma::vec phi=exp(kappa)/exp(log_sum_exp_cpp(kappa));
    
    // OUTPUT
    // DEBUG: If the sum of phi is not 1 (subject to <0.001 error margin), reject proposal (to avoid numerical innacuracies)
    if(sum(phi) > 0.999 && sum(phi) < 1.001) {return join_rows(join_rows(kappa, arma::conv_to<arma::mat>::from(ind)),n*phi); }
    else {return join_rows(join_rows(kappa0, arma::zeros(n)),n*phi0); }
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
  double const& bs) /* Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
{
   // Allowing the use of rgig within C++ code (transform arma::vec elements into NumericVector elements)
   NumericVector s0 = wrap(s0_arma);
   NumericVector nu = wrap(nu_arma);
   
   // Creating a vector where draws will be stored
   int n = s0.size();
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
  arma::vec const& sum_bygene_all) /* Sum of expression counts by gene (all genes) */
{
    using arma::span;
    
    // CREATING VARIABLES WHERE TO STORE DRAWS
    int q = size(mu,0);
    int q_bio = arma::size(delta,0);
    int n = size(nu0,0);
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
    if(size(find_nonfinite(log_aux),0)>0) {Rcpp::Rcout << size(find_nonfinite(log_aux),0) << std::endl;}

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
  double const& b_theta) /* Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  {
    using arma::span;

    // CREATING VARIABLES WHERE TO STORE DRAWS
    int n = size(nu,0);
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
  NumericVector kappa0, // Starting value of $\kappa=(\kappa_1,...,\kappa_n)$'$ 
  NumericVector s0, // Starting value of $s=(s_1,...,s_n)$'$ 
  NumericVector nu0, // Starting value of $\nu=(\nu_1,...,\nu_n)$'$   
  double theta0, // Starting value of $\theta$ 
  double adelta, // Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  double bdelta, // Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ 
  double s2kappa, // Variance hyper-parameter of the Normal($0$,$\sigma^2_{\kappa}$) prior assigned to each $\kappa_j$ 
  double as, // Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double bs, // Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$  
  double atheta, // Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
  double btheta, // Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
  double ar, // Optimal acceptance rate for adaptive Metropolis-Hastings updates
  NumericVector LSmu0, // Starting value of adaptive proposal variance of $\mu=(\mu_1,...,\mu_q)'$ (log-scale)
  NumericVector LSdelta0, // Starting value of adaptive proposal variance of $\delta=(\delta_1,...,\delta_{q_0})'$ (log-scale)
  NumericVector LSkappa0, // Starting value of adaptive proposal variance of $\kappa=(\kappa_1,...,\kappa_n)'$ (log-scale)
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
  arma::mat nu = arma::zeros(Naux,n); 
  arma::vec theta = arma::zeros(Naux); 
  arma::mat delta = arma::zeros(Naux,qbio); 
  arma::mat kappa = arma::zeros(Naux,n); 
  arma::mat s = arma::zeros(Naux,n);  
  arma::mat LSmu;
  arma::mat LSnu;
  arma::vec LStheta; 
  arma::mat LSdelta;
  arma::mat LSkappa;
  // LOG-PROPOSAL VARIANCES 
  if(StoreAdapt == 1)
  {
    LSmu = arma::zeros(Naux,q); 
    LSnu = arma::zeros(Naux,n); 
    LStheta = arma::zeros(Naux); 
    LSdelta = arma::zeros(Naux,qbio); 
    LSkappa = arma::zeros(Naux,n);     
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(q); arma::vec PmuAux = arma::zeros(q);
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  double thetaAccept=0; double PthetaAux=0;
  arma::vec deltaAccept = arma::zeros(qbio); arma::vec PdeltaAux = arma::zeros(qbio);
  arma::vec kappaAccept = arma::zeros(n); arma::vec PkappaAux = arma::zeros(n);
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(q,2); muAux.col(0)=as_arma(mu0); arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::vec thetaAux = arma::zeros(2); thetaAux(0) = theta0; double LSthetaAux = LStheta0;
  arma::mat deltaAux = arma::zeros(qbio,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0); 
  arma::mat kappaAux = arma::zeros(n,3); kappaAux.col(0)=as_arma(kappa0); kappaAux.col(2)=n * exp(kappaAux.col(0))/exp(log_sum_exp_cpp(kappaAux.col(0))); arma::vec LSkappaAux = as_arma(LSkappa0);
  arma::vec sAux = as_arma(s0); arma::vec phiAux;
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q); arma::vec PnuAux0 = arma::zeros(n); double PthetaAux0=0; arma::vec PdeltaAux0 = arma::zeros(qbio); arma::vec PkappaAux0 = arma::zeros(n); 
  
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
    
    // UPDATE OF KAPPA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR, 3rd COLUMN IS PHI (AS A FUNCTION OF KAPPA)
    kappaAux = kappaUpdate(kappaAux.col(0),exp(LSkappaAux),Counts_arma,muAux.col(0),deltaAux.col(0),nuAux.col(0),s2kappa,sumByGeneBio_arma); 
    PkappaAux += kappaAux.col(1); if(i>=burn) {kappaAccept += kappaAux.col(1);}
    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
    thetaAux = thetaUpdate(thetaAux(0),exp(LSthetaAux),sAux,nuAux.col(0),atheta,btheta); 
    PthetaAux += thetaAux(1); if(i>=burn) {thetaAccept += thetaAux(1);}
    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR    
    muAux = muUpdate(muAux.col(0),exp(LSmuAux),Counts_arma,deltaAux.col(0),kappaAux.col(2),nuAux.col(0),sumByCellBio_arma);     
    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}
    // UPDATE OF S
    sAux = sUpdate(sAux,nuAux.col(0),thetaAux(0),as,bs);    
    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    deltaAux = deltaUpdate(deltaAux.col(0),exp(LSdeltaAux),Counts_arma,muAux.col(0),nuAux.col(0),kappaAux.col(2),adelta,bdelta);  
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);}    
    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdate(nuAux.col(0),exp(LSnuAux),Counts_arma,muAux.col(0),deltaAux.col(0),kappaAux.col(2),sAux,thetaAux(0),sumByGeneAll_arma); 
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/50; PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
        LSmuAux=LSmuAux+PmuAux*std::min(0.01,1/sqrt(i)); 
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
        LSnuAux=LSnuAux+PnuAux*std::min(0.01,1/sqrt(i)); 
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*(PthetaAux>ar); 
        LSthetaAux=LSthetaAux+PthetaAux*std::min(0.01,1/sqrt(i)); 
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
        LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.01,1/sqrt(i)); 
        PkappaAux=PkappaAux/50; PkappaAux = -1+2*arma::conv_to<arma::mat>::from(PkappaAux>ar);
        LSkappaAux=LSkappaAux+PkappaAux*std::min(0.01,1/sqrt(i)); 
        
        Ibatch=0; PmuAux=PmuAux0; 
        PnuAux=PnuAux0; PthetaAux=PthetaAux0; PdeltaAux=PdeltaAux0; 
        PkappaAux=PkappaAux0; 
      }
      
    }
        
    // STORAGE OF DRAWS
    if(i%thin==0 & i>=burn)
    {      
      mu.row(i/thin - burn/thin) = muAux.col(0).t();       
      nu.row(i/thin - burn/thin) = nuAux.col(0).t();       
      theta(i/thin - burn/thin) = thetaAux(0);       
      delta.row(i/thin - burn/thin) = deltaAux.col(0).t();      
      kappa.row(i/thin - burn/thin) = kappaAux.col(0).t(); 
      s.row(i/thin - burn/thin) = sAux.t();    
      
      if(StoreAdapt == 1)
      {
        LSmu.row(i/thin - burn/thin) = LSmuAux.t();
        LSnu.row(i/thin - burn/thin) = LSnuAux.t();
        LStheta(i/thin - burn/thin) = LSthetaAux; 
        LSdelta.row(i/thin - burn/thin) = LSdeltaAux.t();  
        LSkappa.row(i/thin - burn/thin) = LSkappaAux.t();
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
        Rcpp::Rcout << "phi (cell 1): " << kappaAux(0,2) << std::endl;
        Rcpp::Rcout << "s (cell 1): " << sAux(0) << std::endl;
        Rcpp::Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
        Rcpp::Rcout << "theta: " << thetaAux(0) << std::endl;
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
        Rcpp::Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
        Rcpp::Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl; 
        Rcpp::Rcout << "LSphi (cell 2): " << LSkappaAux(1) << std::endl;
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
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << max(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among kappa[j]'s: " << min(kappaAccept(arma::span(1, n - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among kappa[j]'s: " << mean(kappaAccept(arma::span(1, n - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among kappa[j]'s: " << max(kappaAccept(arma::span(1, n - 1))/(N-burn)) << std::endl;
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
           Rcpp::Named("nu")=nu,
           Rcpp::Named("theta")=theta,
           Rcpp::Named("delta")=delta,
           Rcpp::Named("kappa")=kappa,
           Rcpp::Named("s")=s,
           Rcpp::Named("ls.mu")=LSmu,
           Rcpp::Named("ls.nu")=LSnu,
           Rcpp::Named("ls.theta")=LStheta,
           Rcpp::Named("ls.delta")=LSdelta,
           Rcpp::Named("ls.kappa")=LSkappa)); 
  }
  
  else
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu")=mu,
           Rcpp::Named("nu")=nu,
           Rcpp::Named("theta")=theta,
           Rcpp::Named("delta")=delta,
           Rcpp::Named("kappa")=kappa,
           Rcpp::Named("s")=s)); 
  }

}
