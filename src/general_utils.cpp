#include "utils.h"

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
  //register 
  double t1,cb,t2;
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
  if ((chi < ZTOL) & (lambda > 0.0)) {
    int i;
    for (i = 0; i < n; i++) {
      samps(i) = R::rgamma(lambda, 2.0/psi);
    }
    return samps;
  }
  
  /* special cases which is basically an inverse gamma distribution */
  if ((psi < ZTOL) & (lambda < 0.0)) { 
    int i;
    for (i = 0; i < n; i++) {
      samps(i) = 1.0/R::rgamma(0.0-lambda, 2.0/chi);
    }
    return samps;
  }
  
  /* special case which is basically an inverse gaussian distribution */
  if (lambda == -0.5) {
    double alpha;
    int i;
    alpha = sqrt(chi/psi);
    for (i = 0; i < n; i++) {
      samps(i) = rinvgauss(alpha, chi);
    }
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
// double RgigDouble(const double lambda, const double chi, const double psi)
// {
//   double samps;
//   /* special case which is basically a gamma distribution */
//   if((chi < ZTOL) & (lambda > 0.0)) {
//     samps = R::rgamma(lambda, 2.0/psi);
//     return samps;
//   }
  
//   /* special cases which is basically an inverse gamma distribution */
//   if((psi < ZTOL) & (lambda < 0.0)) { 
//     samps = 1.0/R::rgamma(0.0-lambda, 2.0/chi);
//     return samps;
//   }
  
//   /* special case which is basically an inverse gaussian distribution */
//   if(lambda == -0.5) {
//     double alpha;
//     alpha = sqrt(chi/psi);
//     samps = rinvgauss(alpha, chi);
//     return samps;
//   }  
  
//   /*
//   * begin general purpose rgig code, which was basically 
//   * translated from the R function rgig in the ghyp package v_1.5.2
//   */
  
//   double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2, Y;
//   int need;
  
//   alpha = sqrt(chi/psi);
//   beta2 = psi*chi;
//   beta = sqrt(psi*chi);
//   lm1 = lambda - 1.0;
//   lm12 = lm1*lm1;
//   m = (lm1 + sqrt(lm12 + beta2))/beta;
//   m1 = m + 1.0/m;
  
//   upper = m;
//   while (gig_y_gfn(upper, m, beta, lambda) <= 0) { upper *= 2.0; }
  
//   yM = zeroin_gig(0.0, m, gig_y_gfn, ZTOL, m, beta, lambda);
//   yP = zeroin_gig(m, upper, gig_y_gfn, ZTOL, m, beta, lambda);
  
//   a = (yP - m) * pow(yP/m, 0.5 * lm1);
//   a *=  exp(-0.25 * beta * (yP + 1.0/yP - m1));
//   b = (yM - m) * pow(yM/m, 0.5 * lm1);
//   b *= exp(-0.25 * beta * (yM + 1/yM - m1));
//   c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);
  
//   need = 1;
//   while (need) 
//   {
//     R1 = unif_rand();
//     R2 = unif_rand();
    
//     Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
//     if (Y > 0.0) 
//     {
//       if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) 
//       {
//         need = 0;
//       }
//     }
//   }
//   samps = Y*alpha;
//   return(samps);
// }



/* END OF ADAPTED CODE*/

void StartSampler(int const& N)
{
  Rcout << "-----------------------------------------------------" << std::endl;  
  Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;  
}

void EndBurn()
{
  Rcout << "-----------------------------------------------------" << std::endl; 
  Rcout << "End of Burn-in period."<< std::endl;
  Rcout << "-----------------------------------------------------" << std::endl; 
}

void CurrentIter(int const& i, int const& N)
{
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << "Iteration " << i << " out of " << N << " has been completed." << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << "Current draws for selected parameters are displayed below." << std::endl;
}  

void EndSampler(int const& N)
{
  Rcout << " " << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << "All " << N << " MCMC iterations have been completed." << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << " " << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
  Rcout << "Please see below a summary of the overall acceptance rates." << std::endl;
  Rcout << "-----------------------------------------------------" << std::endl;
}

void ReportAR(arma::vec const& AR, 
              std::string const& Param)
{
  Rcout << " " << std::endl;  
  Rcout << "Minimum acceptance rate among " << 
    Param << ": " << min(AR) << std::endl;
  Rcout << "Average acceptance rate among " << 
    Param << ": " << mean(AR) << std::endl;
  Rcout << "Maximum acceptance rate among " << 
    Param << ": " << max(AR) << std::endl;
  Rcout << " " << std::endl;  
}

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
  for (unsigned int i=0; i<arma::size(x,0); i++) {
    for(unsigned int j=0; j<arma::size(x,1); j++) {
      output(i,j) = std::lgamma(x(i,j));
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
  for (unsigned int i=0; i < output.size(); i++) {
    output(i) = std::lgamma(x(i));
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

arma::vec DegubInd(arma::vec ind,
                   int const q,
                   arma::vec const& u, 
                   arma::vec const& log_aux,
                   arma::vec const& y,
                   double const& threshold,
                   std::string const& param)
{
  for (int i=0; i < q; i++) {
    if( arma::is_finite(log_aux(i)) & arma::is_finite(y(i)) ) {
      if((log(u(i)) < log_aux(i)) & (y(i) > threshold)) { ind(i) = 1; }
      else{ind(i) = 0;}            
    }
    else {
      ind(i) = 0; 
      Rcpp::Rcout << "Error when updating " << param << " " << i+1 << std::endl;
      Rcpp::Rcout << "Consider applying additional quality control" << std::endl; 
      Rcpp::Rcout << "to remove genes/cells when low total counts." << std::endl; 
    }
  }  
  return ind;
}

/*Dirichlet sampler*/
arma::vec rDirichlet(arma::vec alpha) {
  arma::vec aux = arma::ones(alpha.size());
  unsigned int i;
  for(i=0; i<alpha.size(); i++) {
    // shape = alpha(i) ; scale = 1
    // See https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
    aux(i) = R::rgamma(alpha(i),1);   
  }
  return aux / sum(aux);  
}

// Multivariate normal distribution
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
