// #include "delta.h"
//
// #include <cmath>
//
// using vec;
// using Col;
// using is_finite;
//
// void clear_variable(vec const& old_var, vec const& new_var,
//                     Col<int>& indicator, vec const& unif,
//                     vec const& log_aux)
// {
//     for (int i=0; i < new_var.size(); i++)
//     {
//         if(is_finite(log_aux(i)))
//         {
//             if(log(unif(i)) < log_aux(i))
//             {
//                 // DEBUG: Reject very small values to avoid numerical issues
//                 if (is_finite(y(i)) && y(i) > 1e-3)
//                 {
//                     indicator(i) = 1;
//                     new_var(i) = y(i);
//                 }
//                 else
//                 {
//                     indicator(i) = 0;
//                     new_var(i) = old_var(i);
//                 }
//             } else
//             {
//                 indicator(i) = 0;
//                 new_var(i) = old_var(i);
//             }
//         } else
//         {
//             // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
//             // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
//             indicator(i) = 0;
//             new_var(i) = old_var(i);
//             // Rcpp::Rcout << "Something went wrong when updating mu " << i+1 << std::endl;
//             // Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");
//         }
//     }
// }
//
// void DeltaVar::update(int q_bio, vec const& sum_bycell_bio)
// {
//     arma::span range = = arma::span(0, q_bio - 1);
//     vec const& prop_var = this->proposal_variances;
//
//     // PROPOSAL STEP
//     vec y = arma::randn(q_bio) % sqrt(prop_var(range));
//     y += log(mu0(range));
//     y = exp(y);
//     vec u = arma::randu(q_bio);
//
//     // ACCEPT/REJECT STEP
//     vec log_aux = (log(y) - log(mu0(range))) % sum_bycell_bio;
//
//     // Loop to replace matrix operations, through genes and cells
//     for (int i=0; i < q_bio; i++)
//     {
//       for (int j=0; j < n; j++)
//       {
//         log_aux(i) -= ( Counts(i,j) + (1/delta(i)) ) *
//                       log( ( phi(j)*nu(j)*y(i)+(1/delta(i)) ) / ( phi(j)*nu(j)*mu0(i)+(1/delta(i)) ));
//       }
//     }
//
//     log_aux -= (0.5/s2_mu) * (pow(log(y),2) - pow(log(mu0(range)),2));
//
//     // CREATING OUTPUT VARIABLE & DEBUG
//     for (int i=0; i < q_bio; i++)
//     {
//       if(is_finite(log_aux(i)))
//       {
//         if(log(u(i)) < log_aux(i))
//         {
//           // DEBUG: Reject very small values to avoid numerical issues
//           if(is_finite(y(i)) & (y(i) > 1e-3)) {ind(i) = 1; mu(i) = y(i);}
//           else{ind(i) = 0; mu(i) = mu0(i);}
//         }
//         else{ind(i) = 0; mu(i) = mu0(i);}
//       }
//       // DEBUG: Reject values such that acceptance rate cannot be computed (due no numerical innacuracies)
//       // DEBUG: Reject values such that the proposed value is not finite (due no numerical innacuracies)
//       else
//       {
//         ind(i) = 0; mu(i) = mu0(i);
//         Rcpp::Rcout << "Something went wrong when updating mu " << i+1 << std::endl;
//         Rcpp::warning("If this error repeats often, please consider additional filter of the input dataset or use a smaller value for s2mu.");
//       }
//     }
//
//     // OUTPUT
//     return join_rows(mu, ind);
//
// }
