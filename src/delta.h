#ifndef DELTA_H
#define DELTA_H

#include "variable.h"
// #include "util.h"

using R::lgammafn;
using arma::mat;

class Delta : public Variable
{
    const mat& counts;
    int has_prior;
public:
    double s2;
    double gamma_shape;
    double gamma_rate;

    Delta(int n, const mat& counts) :
        Variable(n),
        counts(counts),
        has_prior(has_prior),
        s2(0),
        gamma_shape(0),
        gamma_rate(0)
    {
        // pass
    }

    void sample_y()
    {
        int n = this->y.size();
        vec& pv = this->proposal_variances;
        this->y = exp(arma::randn(n) % sqrt(pv) + log(this->previous_value));
    }

    void sample_u()
    {
        int n = this->y.size();
        this->u = arma::randu(n);
    }

    void propose(const vec& phi,
                 const vec& nu,
                 const vec& mu)
    {
        int n = this->y.size();
        vec& p = this->logproposal;
        vec& y = this->y;
        vec& delta0 = this->previous_value;
        const mat& C = this->counts;

        p = vec(1/y).for_each(lgammafn) - vec(1/delta0).for_each(lgammafn);
        p *= -n;

        for (int i=0; i < n; i++)
        {
            for (int j=0; j < n; j++)
            {
                p(i) += lgammafn(C(i,j) + (1/y(i))) - lgammafn(C(i,j) + (1/delta0(i)));
                p(i) -= ( C(i,j) + (1/y(i)) ) *  log(phi(j)*nu(j)*mu(i)+(1/y(i)));
                p(i) += ( C(i,j) + (1/delta0(i)) ) *  log( phi(j)*nu(j)*mu(i)+(1/delta0(i)) );
            }
        }

        // +1 should appear because we update log(delta) not delta. However, it cancels out with the prior.
        p -= n * ( (log(y)/y) - (log(delta0)/delta0) );
        if (this->has_prior == 1)
        {
            p += (log(y) - log(delta0)) * this->gamma_shape - this->gamma_rate * (y - delta0);
        } else
        {
            p -= (0.5/this->s2) * (pow(log(y),2) - pow(log(delta0),2));
        }
    }
};

#endif
