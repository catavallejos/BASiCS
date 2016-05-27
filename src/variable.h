#ifndef VARIABLE_H
#define VARIABLE_H

#include <RcppArmadillo.h>

using arma::vec;

class Variable
{
public:
    vec u;
    vec y;
    vec logproposal;
    vec proposal_variances;
    vec current_value;
    vec previous_value;
    arma::Col<int> valid;

    Variable(int n) : u(n),
                      y(n),
                      logproposal(n),
                      proposal_variances(n),
                      current_value(n),
                      previous_value(n)
    {
        // pass
    }
};

#endif
