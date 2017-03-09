#ifndef PROPOSAL_H
#define PROPOSAL_H

#include <RcppArmadillo.h>
#include "util.h"

class Proposal
{
public:
    arma::vec u;
    arma::vec y;
    arma::vec log_aux;
    arma::Col<int> indicator;
    arma::vec proposal_variances;
    arma::vec current_variable;
    arma::vec previous_variable;

    Proposal(int n) : u(n), y(n), log_aux(n) {}
};

#endif
