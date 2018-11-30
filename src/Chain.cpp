#include <list>
#include <RcppArmadillo.h>
// Abstract class to define a chain
// Chain should hold all settings
// Parameter settings (including priors) should be held in Parameter objects
// Probably worthwhile adding a serialise() method (or similar)
// to create output suitable for R
class Chain {
   public:
    // pure virtual function
    virtual void iterate() = 0;

    std::list<arma::mat> getParameters() {
      return(parameters);
    }

    private:
      std::list<arma::mat> parameters;
};
