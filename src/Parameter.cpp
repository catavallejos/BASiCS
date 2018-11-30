// Abstract class to define parameter behaviour
// All we need these to do is to have an update method
// And methods to return the values and the acceptances.
// Constructor should initialise matrix for storing values
// Chain has to ask it to store values, since we don't care about Burn otherwise
#include <RcppArmadillo.h>
class Parameter {
  public:
    // pure virtual function
    virtual void update() {
      arma::vec proposal = propose();
      arma::vec acceptance = accept(proposal);
      // This should be the DegubInd function more or less
      // Will need to change acceptance(i) under some conditions?
      for (unsigned int i = 0; i < proposal.size(); i ++) {
        if (acceptance(i) == 1) {
          this->value(i) = proposal(i);
        }
      }
      this->acceptance += acceptance;
    };

    arma::vec getValue() {
      return(value);
    }

    arma::vec getAcceptance() {
      return(acceptance);
    }

    virtual void reportAR() = 0;

  private:
    virtual arma::vec propose() = 0;
    virtual arma::vec accept(arma::vec proposal) = 0;
    arma::mat values;
    arma::vec value;
    int valueIndex;
    arma::vec acceptance;
};
