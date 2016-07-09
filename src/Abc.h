#ifndef _ABC_H_
#define _ABC_H_

#include <Rcpp.h>

class Abc {
    private:
        // MH parameters
        double _delta;

        // ABC parameters
        double _epsilon;

        // value 
        Rcpp::NumericVector _theta(1);

    public:
        // propose new theta
        void propose();

        // get theta
        Rcpp::NumericVector theta();

};

#endif
