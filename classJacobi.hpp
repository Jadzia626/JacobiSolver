/*
**  Class : Jacobi
** ~~~~~~~~~~~~~~~~
*/

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <armadillo>

#define ZERO     1e-15
#define CONVERGE 1e-4
#define MAX_IT   1e6

namespace eigensolvers {

class Jacobi {

    public:

    void   SimpleRotate(arma::Mat<double>*, arma::Col<double>*, int*, int*);
    void   Simple(arma::Mat<double>*, arma::Col<double>*, int*, int*);
    void   Cyclic(arma::Mat<double>*, arma::Col<double>*, int*, int*);
    void   Parallel(arma::Mat<double>*, arma::Col<double>*, int*, int*);
    double Sparsity(arma::Mat<double>*);

    private:

    double fLFrobenius(arma::Mat<double>*);
    void   fSearch(arma::Mat<double>*, int*, int*);
    double fLFrobeniusNSearch(arma::Mat<double>*, int*, int*);

};

} // End namespace eigensolvers
