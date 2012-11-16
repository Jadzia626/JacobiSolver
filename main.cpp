/*
**  Example : Jacobi
** ~~~~~~~~~~~~~~~~~~
*/

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <ctime>
#include <armadillo>
#include "classJacobi.hpp"

#define ONLY_LOWEST

using namespace std;
using namespace arma;
using namespace eigensolvers;

int main(int argc, char* argv[]) {

    Mat<double> mH;
    Mat<double> mT;
    Col<double> mV;
    int         iN;
    int         iIt, iRot;
    Jacobi      *oJacobi = new Jacobi();
    time_t      tStart, tStop;

    // Default Matrix

    //~ mH.load("Data/Binary-Eff/Ham-2P-2R-Dim3x3.arma");
    //~ mH.load("Data/Binary-Std/Ham-2P-3R-Dim8x8.arma");
    //~ mH.load("Data/Binary-Eff/Ham-2P-3R-Dim8x8.arma");
    //~ mH.load("Data/Binary-Eff-ECut/Ham-2P-3R-Dim8x8.arma");
    //~ mH.load("Data/Binary-Eff/Ham-2P-6R-Dim47x47.arma");
    //~ mH.load("Data/Binary-Eff/Ham-2P-8R-Dim104x104.arma");
    //~ mH.load("Data/Binary-Eff/Ham-2P-10R-Dim195x195.arma");
    mH.load("Data/Binary-Eff/Ham-2P-14R-Dim511x511.arma");
    //~ mH.load("Data/Binary-Eff/Ham-2P-20R-Dim1440x1440.arma");

    iN = (int)mH.n_cols;
    mV = zeros(iN);

    cout << endl;
    #ifndef ONLY_LOWEST
        cout << "Initial Matrix:" << endl;
        cout << mH << endl;
    #endif
    cout << "Dimension: " << iN << "x" << iN << endl;
    cout << "Sparsity:  " << oJacobi->Sparsity(&mH) << endl;
    cout << endl;

    // Aramdillo's Eigenvalue Solver

    Mat<cx_double> mEigVec;
    Col<cx_double> mEigVal;

    tStart = clock();
    eig_gen(mEigVal, mEigVec, mH);
    tStop  = clock();
    mV = sort(real(mEigVal));

    cout << "Armadillo's eig_gen():" << endl;
    cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
    cout << "   Iterations:   N/A" << endl;
    cout << "   Rotations:    N/A" << endl;
    #ifdef ONLY_LOWEST
        cout << "   Eigenvalue 1: " << mV(0) << endl;
        cout << "   Eigenvalue 2: " << mV(1) << endl;
        cout << "   Eigenvalue 3: " << mV(2) << endl;
        cout << "   Eigenvalue 4: " << mV(3) << endl;
        cout << "   Eigenvalue 5: " << mV(4) << endl;
    #else
        cout << "   Eigenvalues:" << endl << mV << endl;
    #endif
    cout << endl;

    // Simple Jacobi w/Rotate
/*
    mT   = mH;
    iIt  = 0;
    iRot = 0;
    tStart = clock();
    oJacobi->SimpleRotate(&mT, &mV, &iIt, &iRot);
    tStop  = clock();
    mV = sort(mV);
    cout << "Simple Jacobi w/Rotations:" << endl;
    cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
    cout << "   Iterations:   " << iIt << endl;
    cout << "   Rotations:    " << iRot << endl;
    #ifdef ONLY_LOWEST
        cout << "   Eigenvalue 1: " << mV(0) << endl;
        cout << "   Eigenvalue 2: " << mV(1) << endl;
        cout << "   Eigenvalue 3: " << mV(2) << endl;
        cout << "   Eigenvalue 4: " << mV(3) << endl;
        cout << "   Eigenvalue 5: " << mV(4) << endl;
    #else
        cout << "   Eigenvalues:" << endl << mV << endl;
    #endif
    cout << endl;
*/
    // Simple Jacobi

    mT   = mH;
    iIt  = 0;
    iRot = 0;
    tStart = clock();
    oJacobi->Simple(&mT, &mV, &iIt, &iRot);
    tStop  = clock();
    mV = sort(mV);

    cout << "Simple Jacobi:" << endl;
    cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
    cout << "   Iterations:   " << iIt << endl;
    cout << "   Rotations:    " << iRot << endl;
    #ifdef ONLY_LOWEST
        cout << "   Eigenvalue 1: " << mV(0) << endl;
        cout << "   Eigenvalue 2: " << mV(1) << endl;
        cout << "   Eigenvalue 3: " << mV(2) << endl;
        cout << "   Eigenvalue 4: " << mV(3) << endl;
        cout << "   Eigenvalue 5: " << mV(4) << endl;
    #else
        cout << "   Eigenvalues:" << endl << mV << endl;
    #endif
    cout << endl;

    // Cyclic Jacobi

    mT   = mH;
    iIt  = 0;
    iRot = 0;
    tStart = clock();
    oJacobi->Cyclic(&mT, &mV, &iIt, &iRot);
    tStop  = clock();
    mV = sort(mV);

    cout << "Cyclic Jacobi:" << endl;
    cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
    cout << "   Iterations:   " << iIt << endl;
    cout << "   Rotations:    " << iRot << endl;
    #ifdef ONLY_LOWEST
        cout << "   Eigenvalue 1: " << mV(0) << endl;
        cout << "   Eigenvalue 2: " << mV(1) << endl;
        cout << "   Eigenvalue 3: " << mV(2) << endl;
        cout << "   Eigenvalue 4: " << mV(3) << endl;
        cout << "   Eigenvalue 5: " << mV(4) << endl;
    #else
        cout << "   Eigenvalues:" << endl << mV << endl;
    #endif
    cout << endl;

    return 0;
}
