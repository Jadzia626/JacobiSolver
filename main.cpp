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

    Mat<double>    mH;
    Col<double>    mV;
    Mat<cx_double> mEigVec;
    Col<cx_double> mEigVal;

    int          iN, iIt, iRot;
    int          iMethod = 0;
    Jacobi       *oJacobi = new Jacobi();
    time_t       tStart, tStop;
    stringstream ssFile;

    // Default Matrix

    if(argc == 4) {
        iMethod = atoi(argv[1]);
        ssFile << "Data/";
        if(atoi(argv[3]) == 0) ssFile << "Std/";
        if(atoi(argv[3]) == 1) ssFile << "Eff/";
        if(atoi(argv[3]) == 2) ssFile << "EffECut/";
        ssFile << "Ham-2P-" << argv[2] << "R.arma";
        mH.load(ssFile.str());
    } else {
        cout << endl;
        cout << "Jacobi takes three arguments: " << endl;
        cout << " #1 : Method" << endl;
        cout << "      0 = Armadillo eig_gen()" << endl;
        cout << "      1 = Simple w/Rotate" << endl;
        cout << "      2 = Simple" << endl;
        cout << "      3 = Cyclic" << endl;
        cout << "      4 = Parallel" << endl;
        cout << " #2 : Number of shells for the Hamiltonian (2-20)" << endl;
        cout << " #3 : Interaction type" << endl;
        cout << "      0 = Standard" << endl;
        cout << "      1 = Effective" << endl;
        cout << "      2 = Effective w/energy cut" << endl;
        cout << endl;
        return 0;
    }

    iN   = (int)mH.n_cols;
    iIt  = 0;
    iRot = 0;

    mV   = zeros(iN);

    cout << endl;
    #ifndef ONLY_LOWEST
        cout << "Initial Matrix:" << endl;
        cout << mH << endl;
    #endif
    cout << "Dimension: " << iN << "x" << iN << endl;
    cout << "Sparsity:  " << oJacobi->Sparsity(&mH) << endl;
    cout << endl;

    switch(iMethod) {

        case 0: // Aramdillo's Eigenvalue Solver

            tStart = clock();
            eig_gen(mEigVal, mEigVec, mH);
            tStop  = clock();

            mV = sort(real(mEigVal));

            cout << "Armadillo's eig_gen():" << endl;
            cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
            cout << "   Iterations:   N/A" << endl;
            cout << "   Rotations:    N/A" << endl;

            break;

        case 1: // Simple Jacobi w/Rotate

            tStart = clock();
            oJacobi->SimpleRotate(&mH, &mV, &iIt, &iRot);
            tStop  = clock();

            mV = sort(mV);

            cout << "Simple Jacobi w/Rotations:" << endl;
            cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
            cout << "   Iterations:   " << iIt << endl;
            cout << "   Rotations:    " << iRot << endl;

            break;

        case 2: // Simple Jacobi

        tStart = clock();
        oJacobi->Simple(&mH, &mV, &iIt, &iRot);
        tStop  = clock();

        mV = sort(mV);

        cout << "Simple Jacobi:" << endl;
        cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
        cout << "   Iterations:   " << iIt << endl;
        cout << "   Rotations:    " << iRot << endl;

        break;

    case 3: // Cyclic Jacobi

        tStart = clock();
        oJacobi->Cyclic(&mH, &mV, &iIt, &iRot);
        tStop  = clock();

        mV = sort(mV);

        cout << "Cyclic Jacobi:" << endl;
        cout << "   Time:         " << ((double)tStop-tStart)/CLOCKS_PER_SEC << " sec" << endl;
        cout << "   Iterations:   " << iIt << endl;
        cout << "   Rotations:    " << iRot << endl;

        break;
    }

    #ifdef ONLY_LOWEST
        if(mV.n_elem > 0) cout << "   Eigenvalue 1: " << mV(0) << endl;
        if(mV.n_elem > 1) cout << "   Eigenvalue 2: " << mV(1) << endl;
        if(mV.n_elem > 2) cout << "   Eigenvalue 3: " << mV(2) << endl;
        if(mV.n_elem > 3) cout << "   Eigenvalue 4: " << mV(3) << endl;
        if(mV.n_elem > 4) cout << "   Eigenvalue 5: " << mV(4) << endl;
    #else
        cout << "   Eigenvalues:" << endl << mV << endl;
    #endif
    cout << endl;

    return 0;
}
