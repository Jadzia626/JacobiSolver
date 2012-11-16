/*
**  Class : Jacobi
** ~~~~~~~~~~~~~~~~
*/

#include "classJacobi.hpp"

using namespace std;
using namespace arma;
using namespace eigensolvers;

/*
** Public Functions
*/

void Jacobi::SimpleRotate(Mat<double> *mA, Col<double> *mV, int *iIt, int *iRot) {

    int p, q, k=0;
    int iDim  = (int)mA->n_cols;
    double dTau, dT, dSin, dCos;
    Col<int> mS;
    mS.zeros(2);

    Mat<double> mJ;

    while(fLFrobenius(mA) >= CONVERGE && k <= MAX_IT) {

        fSearch(mA, &p, &q);

        if(abs(mA->at(p,q)) > 0) {
            dTau = (mA->at(q,q) - mA->at(p,p))/(2*mA->at(p,q));
            if(dTau >= 0) {
                dT = 1/(dTau + sqrt(1+dTau*dTau));
            } else {
                dT = -1/(-dTau + sqrt(1+dTau*dTau));
            }
            dCos = 1/sqrt(1+dT*dT);
            dSin = dT*dCos;
        } else {
            break;
            dCos = 1.0;
            dSin = 0.0;
        }

        mJ.zeros(iDim,iDim);
        mJ.diag() = ones(iDim);

        mJ(p,p) = dCos;
        mJ(q,q) = dCos;
        mJ(p,q) = dSin;
        mJ(q,p) = -dSin;

        *mA = trans(mJ)*(*mA)*mJ;
        k++;
    }

    *iIt  = k;
    *iRot = k;
    *mV   = mA->diag();

    return;
}

void Jacobi::Simple(Mat<double> *mA, Col<double> *mV, int *iIt, int *iRot) {

    int p=0, q=0, k=0;
    int iDim  = (int)mA->n_cols;
    double dTau, dT, dSin, dCos;
    double dApp, dAqq, dApq, dAip, dAiq;
    Col<int> mS;
    mS.zeros(2);

    Mat<double> mJ;

    while(fLFrobeniusNSearch(mA,&p,&q) >= CONVERGE && k < MAX_IT) {

        if(abs(mA->at(p,q)) > 0) {
            dTau = (mA->at(q,q) - mA->at(p,p))/(2*mA->at(p,q));
            if(dTau >= 0) {
                dT = 1/(dTau + sqrt(1+dTau*dTau));
            } else {
                dT = -1/(-dTau + sqrt(1+dTau*dTau));
            }
            dCos = 1/sqrt(1+dT*dT);
            dSin = dT*dCos;
        } else {
            break;
            dCos = 1.0;
            dSin = 0.0;
        }

        dApp = mA->at(p,p);
        dApq = mA->at(p,q);
        dAqq = mA->at(q,q);

        for(int i=0; i<iDim; i++) {
            if(i != p && i != q) {
                dAip = mA->at(i,p);
                dAiq = mA->at(i,q);
                mA->at(i,p) = dAip*dCos - dAiq*dSin;
                mA->at(i,q) = dAiq*dCos + dAip*dSin;
                mA->at(p,i) = mA->at(i,p);
                mA->at(q,i) = mA->at(i,q);
            }
        }

        mA->at(p,p) = dApp*dCos*dCos - 2*dApq*dCos*dSin + dAqq*dSin*dSin;
        mA->at(q,q) = dAqq*dCos*dCos + 2*dApq*dCos*dSin + dApp*dSin*dSin;
        mA->at(p,q) = (dApp - dAqq)*dCos*dSin + dApq*(dCos*dCos - dSin*dSin);
        mA->at(q,p) = mA->at(p,q);

        k++;
    }

    *iIt  = k;
    *iRot = k;
    *mV   = mA->diag();

    return;
}

void Jacobi::Cyclic(Mat<double> *mA, Col<double> *mV, int *iIt, int *iRot) {

    int r=0, k=0;
    int iDim  = (int)mA->n_cols;
    double dTau, dT, dSin, dCos;
    double dApp, dAqq, dApq, dAip, dAiq;
    Col<int> mS;
    mS.zeros(2);

    Mat<double> mJ;

    while(fLFrobenius(mA) >= CONVERGE && r < 3*MAX_IT ) {

        for(int p=0; p<iDim-1; p++) {
            for(int q=p+1; q<iDim; q++) {

                if(abs(mA->at(p,q)) > 0) {
                    dTau = (mA->at(q,q) - mA->at(p,p))/(2*mA->at(p,q));
                    if(dTau >= 0) {
                        dT = 1/(dTau + sqrt(1+dTau*dTau));
                    } else {
                        dT = -1/(-dTau + sqrt(1+dTau*dTau));
                    }
                    dCos = 1/sqrt(1+dT*dT);
                    dSin = dT*dCos;
                } else {
                    break;
                    dCos = 1.0;
                    dSin = 0.0;
                }

                dApp = mA->at(p,p);
                dApq = mA->at(p,q);
                dAqq = mA->at(q,q);

                for(int i=0; i<iDim; i++) {
                    if(i != p && i != q) {
                        dAip = mA->at(i,p);
                        dAiq = mA->at(i,q);
                        mA->at(i,p) = dAip*dCos - dAiq*dSin;
                        mA->at(i,q) = dAiq*dCos + dAip*dSin;
                        mA->at(p,i) = mA->at(i,p);
                        mA->at(q,i) = mA->at(i,q);
                    }
                }

                mA->at(p,p) = dApp*dCos*dCos - 2*dApq*dCos*dSin + dAqq*dSin*dSin;
                mA->at(q,q) = dAqq*dCos*dCos + 2*dApq*dCos*dSin + dApp*dSin*dSin;
                mA->at(p,q) = (dApp - dAqq)*dCos*dSin + dApq*(dCos*dCos - dSin*dSin);
                mA->at(q,p) = mA->at(p,q);

                r++;
            }
        }

        k++;
    }

    *iIt  = k;
    *iRot = r;
    *mV   = mA->diag();

    return;
}

void Jacobi::Parallel(Mat<double> *mA, Col<double> *mV, int *iIt, int *iRot) {
}

double Jacobi::Sparsity(Mat<double> *mA) {

    int iDim  = (int)mA->n_cols;
    int iZero = 0;

    for(int i=0; i<iDim; i++) {
        for(int j=0; j<iDim; j++) {
            if(abs(mA->at(i,j)) <= ZERO) iZero++;
        }
    }

    return iZero/(double)(iDim*iDim);
}

/*
** Private Functions
*/

// Calculates the Frobenius norm of the lower triangular matrix
double Jacobi::fLFrobenius(Mat<double> *mA) {

    int    iDim  = (int)mA->n_cols;
    double dF    = 0.0;

    for(int i=1; i<iDim; i++) {
        for(int j=0; j<i; j++) {
            dF += pow(mA->at(i,j),2);
        }
    }

    return sqrt(dF);
}

// Finds the rows with the two largest elements
void Jacobi::fSearch(Mat<double> *mA, int *p, int *q) {

    int    iDim = (int)mA->n_cols;
    double dMax = 0.0;

    for(int i=1; i<iDim; i++) {
        for(int j=0; j<i; j++) {
            if(abs(mA->at(i,j)) > dMax) {
                dMax = abs(mA->at(i,j));
                if(i>j) {
                    *p = j;
                    *q = i;
                } else {
                    *p = i;
                    *q = j;
                }
            }
        }
    }

    return;
}

// Combines Frobenius norm and search
double Jacobi::fLFrobeniusNSearch(Mat<double> *mA, int *p, int *q) {

    int    iDim  = (int)mA->n_cols;
    double dF    = 0.0;
    double dMax  = 0.0;

    for(int i=1; i<iDim; i++) {
        for(int j=0; j<i; j++) {
            dF += pow(mA->at(i,j),2);
            if(abs(mA->at(i,j)) > dMax) {
                dMax = abs(mA->at(i,j));
                if(i>j) {
                    *p = j;
                    *q = i;
                } else {
                    *p = i;
                    *q = j;
                }
            }
        }
    }

    return sqrt(dF);
}






