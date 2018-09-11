#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>

using namespace Eigen;
using namespace std;

double QForce(const VectorXd &Xa, const VectorXd &e, const MatrixXd &W, double sigma_sqrd, int i) {
    return 2*(W.row(i)*e - Xa(i))/sigma_sqrd;
}

double GreenFuncSum(const VectorXd &X, const VectorXd &X_new, const VectorXd &X_newa, const VectorXd &Xa, const VectorXd &e,\
                    const MatrixXd &W, double sigma_sqrd, double timestep, int D, double Diff) {
    double GreenSum  = 0;
    int P = X.size()/D;

    for(int i=0; i<P; i++) {
        double GreenFunc = 0;
        for(int j=0; j<D; j++) {
            int index = D*i+j;
            double QForceOld = QForce(Xa, e, W, sigma_sqrd, index);
            double QForceNew = QForce(X_newa, e, W, sigma_sqrd, index);
            GreenFunc += 0.5*(QForceOld + QForceNew) * (0.5*Diff*timestep*(QForceOld - QForceNew)-X_new(index)+X(index));
        }
        GreenSum += exp(GreenFunc);
    }
    return GreenSum;
}
