#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>

using namespace Eigen;
using namespace std;

double QForce(VectorXd X, VectorXd a, VectorXd b, MatrixXd W, int N, double sigma) {

    double sigma_sqrd = sigma * sigma;

    VectorXd v = b + (X.transpose() * W).transpose()/(sigma_sqrd);
    VectorXd Xa = X - a;

    double QF = 0;
    QF += Xa.sum();
    for(int i=0; i<N; i++) {
        QF += (W.col(i)).sum()/(1 + exp(v(i)));
    }
    return QF*(2/sigma_sqrd);
}

double GreenFuncSum(VectorXd X, VectorXd X_new, VectorXd a, VectorXd b, MatrixXd W, int N, double sigma, double timestep, int D, double Diff) {
    double GreenSum = 0;

    int P = X.size()/D;
    for(int i=0; i<P; i++) {
        double GreenOld = 0;
        double GreenNew = 0;
        for(int j=0; j<D; j++) {
            GreenOld += (X_new(D*i+j) - X(D*i+j) - Diff*timestep*QForce(X, a, b, W, N, sigma)) * \
                        (X_new(D*i+j) - X(D*i+j) - Diff*timestep*QForce(X, a, b, W, N, sigma));
            GreenNew += (X(D*i+j) - X_new(D*i+j) - Diff*timestep*QForce(X_new, a, b, W, N, sigma)) * \
                        (X(D*i+j) - X_new(D*i+j) - Diff*timestep*QForce(X_new, a, b, W, N, sigma));
        }
        GreenSum += exp(GreenOld/GreenNew);
    }
    return GreenSum;
}
