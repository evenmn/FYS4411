#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>

using namespace Eigen;
using namespace std;

double QForce(const VectorXd &X, const VectorXd &a, const VectorXd &b, const MatrixXd &W, int N, double sigma, int i) {

    double sigma_sqrd = sigma * sigma;

    VectorXd v = b + (X.transpose() * W).transpose()/(sigma_sqrd);
    VectorXd Xa = X - a;


    double QF = 0;

    QF -= Xa(i);
    for(int j=0; j<N; j++) {
        QF += W(i,j)/(1 + exp(-v(j)));
    }


    //cout << QF*(2/sigma_sqrd) << endl;
    return QF*(2/sigma_sqrd);

}

double GreenFuncSum(const VectorXd &X, const VectorXd &X_new, const VectorXd &a, const VectorXd &b, const MatrixXd &W, int N, double sigma, double timestep, int D, double Diff) {
    double GreenSum = 0;

    int P = X.size()/D;

    /*
    for(int i=0; i<P; i++) {
        double GreenOld = 0;
        double GreenNew = 0;
        for(int j=0; j<D; j++) {
            GreenOld += (X_new(D*i+j) - X(D*i+j) - Diff*timestep*QForce(X, a, b, W, N, sigma)) * \
                        (X_new(D*i+j) - X(D*i+j) - Diff*timestep*QForce(X, a, b, W, N, sigma));
            GreenNew += (X(D*i+j) - X_new(D*i+j) - Diff*timestep*QForce(X_new, a, b, W, N, sigma)) * \
                        (X(D*i+j) - X_new(D*i+j) - Diff*timestep*QForce(X_new, a, b, W, N, sigma));
            //cout << "GO " << GreenOld << endl;
            //cout << "GN " << GreenNew << endl;
        }
        GreenSum += exp(GreenOld/GreenNew);
        //cout << GreenSum << endl;
    }
    */
    for(int i=0; i<P; i++) {
        double GreenFunc = 0;
        for(int j=0; j<D; j++) {
            GreenFunc += 0.5*(QForce(X, a, b, W, N, sigma, D*i+j) + QForce(X_new, a, b, W, N, sigma, D*i+j)) * \
                    (Diff*timestep*0.5*(QForce(X, a, b, W, N, sigma, D*i+j) -\
                    QForce(X_new, a, b, W, N, sigma, D*i+j))-X_new(D*i+j)+X(D*i+j));
        }
        GreenSum += exp(GreenFunc);
    }

    return GreenSum;
}
