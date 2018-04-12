#include "wavefunction.h"
#include <cmath>
#include <iostream>
#include <string.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

int WaveFunction::setTrialWF(int N, int M, double sigma, double omega)
{
    m_N = N;
    m_M = M;
    m_sigma_sqrd = sigma*sigma;
    m_omega_sqrd = omega*omega;
}

double WaveFunction::Psi_value_sqrd(VectorXd a, VectorXd b, VectorXd X, MatrixXd W)
{
    //Unnormalized wave function

    VectorXd v = b + (W.transpose()*X)/m_sigma_sqrd;
    VectorXd Xa = X - a;

    double prod = 1;
    for(int i=0; i<m_N; i++) {
        prod *= 1 + exp(v(i));
    }

    double exp_ret = exp((Xa.transpose() * Xa));

    return pow(1/(m_sigma_sqrd), exp_ret) * prod * prod;       //Weird pow
}

double WaveFunction::EL_calc(VectorXd X, VectorXd a, VectorXd b, MatrixXd W) {
    // Local energy calculations

    VectorXd v = b + (W.transpose()*X)/m_sigma_sqrd;
    VectorXd Xa = X - a;

    // Kinetic energy
    VectorXd e = VectorXd::Zero(m_N);
    for(int i=0; i<m_N; i++) {
        e(i) = 1/(1 + exp(-v(i)));
    }

    double E = 0;
    for(int i=0; i<m_N; i++) {
        double weird_E = (Xa.transpose() * W.col(i));// * e(i);
        E += weird_E * e(i); //(Xa.transpose() * W.col(i)) * e(i);
        double weird_E2 =((W.col(i)).transpose() * W.col(i));// * e(i) * e(i);
        E += weird_E2 * e(i) * e(i); //((W.col(i)).transpose() * W.col(i)) * e(i) * e(i);
        for(int j=0; j<m_N; j++) {
            double weird_E3 = ((W.col(i)).transpose() * W.col(j));// * e(i) * e(j);
            E += weird_E3 * e(i) * e(j);
        }
    }


    E -= m_N * m_sigma_sqrd;
    E += Xa.transpose() * Xa;
    E = E/(2 * m_sigma_sqrd * m_sigma_sqrd);

    // External potential
    double weird_scalar = (X.transpose() * X);          //Remove weird scalar
    E += weird_scalar * m_omega_sqrd/ 2;

    // Interaction energy
    //MatrixXd norm;
    //rij(X, norm, D);        // Create distance matrix

    //E += sum(sum(norm));

    return E;
}

void WaveFunction::Gradient_a(VectorXd X, VectorXd a, VectorXd &da) {

    VectorXd Xa = X - a;

    da = Xa/m_sigma_sqrd;
}

void WaveFunction::Gradient_b(VectorXd b, VectorXd X, MatrixXd W, VectorXd &db) {

    VectorXd v = b + (W.transpose()*X)/m_sigma_sqrd;

    db.resize(m_N);
    for(int i=0; i<m_N; i++) {
        db(i) = 1/(1 + exp(-v(i)));
    }
}

void WaveFunction::Gradient_W(VectorXd X, VectorXd b, MatrixXd W, MatrixXd &dW) {

    VectorXd v = b + (W.transpose()*X)/m_sigma_sqrd;

    dW.resize(m_M, m_N);
    for(int i=0; i<m_N; i++) {
        for(int j=0; j<m_N; j++) {
            dW(i,j) = X(i)/(1 + exp(-v(j)));
        }
    }
}
