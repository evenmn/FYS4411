#include "wavefunction.h"
#include <cmath>
#include <iostream>
#include "eigen3/Eigen/Dense"

using namespace std;
using namespace Eigen;

double rij(VectorXd X, int D) {
    double Ep = 0;              // Sum 1/rij
    int P = X.size()/D;
    for(int i=0; i<P; i++) {
        for(int j=0; j<i; j++) {
            double dist = 0;
            for(int d=0; d<D; d++) {
                dist += (X(D*i+d)-X(D*j+d))*(X(D*i+d)-X(D*j+d));
            }
            Ep += 1/sqrt(dist);
        }
    }
    return Ep;
}

int WaveFunction::setTrialWF(int N, int M, int sampling, double sigma_sqrd, double omega)
{
    m_N          = N;
    m_M          = M;
    m_sampling   = sampling;
    m_sigma = sqrt(sigma_sqrd);
    m_sigma_sqrd = sigma_sqrd;
    m_omega_sqrd = omega*omega;
}

double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v)
{
    //Unnormalized wave function
    double prod = 1;
    for(int i=0; i<m_N; i++) {
        prod *= (1 + exp(v(i)));
    }

    return exp(-(double) (Xa.transpose() * Xa)/(m_sigma_sqrd)) * prod * prod;
}


double WaveFunction::EL_calc(VectorXd X, VectorXd Xa, VectorXd v, MatrixXd W, int D, int interaction, \
                             double &E_kin, double &E_ext, double &E_int) {
    // Local energy calculations

    double E = 0;
    E_kin = 0;
    E_ext = 0;
    E_int = 0;

    double E_k = 0;
    double E_e = 0;
    double E_i = 0;

    VectorXd e_n = VectorXd::Zero(m_N);
    VectorXd e_p = VectorXd::Zero(m_N);
    for(int i=0; i<m_N; i++) {
        e_p(i) = 1/(1 + exp(v(i)));
        e_n(i) = 1/(1 + exp(-v(i)));
    }

    // Kinetic energy
    if(m_sampling==2) {
        for(int i=0; i<m_N; i++) {
            E_k += 0.5*(double) ((W.col(i)).transpose() * W.col(i)) * e_p(i) * e_n(i);
        }

        E_k -= (0.5/m_sigma_sqrd)*(Xa.transpose()*W)*e_n;
        E_k += 0.25*((W.transpose()*W).cwiseProduct(e_n*e_n.transpose())).sum();

        E_k -= 0.5*m_M * m_sigma_sqrd;
        E_k += 0.25*Xa.transpose() * Xa;
        E_k = -E_k/(2 * m_sigma_sqrd * m_sigma_sqrd);
    }

    else {
        for(int i=0; i<m_N; i++) {
            E_k += (double) ((W.col(i)).transpose() * W.col(i)) * e_p(i)*e_n(i);
        }

        E_k -= (2/m_sigma_sqrd)*(Xa.transpose()*W)*e_n;
        E_k += ((W.transpose()*W).cwiseProduct(e_n*e_n.transpose())).sum();

        E_k -= m_M * m_sigma_sqrd;
        E_k += Xa.transpose() * Xa;
        E_k = -E_k/(2 * m_sigma_sqrd * m_sigma_sqrd);

    }

    // Interaction energy
    if(interaction) E_i += rij(X, D);

    // Harmonic oscillator potential
    E_e += (double) (X.transpose() * X) * m_omega_sqrd/ 2;

    E = E_k + E_e + E_i;
    E_kin = E_k;
    E_ext = E_e;
    E_int = E_i;
    return E;
}

void WaveFunction::Gradient_a(const VectorXd &Xa, VectorXd &da) {

    if(m_sampling==2) {
        da = 0.5*Xa/m_sigma_sqrd;
    }
    else{
        da = Xa/m_sigma_sqrd;
    }

}

void WaveFunction::Gradient_b(const VectorXd &e, VectorXd &db) {

    if(m_sampling==2) {
        db = 0.5*e;
    }
    else{
        db = e;
    }
}

void WaveFunction::Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW) {

    if(m_sampling==2) {
        dW = 0.5*X*e.transpose()/m_sigma_sqrd;
    }
    else{
        dW = X*e.transpose()/m_sigma_sqrd;
    }
}


double WaveFunction::Gradient_sigma(const MatrixXd &W, const VectorXd &Xa, const VectorXd &X, const VectorXd &e) {

    double constant = ((W.transpose()*X).transpose()*e);
    double Xa_sqrd  = (double)(Xa.transpose()*Xa);

    if(m_sampling==2) {
        return (0.5*Xa_sqrd + constant)/(m_sigma_sqrd*m_sigma);
    }
    else{
        return (Xa_sqrd + 2*constant)/(m_sigma_sqrd*m_sigma);
    }
}
