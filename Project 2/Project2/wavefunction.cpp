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


double WaveFunction::EL_calc(VectorXd X, VectorXd Xa, VectorXd v, MatrixXd W, int D, int interaction, double &E_k, double &E_ext, double &E_int) {
    // Local energy calculations

    double E = 0;
    // Kinetic energy
    if(m_sampling==2) {
        VectorXd e = VectorXd::Zero(m_N);
        VectorXd eNominator = VectorXd::Zero(m_N);
        for(int i=0; i<m_N; i++) {
            double expi = exp(-v(i));
            eNominator(i) = expi;
            e(i) = 1/(1 + expi);
        }

        for(int i=0; i<m_N; i++) {
            E -= 0.5*(double) (Xa.transpose() * W.col(i)) * e(i);
            E += 0.5*(double) ((W.col(i)).transpose() * W.col(i)) * eNominator(i) * e(i) * e(i);
            for(int j=0; j<m_N; j++) {
                E += 0.25*(double) ((W.col(i)).transpose() * W.col(j)) * e(i) * e(j);
            }
        }

        E -= 0.5*m_M * m_sigma_sqrd;
        E += 0.25*Xa.transpose() * Xa;
        E = -E/(2 * m_sigma_sqrd * m_sigma_sqrd);
        E_k += E;
    }

    else {
        VectorXd e = VectorXd::Zero(m_N);
        VectorXd eNominator = VectorXd::Zero(m_N);
        for(int i=0; i<m_N; i++) {
            double expi = exp(v(i));
            eNominator(i) = expi;
            e(i) = 1/(1 + 1./expi);
        }

        for(int i=0; i<m_N; i++) {
            E -= 2*(double) (Xa.transpose() * W.col(i)) * e(i);
            E += (double) ((W.col(i)).transpose() * W.col(i)) * 1/eNominator(i) * e(i)*e(i);
            for(int j=0; j<m_N; j++) {
                E += (double) ((W.col(i)).transpose() * W.col(j)) * e(i) * e(j);
            }
        }

        E -= m_M * m_sigma_sqrd;
        E += Xa.transpose() * Xa;
        E = -E/(2 * m_sigma_sqrd * m_sigma_sqrd);
        E_k += E;
    }

    // Interaction energy
    double E_p = 0;
    if(interaction) E_p += rij(X, D);
    E += E_p;
    E_int += E_p;

    // Harmonic oscillator potential
    E_p = (double) (X.transpose() * X) * m_omega_sqrd/ 2;
    E += E_p;
    E_ext += E_p;

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

void WaveFunction::Gradient_b(const VectorXd &v, VectorXd &db) {

    if(m_sampling==2) {
        for(int i=0; i<m_N; i++)
            db(i) = 0.5/(1 + exp(-v(i)));
    }
    else{
        for(int i=0; i<m_N; i++)
            db(i) = 1/(1 + exp(-v(i)));
    }
}

void WaveFunction::Gradient_W(const VectorXd &X, const VectorXd &v, MatrixXd &dW) {

    if(m_sampling==2) {
        for(int i=0; i<m_N; i++) {
            for(int j=0; j<m_M; j++) {
                dW(j,i) = 0.5*X(j)/(m_sigma_sqrd*(1 + exp(-v(i))));
            }
        }
    }
    else{
        for(int i=0; i<m_N; i++) {
            for(int j=0; j<m_M; j++) {
                dW(j,i) = X(j)/(m_sigma_sqrd*(1 + exp(-v(i))));
            }
        }
    }
}
