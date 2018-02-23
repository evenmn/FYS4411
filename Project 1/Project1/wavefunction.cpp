#include "wavefunction.h"
#include <cmath>
#include <iostream>
#include <string.h>

using namespace std;

double u_der(double dist, double a) {
    return (a/(dist * dist))/(dist - a);
}

double u_secder(double dist, double a) {
    double C = (1/(1-a/dist));
    return -(a/(dist*dist*dist))*C*(2+(a/dist)*C);
}

double V_ext(double pos[3], bool HO, double omega_HO, double omega_z) {
    if(HO) {
        return 0.5*omega_HO*omega_HO*(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    }
    else {
        return 0.5*(omega_HO*omega_HO*(pos[0]*pos[0] + pos[1]*pos[1]) + omega_z*omega_z*pos[2]*pos[2]);
    }
}


int WaveFunction::setTrialWF(int dim, int N, double a, int n_or_a, bool HO)
{

    m_dim    = dim;     //number of dimensions
    m_N      = N;       //number of particles
    m_n_or_a = n_or_a;  //analytical (0) or numerical (1) E_L calculation
    m_a      = a;       //dimension of trap
    m_HO     = HO;      //spherical (true) or elliptical (false) harmonic oscillator
}

double WaveFunction::Psi_value(double pos_mat[][3], double alpha, double beta)
{

    //Returns the wavefunction
    double sumsqrt = 0;         // x1^2 + y1^2 + B*z1^2 + ... + B*zn^2
    double prodf = 1;
    double distx, disty, distz;

    for(int i=0; i<m_N; i++) {
        distx = pos_mat[i][0]*pos_mat[i][0];
        disty = pos_mat[i][1]*pos_mat[i][1];
        distz = pos_mat[i][2]*pos_mat[i][2];
        sumsqrt += distx + disty + beta * distz;

        for(int j=m_N-1; j>i; j--) {
            double norm = sqrt(distx + disty + distz);

            if(norm > m_a){
                prodf *= 1 - m_a/norm;
            }
            else{
                prodf = 0;
                break;
            }
        }
    }
    return exp(-alpha*sumsqrt)*prodf;

}

double WaveFunction::Psi_value_sqrd(double pos_mat[][3], double alpha, double beta)
{
    //Returns the wavefunction squared
    double sumsqrt = 0;         // x1^2 + y1^2 + B*z1^2 + ... + B*zn^2
    double prodf = 1;
    double distx, disty, distz;

    for(int i=0; i<m_N; i++) {
        distx = pos_mat[i][0]*pos_mat[i][0];
        disty = pos_mat[i][1]*pos_mat[i][1];
        distz = pos_mat[i][2]*pos_mat[i][2];
        sumsqrt += distx + disty + beta * distz;

        for(int j=m_N-1; j>i; j--) {
            double norm = sqrt(distx + disty + distz);

            if(norm > m_a){
                prodf *= 1 - m_a/norm;
            }
            else{
                prodf = 0;
                break;
            }
        }
    }
    return exp(-alpha*sumsqrt*2) * prodf * prodf;
}


double WaveFunction::E_L_ana(double pos_mat[][3], double alpha, double beta, double omega_HO, double omega_z)
{
    double distij;
    double distik;
    double r_ij[3];
    double r_ik[3];
    double EL;
    double E_TOT = 0;

    for(int i=0; i<m_N; i++) {
        EL = 0;
        EL += 4*alpha*alpha * (pos_mat[i][0]*pos_mat[i][0] + pos_mat[i][1]*pos_mat[i][1] + \
              beta*beta * pos_mat[i][2]*pos_mat[i][2]);

        if(m_dim==1) {
            EL += -2*alpha;
        }
        else if(m_dim==2) {
            EL += -4*alpha;
        }

        else if(m_dim==3) {
            EL += -4*alpha - 2*alpha*beta;
        }

        for(int j=m_N-1; j>i; j--) {
            for(int l = 0; l < m_dim; l++)
                r_ij[l] = pos_mat[i][l] - pos_mat[j][l];
            distij = sqrt(r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2]);
            double u_der_ij = u_der(distij, m_a);

            EL += -4*alpha*(pos_mat[i][0] * r_ij[0] + pos_mat[i][1] * r_ij[1] +\
                            pos_mat[i][2] * r_ij[2] * beta) * (u_der_ij);

            for(int k=m_N-1; k>i; k--) {
                for(int l = 0; l < m_dim; l++)
                    r_ik[l] = pos_mat[i][l] - pos_mat[k][l];
                distik = sqrt(r_ik[0]*r_ik[0] + r_ik[1]*r_ik[1] + r_ik[2]*r_ik[2]);
                EL += (r_ik[0] * r_ij[0] + r_ik[1] * r_ij[1] + r_ik[2] * r_ij[2]) * \
                      u_der_ij * u_der(distik, m_a);
            }
            EL += u_secder(distij, m_a) + 2 * u_der_ij;
        }
        E_TOT += -0.5*EL + V_ext(pos_mat[i], m_HO, omega_HO, omega_z);
    }
    return E_TOT;
}

double WaveFunction::E_L_num(double pos_mat[][3], double alpha, double beta, double omega_HO, double omega_z, double h)
{

    // Kinetic energy
    double pos_mat_plus[m_N][3];
    double pos_mat_minus[m_N][3];

    memcpy(pos_mat_plus, pos_mat, sizeof(pos_mat_plus));
    memcpy(pos_mat_minus, pos_mat, sizeof(pos_mat_minus));

    double waveFunctionMinus = 0;
    double waveFunctionPlus  = 0;
    double waveFunctionOld   = Psi_value(pos_mat, alpha, beta);

    double kineticEnergy = 0;
    for(int i = 0; i < m_N; i++) {
        for(int j = 0; j < m_dim; j++) {
            pos_mat_plus[i][j] += h;
            pos_mat_minus[i][j] -= h;
            waveFunctionPlus = Psi_value(pos_mat_plus, alpha, beta);
            waveFunctionMinus = Psi_value(pos_mat_minus, alpha, beta);
            kineticEnergy += waveFunctionPlus + waveFunctionMinus -2*waveFunctionOld;
            pos_mat_plus[i][j] = pos_mat[i][j];
            pos_mat_minus[i][j] = pos_mat[i][j];
        }
    }
    kineticEnergy = -(0.5/(waveFunctionOld*h*h))*kineticEnergy;

    // Potential energy
    double potentialEnergy = 0;
    for(int i=0; i<m_N; i++) {
        potentialEnergy += V_ext(pos_mat[i], m_HO, omega_HO, omega_z);
    }
    return kineticEnergy + potentialEnergy;
}
