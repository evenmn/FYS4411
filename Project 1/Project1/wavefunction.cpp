#include "wavefunction.h"
#include <cmath>
#include <iostream>

using namespace std;

double rij(double pos1[3], double pos2[3]) {
    return sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0]) + (pos1[1]-pos2[1])*(pos1[1]-pos2[1]) + (pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
}

double u_der(double dist, double a) {
    if(dist > a) {
        return (a/(dist*dist))/(1-a/dist);
    }
    else {
        return 0;
    }
}

double u_secder(double dist, double a) {
    double C;
    if(dist > a) {
        C = (1/(1-a/dist));
        return -(a/pow(dist, 3))*C*(2+(a/dist)*C);
    }
    else {
        return 0;
    }
}

double V_ext(double pos[3], int HO, double omega_HO, double omega_z) {
    if(HO == 0) {
        return 0.5*omega_HO*omega_HO*(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    }
    else if(HO == 1) {
        return 0.5*(omega_HO*omega_HO*(pos[0]*pos[0] + pos[1]*pos[1]) + omega_z*omega_z*pos[2]*pos[2]);
    }
    else {
        cout << "HO has unexpected value" << endl;
        return false;
    }
}


int WaveFunction::setTrialWF(int dim, int N, double a, int n_or_a, int HO)
{

    m_dim = dim;        //number of dimensions
    m_N = N;            //number of particles
    m_n_or_a = n_or_a;  //analytical (0) or numerical (1) E_L calculation
    m_a = a;            //dimension of trap
    m_HO = HO;          //spherical (0) or elliptical (1) harmonic oscillator
}

double WaveFunction::Psi_value(double pos_mat[][3], double alpha, double beta)
{

    //Returns the wavefunction
    double sumsqrt = 0;         // x1^2 + y1^2 + B*z1^2 + ... + B*zn^2
    double sumu = 0;            // sum(u(rij))

    for(int i=0; i<m_N; i++) {
        sumsqrt += pos_mat[i][0]*pos_mat[i][0] + \
                   pos_mat[i][1]*pos_mat[i][1] + \
                   beta*pos_mat[i][2]*pos_mat[i][2];

        for(int j=m_N; j>i; j--) {
            double norm = rij(pos_mat[i], pos_mat[j]);
            double f;

            if(norm > m_a){
                f = 1 - m_a/norm;
            }
            else{
                f = 0;
            }
            sumu += log(f);
        }
    }
    return exp(-alpha*sumsqrt)*exp(sumu);

}

double WaveFunction::Psi_value_sqrd(double pos_mat[][3], double alpha, double beta)
{
    //Returns the wavefunction
    double sumsqrt = 0;         // x1^2 + y1^2 + B*z1^2 + ... + B*zn^2
    double sumu = 0;            // sum(u(rij))

    for(int i=0; i<m_N; i++) {
        sumsqrt += pos_mat[i][0]*pos_mat[i][0] + \
                   pos_mat[i][1]*pos_mat[i][1] + \
                   beta*pos_mat[i][2]*pos_mat[i][2];

        for(int j=m_N; j>i; j--) {
            double norm = rij(pos_mat[i], pos_mat[j]);
            double f;

            if(norm > m_a){
                f = 1 - m_a/norm;
            }
            else{
                f = 0;
            }
            sumu += log(f);
        }
    }
    return exp(-alpha*sumsqrt*2)*exp(sumu*2);
}


double WaveFunction::E_L(double pos_mat[][3], double alpha, double beta, double omega_HO, double omega_z)
{
    double distij;
    double distik;
    double EL;
    double E_TOT = 0;

    // This should work for all a's and in all dimensions
    for(int i=0; i<m_N; i++) {
        EL = 0;
        EL += 4*alpha*alpha*(pos_mat[i][0]*pos_mat[i][0] + pos_mat[i][1]*pos_mat[i][1] + \
              beta*beta*pos_mat[i][2]*pos_mat[i][2]);

        if(m_dim==1) {
            EL += -2*alpha;
        }
        else if(m_dim==2) {
            EL += -4*alpha;
        }

        else if(m_dim==3) {
            EL += -4*alpha -2*alpha*beta;
        }

        for(int j=m_N-1; j>i; j--) {
            distij = rij(pos_mat[i], pos_mat[j]);
            EL += -4*alpha*(pos_mat[i][0] * (pos_mat[i][0]-pos_mat[j][0]) +\
                            pos_mat[i][1] * (pos_mat[i][1]-pos_mat[j][1]) +\
                            pos_mat[i][2] * (pos_mat[i][2]-pos_mat[j][2]) * beta) * \
                  (u_der(distij, m_a)/distij);

            for(int k=m_N-1; k>i; k--) {
                distik = rij(pos_mat[i], pos_mat[k]);
                EL += ((pos_mat[i][0] - pos_mat[k][0]) * (pos_mat[i][0] - pos_mat[j][0]) + \
                       (pos_mat[i][1] - pos_mat[k][1]) * (pos_mat[i][1] - pos_mat[j][1]) + \
                       (pos_mat[i][2] - pos_mat[k][2]) * (pos_mat[i][2] - pos_mat[j][2])) * \
                      u_der(distij, m_a) * u_der(distik, m_a) / (distij * distik);
            }

            EL += u_secder(distij, m_a) + (2/distij) * u_der(distij, m_a);
        }
        E_TOT += -0.5*EL + V_ext(pos_mat[i], m_HO, omega_HO, omega_z);
    }
    return E_TOT;
}

double WaveFunction::E_L_num(double pos_mat[][3], double alpha, double beta, double h)
{
    // Kinetic energy
    double psi_old = Psi_value(pos_mat, alpha, beta);

    double pos_mat_plus[m_N][3];
    double pos_mat_minus[m_N][3];

    for(int i=0; i<m_N; i++) {
        for(int j=0; j<m_dim; j++) {
            pos_mat_plus[i][j] = pos_mat[i][j] + h;
            pos_mat_minus[i][j] = pos_mat[i][j] - h;
        }
    }

    double sec_der_psi = (Psi_value(pos_mat_plus, alpha, beta) -2*psi_old + Psi_value(pos_mat_minus, alpha, beta))/(h*h);
    double E_k = -0.5*(sec_der_psi/psi_old);

    // Potential energy
    double E_p = 0;
    for(int i=0; i<m_N; i++) {
        for(int j=m_N-1; j>i; j--) {
            E_p += rij(pos_mat[i], pos_mat[j]);
        }
    }
    return E_k + E_p;
}
