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

double V_ext(vector<double> &pos, double beta, int dim) {
    double VEXT = 0;
    for(int i=0;i<dim;i++){
        if(i==2){
            VEXT += beta*beta*pos[i]*pos[i];
        }
        else{
            VEXT += pos[i]*pos[i];
        }
    }
    return 0.5*VEXT;
}


int WaveFunction::setTrialWF(int dim, int N, double a, int n_or_a)
{

    m_dim    = dim;     //number of dimensions
    m_N      = N;       //number of particles
    m_n_or_a = n_or_a;  //analytical (0) or numerical (1) E_L calculation
    m_a      = a;       //dimension of trap
}

double WaveFunction::Psi_value(vector<vector<double>> &pos_mat, double alpha, double beta)
{
    //Returns the wavefunction squared
    double sumsqrt = 0;         // x1^2 + y1^2 + B*z1^2 + ... + B*zn^2
    double prodf = 1;
    double norm = 0;
    bool leave = false;

    for(int i=0; i<m_N; i++) {
        for(int k=0; k<m_dim; k++){
            if(k==2){
                sumsqrt += beta*pos_mat[i][k]*pos_mat[i][k];
            }
            else{
                sumsqrt += pos_mat[i][k]*pos_mat[i][k];
            }

            norm += pos_mat[i][k]*pos_mat[i][k];
        }
        norm = sqrt(norm);
        for(int j=m_N-1; j>i; j--) {
            if(norm > m_a){
                prodf *= 1 - m_a/norm;
            }
            else{
                prodf = 0;
                leave = true;
                break;
            }
        }
        if(leave)break;
    }
    return exp(-alpha*sumsqrt) * prodf;

}

double WaveFunction::Psi_value_sqrd(vector<vector<double>> &pos_mat, double alpha, double beta)
{
    //Returns the wavefunction squared
    double sumsqrt = 0;         // x1^2 + y1^2 + B*z1^2 + ... + B*zn^2
    double prodf = 1;
    double norm = 0;
    bool leave = false;

    for(int i=0; i<m_N; i++) {
        for(int k=0; k<m_dim; k++){
            if(k==2){
                sumsqrt += beta*pos_mat[i][k]*pos_mat[i][k];
            }
            else{
                sumsqrt += pos_mat[i][k]*pos_mat[i][k];
            }

            norm += pos_mat[i][k]*pos_mat[i][k];
        }
        norm = sqrt(norm);
        for(int j=m_N-1; j>i; j--) {
            if(norm > m_a){
                prodf *= 1 - m_a/norm;
            }
            else{
                prodf = 0;
                break;
                leave = true;
            }
        }
        if(leave)break;
    }
    return exp(-alpha*sumsqrt*2) * prodf * prodf;
}


double WaveFunction::E_L_ana(vector<vector<double>> &pos_mat, double alpha, double beta)
{
    double distij=0;
    double distik=0;
    double r_ij[3];
    double r_ik[3];
    double EL;
    double E_TOT = 0;

    if(m_a==0){
        double r_sqrt = 0;
        double energy = 0;
        for(int i=0; i<m_N; i++){
            for(int j=0; j<m_dim; j++){
                if(j == 2) {
                    r_sqrt += beta*beta*pos_mat[i][j]*pos_mat[i][j];
                    energy += alpha*beta;
                }
                else{
                    r_sqrt += pos_mat[i][j]*pos_mat[i][j];
                    energy += alpha;
                }
            }
        }
        return (0.5 - 2*alpha*alpha)*r_sqrt + energy;
    }

    else{
        for(int i=0; i<m_N; i++) {
            EL = 0;

            for(int j=0; j<m_dim; j++){
                if(j==2){
                    EL += beta*beta*pos_mat[i][j]*pos_mat[i][j] - beta/(2*alpha);
                }
                else{
                    EL += pos_mat[i][j]*pos_mat[i][j] - 1/(2*alpha);
                }
            }

            EL = 4*alpha*alpha*EL;

            double term1 = 0;
            for(int j=m_N-1; j>i; j--) {
                for(int l = 0; l < m_dim; l++){
                    r_ij[l] = pos_mat[i][l] - pos_mat[j][l];
                    distij += r_ij[l]*r_ij[l];
                    if(l==2){
                        term1 += pos_mat[i][l] * r_ij[l] * beta;
                    }
                    else{
                        term1 += pos_mat[i][l] * r_ij[l];
                    }


                }
                distij = sqrt(distij);

                double u_der_ij = u_der(distij, m_a);
                EL -= 4*alpha*term1*u_der_ij;

                double term2 = 0;
                for(int k=m_N-1; k>i; k--) {
                    for(int l = 0; l < m_dim; l++){
                        r_ik[l] = pos_mat[i][l] - pos_mat[k][l];
                        distik += r_ik[l]*r_ik[l];
                        term2 += r_ik[l]*r_ij[l];
                    }
                    distik = sqrt(distik);

                    EL += term2*u_der_ij * u_der(distik, m_a);
                }

                EL += u_secder(distij, m_a) + 2 * u_der_ij;
            }

            E_TOT += -0.5*EL + V_ext(pos_mat[i], beta, m_dim);
        }
    }

    return E_TOT;
}

double WaveFunction::E_L_num(vector<vector<double>> &pos_mat, double alpha, double beta, double h)
{
    //Obs:momentarily only work with a=0
    // Kinetic energy
    vector<vector<double>> pos_mat_plus;
    vector<vector<double>> pos_mat_minus;

    pos_mat_plus = pos_mat;
    pos_mat_minus = pos_mat;

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
        potentialEnergy += V_ext(pos_mat[i], beta, m_dim);
    }
    return kineticEnergy + potentialEnergy;
}

double WaveFunction::Psi_der(vector<vector<double>> &pos_mat, double beta) {

    double psi_der_div_psi = 0;

    for(int i=0; i<m_N; i++){
        for(int j=0;j<m_dim;j++){
            if(j==2){
                psi_der_div_psi += beta*pos_mat[i][j]*pos_mat[i][j];
            }
            else{
                psi_der_div_psi += pos_mat[i][j]*pos_mat[i][j];
            }
        }
    }
    return -psi_der_div_psi;
}


