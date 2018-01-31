#include "wavefunction.h"
#include <cmath>
#include <iostream>

using namespace std;


int WaveFunction::setTrialWF(int dim, int N, double a, int n_or_a)
{

    m_dim = dim;        //number of dimensions
    m_N = N;            //number of particles
    m_n_or_a = n_or_a;  //analytical (0) or numerical (1) E_L calculation
    m_a = a;            //dimension of trap
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
            double rij = sqrt((pos_mat[i][0]-pos_mat[j][0])*(pos_mat[i][0]-pos_mat[j][0]) + \
                              (pos_mat[i][1]-pos_mat[j][1])*(pos_mat[i][1]-pos_mat[j][1]) + \
                              (pos_mat[i][2]-pos_mat[j][2])*(pos_mat[i][2]-pos_mat[j][2]));
            double f;

            if(rij>m_a){
                f = 1 - m_a/rij;
            }
            else{
                f = 0;
            }

            double u = log(f);
            sumu += u;
        }
    }
    //cout << sumu << " " << sumsqrt << endl;
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
            double rij = sqrt((pos_mat[i][0]-pos_mat[j][0])*(pos_mat[i][0]-pos_mat[j][0]) + \
                              (pos_mat[i][1]-pos_mat[j][1])*(pos_mat[i][1]-pos_mat[j][1]) + \
                              (pos_mat[i][2]-pos_mat[j][2])*(pos_mat[i][2]-pos_mat[j][2]));
            double f;

            if(rij>m_a){
                f = 1 - m_a/rij;
            }
            else{
                f = 0;
            }

            double u = log(f);
            sumu += u;
        }
    }
    //cout << sumu << " " << sumsqrt << endl;
    return exp(-alpha*sumsqrt*2)*exp(sumu*2);

}


double WaveFunction::E_L(double pos_mat[][3], double alpha, double omega_HO, double beta)
{
    if(m_n_or_a==0){

        if(m_a==0){
            //OBS! For a=0
            double sum_xixj = 0;
            double sum_yiyj = 0;
            double sum_zizj = 0;
            double sum_xyz_sqrd = 0;

            for(int i=0;i<m_N;i++){
                sum_xyz_sqrd += pos_mat[i][0]*pos_mat[i][0] + pos_mat[i][1]*pos_mat[i][1] + beta*pos_mat[i][2]*pos_mat[i][2];
                cout << sum_xyz_sqrd << endl;
                for(int j=0;j<m_N;j++){
                    sum_xixj += pos_mat[i][0]*pos_mat[j][0];
                    sum_yiyj += pos_mat[i][1]*pos_mat[j][1];
                    sum_zizj += pos_mat[i][2]*pos_mat[j][2];
                }
            }

            double func = -2*m_N*alpha*(sum_xixj + sum_yiyj + beta*beta*sum_zizj) + alpha*m_N*m_N + 0.5*omega_HO*omega_HO*sum_xyz_sqrd;

            if(m_dim==2){
                func += alpha*m_N*m_N;
            }
            else if(m_dim==3){
                func +=  2*alpha*m_N*m_N + beta*alpha*m_N*m_N;
            }
        }
    }

    else if(m_n_or_a==1){
        //Implement numerical solution
    }
}

