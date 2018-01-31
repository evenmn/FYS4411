#include "wavefunction.h"
#include <cmath>
#include <iostream>

using namespace std;


int WaveFunction::setTrialWF(int dim, int N)
{

    m_dim = dim;        //number of dimensions
    m_N = N;            //number of particles
}

double WaveFunction::Psi_value(double pos_mat[][3],  double alpha, double beta)
{

    if(m_dim==1&&m_N==1){
        double x = pos_mat[0][0];
        return exp(-alpha*x*x);
    }

    else if(m_dim==3){
        /* Returns the wavefunction based on an array of
        positions preferably in three dimensions */

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
                sumu += 1/rij;
            }
        }
        //cout << sumu << " " << sumsqrt << endl;
        return exp(-alpha*sumsqrt)*exp(sumu);
    }

}

double WaveFunction::Psi_value_sqrd(double pos_mat[][3], double alpha, double beta)
{

    if(m_dim==1&&m_N==1){
        double x = pos_mat[0][0];
        return exp(-alpha*x*x*2);
    }

}


double WaveFunction::E_L(double pos_mat[][3], double alpha, double omega_HO, double beta)
{

    if(m_dim==1&&m_N==1){
        double x = pos_mat[0][0];
        return -2*x*x*alpha*alpha + alpha + 0.5*omega_HO*omega_HO*x*x;
    }

}

