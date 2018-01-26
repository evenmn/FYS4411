#include "wavefunction.h"
#include <cmath>
#include <iostream>

using namespace std;


int WaveFunction::setTrialWF(int dim, int N)
{

    m_dim = dim;        //number of dimensions
    m_N = N;            //number of particles
}

double WaveFunction::Psi_value(vec3 r_1, double alpha, double beta)
{

    if(m_dim==1&&m_N==1){
        double x = r_1[0];
        return exp(-alpha*x*x);
    }

    else{
        //Here, implement E_L for N particles and 3 dimensions -> can work for any m_dim and N_dim
        //Cannot use vec3... must implement Nx3 matrix, that contains the positions of all the N particles

    }

}

double WaveFunction::Psi_value_sqrd(vec3 r_1, double alpha, double beta)
{

    if(m_dim==1&&m_N==1){
        double x = r_1.length();
        return exp(-alpha*x*x*2);
    }

}


double WaveFunction::E_L(vec3 r_1, double alpha, double omega_HO, double beta)
{

    if(m_dim==1&&m_N==1){
        double x = r_1.length();
        return -2*x*x*alpha*alpha + alpha + 0.5*omega_HO*omega_HO*x*x;
    }

}

