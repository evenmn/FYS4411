#pragma once
#include "vec3.h"


class WaveFunction
{
private:
    int                         m_dim = 1;
    int                         m_N = 1;
    int                         m_n_or_a;
public:
    WaveFunction() {}
    int setTrialWF              (int dim, int N, int n_or_a);
    double Psi_value(double pos_mat[][3], double alpha, double beta);
    double Psi_value_sqrd(double pos_mat[][3], double alpha, double beta);
    double E_L(double pos_mat[][3], double alpha, double omega_HO, double beta);
};
