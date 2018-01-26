#pragma once
#include "vec3.h"


class WaveFunction
{
private:
    int                         m_dim = 1;
    int                         m_N = 1;
public:
    WaveFunction() {}
    int setTrialWF              (int dim, int N);
    double Psi_value(vec3 r_1, double alpha, double beta);
    double Psi_value_sqrd(vec3 r_1, double alpha, double beta);
    double E_L(vec3 r_1, double alpha, double omega_HO, double beta);
};
