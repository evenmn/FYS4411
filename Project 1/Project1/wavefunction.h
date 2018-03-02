#pragma once

class WaveFunction
{
private:
    int                         m_dim;
    int                         m_N;
    double                      m_a;
    int                         m_n_or_a;
    int                         m_HO;
public:
    WaveFunction() {}
    int setTrialWF              (int dim, int N, double a, int n_or_a, bool HO);
    double Psi_value(double pos_mat[][3], double alpha, double beta);
    double Psi_value_sqrd(double pos_mat[][3], double alpha, double beta);
    double E_L_ana(double pos_mat[][3], double alpha, double beta, double omega_HO);
    double E_L_num(double pos_mat[][3], double alpha, double beta, double omega_HO, double h);
    double Psi_der(double pos_mat[][3], double alpha, double beta);
};
