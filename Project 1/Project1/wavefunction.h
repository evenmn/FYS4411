#pragma once
#include <vector>

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
    int setTrialWF              (int dim, int N, double a, int n_or_a);
    double Psi_value(std::vector<std::vector<double>> &pos_mat, double alpha, double beta);
    double Psi_value_sqrd(std::vector<std::vector<double>> &pos_mat, double alpha, double beta);
    double E_L_ana(std::vector<std::vector<double>> &pos_mat, double alpha, double beta);
    double E_L_num(std::vector<std::vector<double>> &pos_mat, double alpha, double beta, double h);
    double Psi_der(std::vector<std::vector<double>> &pos_mat, double beta);
};
