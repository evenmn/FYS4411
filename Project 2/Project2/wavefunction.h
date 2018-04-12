#pragma once
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class WaveFunction
{
private:
    int m_M;
    int m_N;
    double m_sigma_sqrd;
    double m_omega_sqrd;
public:
    WaveFunction() {}
    int setTrialWF              (int N, int M, double sigma, double omega);
    double Psi_value_sqrd(VectorXd a, VectorXd b, VectorXd X, MatrixXd W);
    double EL_calc(VectorXd X, VectorXd a, VectorXd b, MatrixXd W);
};
