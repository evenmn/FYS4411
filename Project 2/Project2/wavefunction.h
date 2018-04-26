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
    int setTrialWF              (int N, int M, double sigma_sqrd, double omega);
    double Psi_value_sqrd(VectorXd a, VectorXd b, VectorXd X, MatrixXd W);
    double Psi_value_sqrd_hastings(VectorXd Xa, VectorXd v);
    double EL_calc(VectorXd X, VectorXd Xa, VectorXd v, MatrixXd W, int D, int interaction, double &E_k, double &E_ext, double &E_int, int sampling);
    void Gradient_a(const VectorXd &a, VectorXd &da, int sampling);
    void Gradient_b(const VectorXd &b, VectorXd &db, int sampling);
    void Gradient_W(const VectorXd &X, const VectorXd &v, MatrixXd &dW, int sampling);
};
