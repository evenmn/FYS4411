#pragma once
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class WaveFunction
{
private:
    int m_M;
    int m_N;
    int m_sampling;
    double m_sigma;
    double m_sigma_sqrd;
    double m_omega_sqrd;
public:
    WaveFunction() {}
    int setTrialWF              (int N, int M, int sampling, double sigma_sqrd, double omega);
    double Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v);
    double Psi_value_sqrd_hastings(VectorXd Xa, VectorXd v);
    double EL_calc(VectorXd X, VectorXd Xa, VectorXd v, MatrixXd W, int D, int interaction, double &E_kin, double &E_ext, double &E_int);
    void Gradient_a(const VectorXd &a, VectorXd &da);
    void Gradient_b(const VectorXd &e, VectorXd &db);
    void Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW);
    double Gradient_sigma(const MatrixXd &W, const VectorXd &Xa, const VectorXd &X, const VectorXd &e);
};
