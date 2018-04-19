#pragma once
#include "eigen3/Eigen/Dense"

double QForce(const Eigen::VectorXd &X, const Eigen::VectorXd &a, const Eigen::VectorXd &b, const Eigen::MatrixXd &W, int N, double sigma, int i);
double GreenFuncSum(const Eigen::VectorXd &X, const Eigen::VectorXd &X_new, const Eigen::VectorXd &a, const Eigen::VectorXd &b, const Eigen::MatrixXd &W, int N, double sigma, double timestep, int D, double Diff);
