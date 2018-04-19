#pragma once
#include "eigen3/Eigen/Dense"

double QForce(Eigen::VectorXd X, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd W, int N, double sigma);
double GreenFuncSum(Eigen::VectorXd X, Eigen::VectorXd X_new, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd W, int N, double sigma, double timestep, int D, double Diff);
