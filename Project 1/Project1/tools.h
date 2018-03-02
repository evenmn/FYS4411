#pragma once
#include <vector>

double QForce(double pos, double alpha, double beta, int dim_rand);

double GreenFuncSum(std::vector<std::vector<double>> &pos_mat, std::vector<std::vector<double>> &pos_mat_new, double D, double timestep, int N, double alpha, double beta);
