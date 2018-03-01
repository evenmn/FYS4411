#pragma once

double QForce(double pos, double alpha, double beta, int dim_rand);

double GreenFuncSum(double pos_mat[][3], double pos_mat_new[][3], double D, double timestep, int N, double alpha, double beta);
