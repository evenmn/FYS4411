#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>
#include "gibbs_tools.h"

using namespace std;

random_device RD;                   //Will be used to obtain a seed for the random number engine
mt19937 Gen(RD());                  //Standard mersenne_twister_engine seeded with rd()

double x_sampling(const VectorXd &a, const VectorXd &h, const MatrixXd &W, double sigma, int i) {
    double mu = a(i);
    int N = h.size();
    for(int j=0; j<N; j++) {
        mu += h(j)*W(i,j);
    }
    normal_distribution<double> d(mu, sigma * sigma);
    return d(Gen);
}
