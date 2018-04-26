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
uniform_real_distribution<> Dis(0, 1);

double Random_position(){
    return Dis(Gen);
}

double x_sampling(const VectorXd &a, const VectorXd &h, const MatrixXd &W, double sigma_sqrd, int i) {
    double mu = a(i);
    int N = h.size();
    for(int j=0; j<N; j++) {
        mu += h(j)*W(i,j);
    }
    normal_distribution<double> d(mu, sigma_sqrd/2);
    return d(Gen);

}

double h_sampling(const VectorXd &b, const VectorXd &X, const MatrixXd &W, double sigma_sqrd, int i){
    //unsure how to implement; should the largest of P(h=1 given x) and P(h=0 given x be chosen?)

    VectorXd v = b + (X.transpose() * W).transpose()/sigma_sqrd;

    double P_h1 = 1 + exp(-2*v(i));
    P_h1 = 1.0/P_h1;


    double P_h0 = 1 + exp(2*v(i));
    P_h0 = 1.0/P_h0;

    if(P_h1>= Random_position()){
        return 1;
    }

    else{
        return 0;
    }

}

