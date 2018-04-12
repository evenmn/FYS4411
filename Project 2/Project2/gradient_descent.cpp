#include <iostream>
#include <wavefunction.h>
#include <random>
#include <eigen3/Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace std;

//Mersenne Twister RNG
random_device rd;                   //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0, 1);

double random_position(){
    return dis(gen);
}

void GradientDescent(int P, int D, int N, double sigma, double omega) {

    double psi_ratio;               //ratio of new and old wave function
    int    N_rand;                  //randomly chosen N
    int    dim_rand;                //randomly chosen dimension
    int M = P*D;

    for(int i=0; i<1; i++) {
        //averages and energies
        double EL_tot      = 0;          //sum of energies of all states
        double EL_tot_sqrd = 0;          //sum of energies of all states squared
        double EL          = 0;          //energy after change in position

        MatrixXd W = MatrixXd::Random(M,N);     //Weights
        VectorXd a = VectorXd::Random(M);       //Visible biases
        VectorXd b = VectorXd::Random(N);       //Hidden biases
        VectorXd X = VectorXd::Random(M);       //Visible nodes (position coordinates)
        VectorXd X_new;

        cout << W << endl;
    }
}
