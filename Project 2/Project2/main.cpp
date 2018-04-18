#include <iostream>
#include "eigen3/Eigen/Dense"
#include <random>
#include <gradient_descent.h>

using namespace std;
using namespace Eigen;


int main()
{
    int P = 2;               //Number of particles
    int D = 2;               //Number of dimensions
    int N = 2;               //Number of hidden nodes
    int MC = 100000;        //Number of Monte Carlo cycles
    int iterations = 100;      //Number of gradient decent cycles
    double sigma = 1.0;      //Width of Gaussian distribution
    double omega = 1.0;      //Frequency
    double steplength = 1.0; //Steplength for Metropolis
    double eta = 0.01;      //Learning rate for gradient decent

    GradientDescent(P, D, N, MC, iterations, sigma, omega, steplength, eta);
    return 0;
}


