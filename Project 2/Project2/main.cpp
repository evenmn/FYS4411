#include <iostream>
#include "eigen3/Eigen/Dense"
#include <random>
#include <gradient_descent.h>

using namespace std;
using namespace Eigen;


int main()
{
    int P = 2;           //Number of particles
    int D = 2;           //Number of diemnsions
    int N = 2;           //Number of hidden nodes
    int M = P*D;         //Number of visible nodes, free dimensions
    int MC = 1;          //Number of Monte Carlo cycles
    double sigma = 1.0;  //Width of Gaussian distribution
    double omega = 1.0;  //Frequency
    double steplength = 1.0;

    GradientDescent(P, D, N, MC, sigma, omega, steplength);
    return 0;
}


