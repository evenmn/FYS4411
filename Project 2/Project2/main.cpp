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
    double sigma = 1.0;  //Width of Gaussian distribution
    double omega = 1.0;  //Frequency

    GradientDescent(P, D, N, sigma, omega);
    return 0;
}


