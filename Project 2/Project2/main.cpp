#include <iostream>
#include <random>
#include "eigen3/Eigen/Dense"
#include "gradient_descent.h"

using namespace std;
using namespace Eigen;


int main()
{
    int P = 2;               //Number of particles
    int D = 2;               //Number of dimensions
    int N = 2;               //Number of hidden nodes
    int MC = 100000;//pow(2,20);      //Number of Monte Carlo cycles
    int iterations = 100;      //Number of gradient decent cycles
    double sigma = 1.0;      //Width of Gaussian distribution
    double omega = 1.0;      //Frequency
    double steplength = 1.0; //Steplength for Metropolis
    double timestep = 0.01;  //Timestep used in Hastings algorithm
    double eta = 0.01;      //Learning rate for gradient decent
    bool interaction = 0;   //Interaction on if true
    double Diff = 0.5;         //Diffusion constant
    int sampling = 2;       //Brute force- (0), Hastings- (1) or Gibbs' sampling (2)

    GradientDescent(P, Diff, D, N, MC, iterations, sampling, sigma, omega, steplength, timestep, eta, interaction);
    return 0;
}


