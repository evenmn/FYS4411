#include <iostream>
#include <random>
#include "eigen3/Eigen/Dense"
#include "gradient_descent.h"

using namespace std;
using namespace Eigen;


int main()
{
    int     P           = 2;               //Number of particles
    int     D           = 2;               //Number of dimensions
    int     N           = 2;              //Number of hidden nodes
    int     MC          = 500000;        //Number of Monte Carlo cycles
    int     iterations  = 500;             //Number of gradient decent cycles
    int     sampling    = 2;               //Brute force- (0), Hastings- (1) or Gibbs' sampling (2)
    bool    interaction = 1;               //Interaction on if true
    bool    one_body    = 1;               //Calculating onebody density if true
    double  sigma       = 1.0;             //Width of Gaussian distribution
    double  omega       = 1.0;             //Frequency
    double  steplength  = 1.0;             //Steplength for Metropolis
    double  timestep    = 0.1;             //Timestep used in Hastings algorithm
    double  eta         = 0.01;            //Learning rate for gradient decent
    double  Diff        = 0.5;             //Diffusion constant

    GradientDescent(P, Diff, D, N, MC, iterations, sampling, sigma, omega, steplength, timestep, eta, interaction, one_body);
    return 0;
}


