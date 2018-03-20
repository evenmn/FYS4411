#include <iostream>
#include <metropolis.h>
#include <gd.h>

using namespace std;

int main()
{
    //variables chosen by user
    int    M          = 1000000; //number of MC cycles
    int    N          = 10;      //number of particles
    int    dim        = 3;       //number of dimensions concidered
    double beta       = 2.82843;       //weight parameter along z-axis
    double steplength = 1.0;     //steplength when changing position
    double h          = 0.01;    //Step length for numerical double differentiation
    double timestep   = 1;     //Timestep, to be used in Metropolis-Hastings
    double a          = 0;     //distance parameter

    //Choices of simulation types
    int    BF_H       = 0;       //brute force (0) or hastings(1) metropolis algorithm
    int    one_body   = 1;       //calculate one body density (1) or not calculate (0)
    int    num_or_an  = 0;       //if calculation is to be based on analytical(0) or numerical(1) E_L


    double alpha[]    = {0.5};           //variational parameter
    int    len_alpha  = sizeof(alpha)/sizeof(*alpha);    //length of alpha

    cout << "Running with the following parameters:" << endl;
    cout << "Number of particles: " << N << endl;
    cout << "Number of dimensions: " << dim << endl;
    cout << "a, distance parameter: " << a << "\n" << endl;
    cout << "Analytical(0) or numerical(1) E_L calculation: " << num_or_an << endl;
    cout << "Brute force(0) or Hastings(1) Metropolis algo: " << BF_H << endl;
    cout << "One body calculations active (1): " << one_body << "\n" << endl;

    Met_algo(N, dim, M, a, steplength, alpha, len_alpha, beta, h, num_or_an, BF_H, timestep, one_body);
    //GradientDecent(N, dim, M, a, steplength, beta, h, num_or_an, BF_H, timestep);
    return 0;
}
