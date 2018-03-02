#include <iostream>
#include <metropolis.h>
#include <gd.h>

using namespace std;

int main()
{
    //variables chosen by user
    double beta       = 1;       //weight parameter along z-axis
    double omega_HO   = 1;       //HO frequency in x- and y-direction
    int    M          = 1000000;  //number of MC cycles
    double steplength = 0.1;     //steplength when changing position ->Is this correct?
    int    N          = 10;      //number of particles
    int    dim        = 3;       //number of dimensions concidered
    int    num_or_an  = 0;       //if calculation is to be based on analytical(0) or numerical(1) E_L
    bool   HO         = true;    //spherical (true) or elliptical (false) harmonic oscillator
    int    BF_H       = 0;       //brute force (0) or hastings(1) metropolis algorithm
    int    one_body   = 0;       //calculate one body density (1) or not calculate (0)
    double a          = 0.0; //distance parameter
    double h          = 0.01;    //Step length for numerical double differentiation
    double timestep   = 1.0;    //Timestep, to be used in numerical derivatives

    double alpha[]    = {0.4,0.5,0.6};           //variational parameter
    int    len_alpha  = sizeof(alpha)/sizeof(*alpha);    //length of alpha

    cout << "Running with the following paramteres:" << endl;
    cout << "Number of particles: " << N << endl;
    cout << "Number of dimensions: " << dim << endl;
    cout << "a, distance parameter:" << a << "\n"<< endl;
    cout << "Analytical(0) or numerical(1) E_L calcualtion: " << num_or_an << endl;
    cout << "Spherical(true) or elliptical(false) harmonic oscillator: " << HO << endl;
    cout << "Brute force(0) or Hastings(1) Metropolis algo:" << BF_H << endl;
    cout << "One body calculations active (1): " << one_body << "\n" << endl;

    //Met_algo(N, dim, M, a, steplength, omega_HO, HO, alpha, len_alpha, beta, h, num_or_an, BF_H, timestep, one_body);
    GradientDecent(N, dim, M, a, steplength, omega_HO, HO, beta, h, num_or_an, BF_H, timestep);

    return 0;
}
