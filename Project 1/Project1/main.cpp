#include <iostream>
#include <metropolis.h>
#include <gd.h>

using namespace std;

int main(int argc, char* argv[])
{
    //variables chosen by user
    int    M          = 1000000; //number of MC cycles
    int    N          = 10;      //number of particles
    int    dim        = 3;       //number of dimensions concidered
    double beta       = 2.82843;       //weight parameter along z-axis
    double steplength = 1.0;     //steplength when changing position
    double h          = 0.01;    //Step length for numerical double differentiation
    double timestep   = 1;     //Timestep, to be used in Metropolis-Hastings
    double a          = 0.0043;     //distance parameter

    //Choices of simulation types
    bool   HO         = false;    //spherical (true) or elliptical (false) harmonic oscillator
    int    BF_H       = 0;       //brute force (0) or hastings(1) metropolis algorithm
    int    one_body   = 0;       //calculate one body density (1) or not calculate (0)
    int    num_or_an  = 0;       //if calculation is to be based on analytical(0) or numerical(1) E_L


    double alpha[]    = {0.350636};           //variational parameter
    int    len_alpha  = sizeof(alpha)/sizeof(*alpha);    //length of alpha

    cout << "Running with the following paramteres:" << endl;
    cout << "Number of particles: " << N << endl;
    cout << "Number of dimensions: " << dim << endl;
    cout << "a, distance parameter:" << a << "\n" << endl;
    cout << "Analytical(0) or numerical(1) E_L calcualtion: " << num_or_an << endl;
    cout << "Spherical(true) or elliptical(false) harmonic oscillator: " << HO << endl;
    cout << "Brute force(0) or Hastings(1) Metropolis algo:" << BF_H << endl;
    cout << "One body calculations active (1): " << one_body << "\n" << endl;

    //double alpha_opt = GradientDecent(N, dim, M, a, steplength, HO, beta, h, num_or_an, BF_H, timestep);
    //cout << "Optimal alpha: " << alpha_opt << "\n" << endl;

    //int M_met = 1000000;
    //double alpha[]    = {alpha_opt};           //variational parameter
    //int    len_alpha  = sizeof(alpha)/sizeof(*alpha);    //length of alpha
    Met_algo(N, dim, M, a, steplength, HO, alpha, len_alpha, beta, h, num_or_an, BF_H, timestep, one_body);
    //GradientDecent(N, dim, M, a, steplength, HO, beta, h, num_or_an, BF_H, timestep);

    return 0;
}
