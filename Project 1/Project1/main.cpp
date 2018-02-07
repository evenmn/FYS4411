#include <iostream>
#include <cmath>
#include <wavefunction.h>
#include <string.h>
#include <BF_met.h>
#include <vector>

using namespace std;

int main()
{
    //variables chosen by user
    double beta       = 1;       //weight parameter along z-axis
    double omega_HO   = 1;       //HO frequency in x- and y-direction
    double omega_z    = 1;       //HO frequency in z-direction
    int    M          = 100000;  //number of MC cycles
    double steplength = 1;       //steplength when changing position
    int    N          = 10;      //number of particles
    int    dim        = 3;       //number of dimensions concidered
    int    num_or_an  = 0;       //if calculation is to be based on analytical(0) or numerical(1) E_L
    bool   HO         = 0;       //spherical (true) or elliptical (false) harmonic oscillator
    double a          = 0;       //distance parameter
    double h          = 0.01;    // Step length when numerical double differentiating

    double alpha[]    = {0.25, 0.5, 1.0, 1.5};           //variational parameter
    int    len_alpha  = sizeof(alpha)/sizeof(*alpha);    //length of alpha

    Met_algo(N, dim, M, a, steplength, omega_HO, omega_z, HO, alpha, len_alpha, beta, h, num_or_an);

    return 0;
}
