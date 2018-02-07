#include <iostream>
#include <cmath>
#include <vec3.h>
#include <wavefunction.h>
#include <string.h>
#include <vector>
#include <BF_met.h>

using namespace std;

double random_position(double steplength){
    return ((double)rand() / (double)RAND_MAX)*steplength;
}

int main()
{
    //variables chosen by user
    double beta = 1;            //weight parameter along z-axis
    double omega_HO = 1;        //HO frequency in x- and y-direction
    double omega_z = 1;         //HO frequency in z-direction
    int M = 100000;               //number of MC cycles
    double steplength = 1;      //steplength when changing position
    int N = 4;                 //number of particles
    int dim = 3;                //number of dimensions concidered
    int num_or_an = 0;          //if calculation is to be based on analytical(0) or numerical(1) E_L
    int HO = 0;                 //spherical (0) or elliptical (1) harmonic oscillator
    double a = 0;               //distance parameter
    double alpha[] = {0.25, 0.5, 1.0, 1.5};           //variational parameter
    int length_alpha_1 = sizeof(alpha)/sizeof(*alpha);
    double h = 0.01;            // Step length when numerical double differentiating

    /*
    vector<vector<double>> posmat;
    posmat.resize(N);
    for(int i = 0; i < posmat.size(); i++)
        posmat[i].resize(dim);

    for(auto* i : posmat)
        i.resize(dim);
   */

    Met_algo(N, dim, M, a, steplength, omega_HO, omega_z, HO, alpha, length_alpha_1, beta, h);

    return 0;
}
