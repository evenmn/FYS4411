#include <iostream>
#include <cmath>
#include <vec3.h>
#include <wavefunction.h>
#include <string.h>

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
    int M = 10000;               //number of MC cycles
    double steplength = 1;      //steplength when changing position
    int N = 1;                 //number of particles
    int dim = 1;                //number of dimensions concidered
    int num_or_an = 0;          //if calculation is to be based on analytical(0) or numerical(1) E_L
    int HO = 0;                 //spherical (0) or elliptical (1) harmonic oscillator
    double a = 0; //0.4*pow(10,-10); //distance par        if(psi_ratio > 1) {

    //loop over several alphas
    double alpha[] = {0.1, 0.25, 0.5, 0.75, 1.0};           //variational parameter
    int length_alpha = sizeof(alpha)/sizeof(*alpha);

    //Everything below here could be put in own script
    for(int k=0; k<length_alpha;k++){
        //averages and energies
        double E_tot = 0;           //sum of energies of all states
        double E_tot_sqrd = 0;       //sum of energies of all states squared
        double E = 0;               //energy after change in position
        double E_prev = 0;          //energy before change in position
        double delta_EL;            //change in energy

        double psi_ratio;           //ratio of new and old wave function
        double r;                   //random number
        int N_rand;                 //randomly chosen N
        int dim_rand;               //randomly chosen dimension
        double pos_mat_new[N][3];   //new position with random position

        //Initialize position matrix for N particles in dim dimentions
        double pos_mat[N][3];
        for(int i=0; i<N; i++) {
            for(int j=0; j<dim; j++) {
                pos_mat[i][j] = random_position(1);
            }
            for(int k=dim; k<3; k++){
                pos_mat[i][k] = 0;
            }
        }


        //Initialize wave function
        WaveFunction Psi;
        Psi.setTrialWF(dim, N, a, num_or_an, HO);

        //Add initial energies to averages
        E = Psi.E_L(pos_mat, alpha[k], beta, omega_HO, omega_z);
        E_tot += E;
        E_tot_sqrd += E*E;

        //Start Monte Carlo iterations
        for(int i=0;i<M;i++){
            //Draw random position, for one particle and one dimention
            N_rand = rand()%N;
            dim_rand = rand()%dim;

            //Set new meatrix equal old one
            memcpy(pos_mat_new, pos_mat, sizeof(pos_mat_new));

            //Proposed new position
            pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + (random_position(1)-0.5)*steplength;

            //Metropolis algorithm
            psi_ratio = Psi.Psi_value_sqrd(pos_mat_new, alpha[k], beta)/(Psi.Psi_value_sqrd(pos_mat, alpha[k], beta));
            E_prev = Psi.E_L(pos_mat, alpha[k], beta, omega_HO, omega_z);

            r = ((double)rand() / (double)RAND_MAX);

            if(psi_ratio >= r){
                //accept and update pos_mat
                memcpy(pos_mat, pos_mat_new, sizeof(pos_mat)); //maybe more time efficient to only update the one changed position?
            }

            E = Psi.E_L(pos_mat, alpha[k], beta, omega_HO, omega_z);
            delta_EL = E - E_prev;
            E_tot += E;
            E_tot_sqrd += E*E;
        }

        //Calculate <E_l> and <E_L**2>
        double E_L_avg = E_tot/M;
        double E_L_avg_sqrd = E_tot_sqrd/M;

    cout << "E_L_avg: " << E_L_avg << endl;
    cout << "E_L_avg_tot: " << E_L_avg_sqrd << "\n" << endl;
    }

    return 0;

}

