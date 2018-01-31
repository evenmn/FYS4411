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
    //loop over several alphas
    double alpha = 1;           //variational parameter


    double beta = 1;            //weight parameter for z-axis
    double omega_HO = 1;        //frequency
    int M = 10;                 //number of MC cycles
    double steplength = 1;      //steplength when changing position
    int N = 1;                  //number of particles
    int dim = 3;                //number of dimensions concidered

    //averages and energies
    double E_tot = 0;
    double E_tot_sqrd= 0;
    double E = 0;
    double E_prev = 0;
    double delta_EL;

    double psi_ratio;           //ratio of new and old wave function
    double r;                   //random number
    int N_rand;                 //randomly chosen N
    int dim_rand;               //randomly chosen dimension
    double pos_mat_new[N][3]; //new position with random position

    //Initialize position matrix for N particles in dim dimentions
    double pos_mat [N][3];
    for(int i=0; i<N; i++) {
        for(int j=0; j<3; j++) {
            pos_mat[i][j] = random_position(1);
        }
    }

    WaveFunction Psi;
    Psi.setTrialWF(dim, N);

    //add initial energies to averages
    E = Psi.E_L(pos_mat, alpha, omega_HO, beta);

    E_tot += E;
    E_tot_sqrd += E*E;


    for(int i=0;i<M;i++){

        N_rand = rand()%N;
        dim_rand = 0;

        //Set new meatrix equal old one
        memcpy(pos_mat_new, pos_mat, sizeof(pos_mat_new));

        //Proposed new position
        pos_mat_new[N_rand][dim_rand] = random_position(steplength);


        psi_ratio = Psi.Psi_value_sqrd(pos_mat_new, alpha, beta)/(Psi.Psi_value_sqrd(pos_mat, alpha, beta));
        E_prev = Psi.E_L(pos_mat, alpha, omega_HO, beta);

        r = ((double)rand() / (double)RAND_MAX);

        if(psi_ratio >= r){
            //accept
            memcpy(pos_mat, pos_mat_new, sizeof(pos_mat));
        }

        E = Psi.E_L(pos_mat, alpha, omega_HO, beta);
        delta_EL = E - E_prev;
        E_tot += E;
        E_tot_sqrd += E*E;

    }

    //Calculate <E_l> and <E_L**2>
    double E_L_avg = E_tot/M;
    double E_L_avg_sqrd = E_tot_sqrd/M;

    cout << "E_L_avg: " << E_L_avg << endl;
    cout << "E_L_ang_tot: " << E_L_avg_sqrd << endl;


    return 0;
}

