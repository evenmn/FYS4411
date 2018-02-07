#include <BF_met.h>
#include<iostream>
#include<cmath>
#include<wavefunction.h>
#include <string.h>
#include<vector>

using namespace std;

double random_position_1(double steplength){
    return ((double)rand() / (double)RAND_MAX)*steplength;
}

void BF_calc(int N, int dim, int M, double a, double steplength, double omega_HO, double omega_z, int HO, double alpha[], int length_alpha_1, double beta)
{
    cout << length_alpha_1 << endl;
    cout << alpha[0] << " " << alpha[1] << " "<< alpha[2] << endl;

    for(int k=0; k<length_alpha_1; k++){
        //averages and energies
        double E_tot_1 = 0;           //sum of energies of all states
        double E_tot_sqrd_1 = 0;      //sum of energies of all states squared
        double E_1 = 0;               //energy after change in position
        double E_prev_1 = 0;          //energy before change in position
        double delta_EL_1;            //change in energy

        double psi_ratio_1;           //ratio of new and old wave function
        double r_1;                   //random number
        int N_rand_1;                 //randomly chosen N
        int dim_rand_1;               //randomly chosen dimension
        double pos_mat_new_1[N][3];   //new position with random position

        int num_or_an_1 = 0;

        //Initialize position matrix for N particles in dim dimentions
        double pos_mat_1[N][3];
        for(int i=0; i<N; i++) {
            for(int j=0; j<dim; j++) {
                pos_mat_1[i][j] = random_position_1(1);
            }
            for(int k=dim; k<3; k++){
                pos_mat_1[i][k] = 0;
            }
        }

        //Initialize wave function
        WaveFunction Psi_1;
        Psi_1.setTrialWF(dim, N, a, num_or_an_1, HO);

        //Add initial energies to averages
        E_1 = Psi_1.E_L(pos_mat_1, alpha[k], beta, omega_HO, omega_z);
        E_tot_1 += E_1;
        E_tot_sqrd_1 += E_1*E_1;

        //Start Monte Carlo iterations
        for(int i=0;i<M;i++){
            //Draw random position, for one particle and one dimention
            N_rand_1 = rand()%N;
            dim_rand_1 = rand()%dim;

            //Set new meatrix equal old one
            memcpy(pos_mat_new_1, pos_mat_1, sizeof(pos_mat_new_1));

            //Proposed new position
            pos_mat_new_1[N_rand_1][dim_rand_1] = pos_mat_1[N_rand_1][dim_rand_1] + (random_position_1(1)-0.5)*steplength;

            //Metropolis algorithm
            psi_ratio_1 = Psi_1.Psi_value_sqrd(pos_mat_new_1, alpha[k], beta)/(Psi_1.Psi_value_sqrd(pos_mat_1, alpha[k], beta));
            E_prev_1 = Psi_1.E_L(pos_mat_1, alpha[k], beta, omega_HO, omega_z);

            //cout << Psi.E_L_num(pos_mat, alpha[k], beta, 0.001) << endl;

            r_1 = ((double)rand() / (double)RAND_MAX);

            if(psi_ratio_1 >= r_1){
                //accept and update pos_mat
                memcpy(pos_mat_1, pos_mat_new_1, sizeof(pos_mat_1)); //maybe more time efficient to only update the one changed position?
            }

            E_1 = Psi_1.E_L(pos_mat_1, alpha[k], beta, omega_HO, omega_z);
            delta_EL_1 = E_1 - E_prev_1;
            E_tot_1 += E_1;
            E_tot_sqrd_1 += E_1*E_1;
        }

        //Calculate <E_l> and <E_L**2>
        double E_L_avg_1 = E_tot_1/M;
        double E_L_avg_sqrd_1 = E_tot_sqrd_1/M;

        cout << "E_L_avg: " << E_L_avg_1 << endl;
        cout << "E_L_avg_tot: " << E_L_avg_sqrd_1 << "\n" << endl;
    }
}
