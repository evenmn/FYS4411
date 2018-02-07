#include <BF_met.h>
#include <iostream>
#include <cmath>
#include <wavefunction.h>
#include <string.h>
#include <vector>
#include <ctime>

using namespace std;

double random_position(){
    return (double)rand() / (double)RAND_MAX;
}

void Met_algo(int N, int dim, int M, double a, double steplength, double omega_HO, double omega_z, bool HO, double alpha[], int length_alpha_1, double beta, double h, int num_or_an)
{

    for(int k=0; k<length_alpha_1; k++){
        //averages and energies
        double E_tot      = 0;          //sum of energies of all states
        double E_tot_sqrd = 0;          //sum of energies of all states squared
        double E          = 0;          //energy after change in position

        double psi_ratio;               //ratio of new and old wave function
        int    N_rand;                  //randomly chosen N
        int    dim_rand;                //randomly chosen dimension

        double pos_mat[N][3];           //current position
        double pos_mat_new[N][3];       //new position with random step

        //Initialize position matrix for N particles in dim dimentions
        for(int i=0; i<N; i++) {
            for(int j=0; j<dim; j++) {
                pos_mat[i][j] = random_position();
            }
            for(int k=dim; k<3; k++){
                pos_mat[i][k] = 0;
            }
        }

        /*
        vector<vector<double>> posmat;
        posmat.resize(N);
        for(int i = 0; i < posmat.size(); i++)
            posmat[i].resize(dim);

        for(auto* i : posmat)
            i.resize(dim);
       */

        //Initialize wave function
        WaveFunction Psi;
        Psi.setTrialWF(dim, N, a, num_or_an, HO);

        //Add initial energies to averages
        if(num_or_an == 0) {
            E = Psi.E_L_ana(pos_mat, alpha[k], beta, omega_HO, omega_z);
        }
        else if(num_or_an == 1) {
            E = Psi.E_L_num(pos_mat, alpha[k], beta, omega_HO, omega_z, h);
        }
        else {
            cout << "num_or_an is out of range" << endl;
        }

        E_tot += E;
        E_tot_sqrd += E*E;

        clock_t start_time = clock();
        //Start Monte Carlo iterations
        for(int i=0; i<M; i++){
            //Draw random position, for one particle and one dimention
            N_rand   = rand()%N;
            dim_rand = rand()%dim;

            //Set new meatrix equal old one
            memcpy(pos_mat_new, pos_mat, sizeof(pos_mat_new));

            //Proposed new position
            pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + (random_position()-0.5)*steplength;

            //Metropolis algorithm
            psi_ratio = Psi.Psi_value_sqrd(pos_mat_new, alpha[k], beta)/(Psi.Psi_value_sqrd(pos_mat, alpha[k], beta));

            if(psi_ratio >= random_position()){
                //accept and update pos_mat
                memcpy(pos_mat, pos_mat_new, sizeof(pos_mat)); //maybe more time efficient to only update the one changed position?
            }

            if(num_or_an == 0) {
                E = Psi.E_L_ana(pos_mat, alpha[k], beta, omega_HO, omega_z);
            }
            else if(num_or_an == 1) {
                E = Psi.E_L_num(pos_mat, alpha[k], beta, omega_HO, omega_z, h);
            }

            E_tot += E;
            E_tot_sqrd += E*E;
        }
        clock_t end_time = clock();

        //Calculate <E_l> and <E_L**2>
        double E_L_avg = E_tot/M;
        double E_L_avg_sqrd = E_tot_sqrd/M;
        double CPU_time = 1.0*(end_time - start_time)/CLOCKS_PER_SEC;

        cout << "--- ALPHA: " << alpha[k] << " ---" << endl;
        cout << "E_L_avg: " << E_L_avg << endl;
        cout << "E_L_avg_tot: " << E_L_avg_sqrd << endl;
        cout << "CPU time: " << CPU_time << "\n" << endl;
    }
}