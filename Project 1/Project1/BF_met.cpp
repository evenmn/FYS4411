#include <BF_met.h>
#include <iostream>
#include <cmath>
#include <wavefunction.h>
#include <string.h>
#include <vector>
#include <ctime>
#include <random>

using namespace std;

double random_position(){
    random_device rd;                   //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd());                  //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}

double QForce(double pos, double alpha, double beta, int dim_rand){
    double QF = -4*alpha*pos;
    if (dim_rand==2){
        QF = QF*beta;
    }
    return QF;
}

double GreenFuncSum(double pos_mat[][3], double pos_mat_new[][3], double D, double timestep, int N, double alpha, double beta){
    double GreenSum = 0; //Sum of ratios for all particles

    for(int i=0; i<N; i++) {
        double GreenOld = 0;
        double GreenNew = 0;
        for(int j=0; j<3; j++) {
            GreenOld += (pos_mat_new[i][j] - pos_mat[i][j] - D*timestep*QForce(pos_mat[i][j], alpha, beta, j))*(pos_mat_new[i][j] - pos_mat[i][j] - D*timestep*QForce(pos_mat[i][j], alpha, beta, j));
            GreenNew += (pos_mat[i][j] - pos_mat_new[i][j] - D*timestep*QForce(pos_mat_new[i][j], alpha, beta, j))*(pos_mat[i][j] - pos_mat_new[i][j] - D*timestep*QForce(pos_mat_new[i][j], alpha, beta, j));
        }
        GreenOld = exp(-GreenOld/(4*D*timestep));
        GreenNew = exp(-GreenNew/(4*D*timestep));

        GreenSum += GreenOld/GreenNew;
        //cout << GreenOld << " " << GreenNew << " " << GreenSum << endl;
    }
    return GreenSum;
}

void Met_algo(int N, int dim, int M, double a, double steplength, double omega_HO, double omega_z, bool HO, double alpha[], int length_alpha_1, double beta, double h, int num_or_an, int BF_H, double timestep)
{
    //Gaussian distr random number generator
    default_random_engine generator;
    normal_distribution<double> eps_gauss(0,1);

    //Marsenne Twister Random Number Generator
    random_device rd;                               //Will be used to obtain a seed for the random number engine
    mt19937 seed(rd());                              //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> nrand(0, N);         //Random number between 0 and N
    uniform_int_distribution<> dimrand(0, dim);     //Random number between 0 and dim

    for(int k=0; k<length_alpha_1; k++){

        double D = 0.5;                 //Diffusion coeff, to be used in Hastings met.algo

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

        int accept = 0;

        clock_t start_time = clock();
        //Start Monte Carlo iterations
        for(int i=0; i<M; i++){
            //Draw random position, for one particle and one dimention
            N_rand   = nrand(seed);
            dim_rand = dimrand(seed);

            //Set new meatrix equal old one
            memcpy(pos_mat_new, pos_mat, sizeof(pos_mat_new));

            //Choose between Brute force and Hastings
            if(BF_H == 0){
                //Brute force metropolis

                //Proposed new position
                pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + (2*random_position()-1.0)*steplength;
                psi_ratio = Psi.Psi_value_sqrd(pos_mat_new, alpha[k], beta)/(Psi.Psi_value_sqrd(pos_mat, alpha[k], beta));
            }
            else if(BF_H == 1){
                //Hastings metropolis

                //Proposed new position
                pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + D*QForce(pos_mat[N_rand][dim_rand], alpha[k], beta, dim_rand)*timestep + eps_gauss(generator)*sqrt(timestep);
                psi_ratio = GreenFuncSum(pos_mat, pos_mat_new, D, timestep, N, alpha[k], beta)*Psi.Psi_value_sqrd(pos_mat_new, alpha[k], beta)/(Psi.Psi_value_sqrd(pos_mat, alpha[k], beta));
            }

            if(psi_ratio >= random_position()){
                //accept and update pos_mat
                memcpy(pos_mat, pos_mat_new, sizeof(pos_mat)); //maybe more time efficient to only update the one changed position?
                accept += 1;
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
        double accept_ratio = accept*1.0/M;
        double CPU_time = 1.0*(end_time - start_time)/CLOCKS_PER_SEC;
        double variance = E_L_avg_sqrd - E_L_avg*E_L_avg;

        cout << "--- ALPHA: " << alpha[k] << " ---" << endl;
        cout << "E_L_avg: " << E_L_avg << endl;
        cout << "E_L_avg_tot: " << E_L_avg_sqrd << endl;
        cout << "Acceptance ratio: " << accept_ratio << endl;
        cout << "Variance: " << variance << endl;
        cout << "CPU time: " << CPU_time << "\n" << endl;
    }
}
