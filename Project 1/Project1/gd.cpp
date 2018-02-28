#include <gd.h>
#include <iostream>
#include <cmath>
#include <wavefunction.h>
#include <string.h>
#include <vector>
#include <ctime>
#include <random>
#include <fstream>

using namespace std;

//Mersenne Twister RNG
random_device rdd;                   //Will be used to obtain a seed for the random number engine
mt19937 seed(rdd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> diss(0, 1);

double Random_position(){
    return diss(seed);
}

double qForce(double pos, double alpha, double beta, int dim_rand){
    double QF = -4*alpha*pos;
    if (dim_rand==2){
        QF = QF*beta;
    }
    return QF;
}

double Greenfuncsum(double pos_mat[][3], double pos_mat_new[][3], double D, double timestep, int N, double alpha, double beta){
    double GreenSum = 0; //Sum of ratios for all particles

    for(int i=0; i<N; i++) {
        double GreenOld = 0;
        double GreenNew = 0;
        for(int j=0; j<3; j++) {
            GreenOld += (pos_mat_new[i][j] - pos_mat[i][j] - D*timestep*qForce(pos_mat[i][j], alpha, beta, j))*(pos_mat_new[i][j] - pos_mat[i][j] - D*timestep*qForce(pos_mat[i][j], alpha, beta, j));
            GreenNew += (pos_mat[i][j] - pos_mat_new[i][j] - D*timestep*qForce(pos_mat_new[i][j], alpha, beta, j))*(pos_mat[i][j] - pos_mat_new[i][j] - D*timestep*qForce(pos_mat_new[i][j], alpha, beta, j));
        }
        //GreenOld = exp(-GreenOld/(4*D*timestep));
        //GreenNew = exp(-GreenNew/(4*D*timestep));

        GreenSum += exp(GreenOld/GreenNew);
        //cout << GreenOld << " " << GreenNew << " " << GreenSum << endl;
    }
    return GreenSum;
}

void Metropolis(int N, int dim, int M, double a, double steplength, double omega_HO, double omega_z, bool HO, double beta, double h, int num_or_an, int BF_H, double timestep)
{
    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);
    uniform_int_distribution<> nrand(0, N-1);         //Random number between 0 and N
    uniform_int_distribution<> dimrand(0, dim-1);     //Random number between 0 and dim

    double alpha = 0.7;          //Initial guess
    double alpha_old;
    double eps = 0.001;
    double eta0 = 0.001;              //Learning rate
    double D = 0.1;               //Diffusion coeff, to be used in Hastings met.algo
    int T = 50;                     //Number of iterations (alphas)


    //Open file for writing
    //ofstream myfile;
    //myfile.open ("data/local_energy.dat");

    for(int iter=0; iter<T; iter++){

        //averages and energies
        double E_tot       = 0;          //sum of energies of all states
        double E_tot_sqrd  = 0;          //sum of energies of all states squared
        double E           = 0;          //energy after change in position
        double psi_E_tot   = 0;         //used for calculating derv of local energy
        double psi_tot     = 0;         //used for calculating derv of local energy

        double psi_ratio;               //ratio of new and old wave function
        int    N_rand;                  //randomly chosen N
        int    dim_rand;                //randomly chosen dimension

        double pos_mat[N][3];           //current position
        double pos_mat_new[N][3];       //new position with random step

        //Initialize position matrix for N particles in dim dimentions
        for(int i=0; i<N; i++) {
            for(int j=0; j<dim; j++) {
                pos_mat[i][j] = Random_position();
            }
            for(int k=dim; k<3; k++){
                pos_mat[i][k] = 0;
            }
        }


        //Initialize wave function
        WaveFunction Psi;
        Psi.setTrialWF(dim, N, a, num_or_an, HO);

        //Add initial energies to averages
        if(num_or_an == 0) {
            E = Psi.E_L_ana(pos_mat, alpha, beta, omega_HO, omega_z);
        }
        else if(num_or_an == 1) {
            E = Psi.E_L_num(pos_mat, alpha, beta, omega_HO, omega_z, h);

        }
        else {
            cout << "num_or_an is out of range" << endl;
        }

        E_tot += E;
        E_tot_sqrd += E*E;
        psi_tot += Psi.Psi_der(pos_mat, alpha, beta);
        psi_E_tot += Psi.Psi_der(pos_mat, alpha, beta)*E;

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
                pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + (2*Random_position()-1.0)*steplength;
                psi_ratio = Psi.Psi_value_sqrd(pos_mat_new, alpha, beta)/(Psi.Psi_value_sqrd(pos_mat, alpha, beta));
            }
            else if(BF_H == 1){
                //Hastings metropolis

                //Proposed new position
                pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + D*qForce(pos_mat[N_rand][dim_rand], alpha, beta, dim_rand)*timestep + eps_gauss(seed)*sqrt(timestep);
                psi_ratio = Greenfuncsum(pos_mat, pos_mat_new, D, timestep, N, alpha, beta)*Psi.Psi_value_sqrd(pos_mat_new, alpha, beta)/(Psi.Psi_value_sqrd(pos_mat, alpha, beta));
            }

            if(psi_ratio >= Random_position()){
                //accept and update pos_mat
                memcpy(pos_mat, pos_mat_new, sizeof(pos_mat)); //maybe more time efficient to only update the one changed position?
                accept += 1;

                if(num_or_an == 0) {
                    E = Psi.E_L_ana(pos_mat, alpha, beta, omega_HO, omega_z);
                }
                else if(num_or_an == 1) {
                    E = Psi.E_L_num(pos_mat, alpha, beta, omega_HO, omega_z, h);
                }
            }

            E_tot += E;
            E_tot_sqrd += E*E;
            psi_tot += Psi.Psi_der(pos_mat, alpha, beta);
            psi_E_tot += Psi.Psi_der(pos_mat, alpha, beta)*E;
        }
        clock_t end_time = clock();

        //Calculate <E_L> and <E_L**2>
        double E_L_avg = E_tot/M;
        double E_L_avg_sqrd = E_tot_sqrd/M;
        double accept_ratio = accept/M; //too small for wrong alpha; can say 0, but is not exactly 0
        double CPU_time = 1.0*(end_time - start_time)/CLOCKS_PER_SEC;
        double variance = E_L_avg_sqrd - E_L_avg*E_L_avg;

        double psi_E_avg = psi_E_tot/M;
        double psi_avg = psi_tot/M;
        double E_L_der = 2*(psi_E_avg - psi_avg*E_L_avg);

        cout << "E_L_der: " << E_L_der << endl;
        cout << "--- ALPHA: " << alpha << " ---" << endl;
        cout << "E_L_avg: " << E_L_avg << endl;
        cout << "E_L_avg_tot: " << E_L_avg_sqrd << endl;
        cout << "Acceptance ratio: " << accept_ratio << endl;
        cout << "Variance: " << variance << endl;
        cout << "CPU time: " << CPU_time << "\n" << endl;

        alpha_old = alpha;
        //Update alpha
        alpha = alpha - eta0 * E_L_der*sqrt(iter + 1);


        if(abs(E_L_der)<eps&&abs(alpha-alpha_old)<eps){
            cout <<"FINAL VALUES" << endl;
            cout << "alpha: " << alpha << endl;
            cout << "iteration alpha nr: " << iter << endl;
            break;
        }

        //Write to file
        //myfile << alpha << " " << E_L_avg << " " << variance << "\n";
    }

    //Close myfile
    //myfile.close();
}
