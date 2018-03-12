#include <iostream>
#include <cmath>
#include <wavefunction.h>
#include <vector>
#include <random>
#include <tools.h>

using namespace std;

//Mersenne Twister RNG
random_device rdd;                   //Will be used to obtain a seed for the random number engine
mt19937 seed(rdd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> diss(0, 1);

double Random_position(){
    return diss(seed);
}

void GradientDecent(int N, int dim, int M, double a, double steplength, bool HO, \
                    double beta, double h, int num_or_an, int BF_H, double timestep)
{
    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);       //Random gaussian
    uniform_int_distribution<> nrand(0, N-1);         //Random number between 0 and N
    uniform_int_distribution<> dimrand(0, dim-1);     //Random number between 0 and dim

    //Parameters
    double alpha        = 0.4;          //Initial guess
    double eps          = 0.001;        //Tolerance
    double eta0         = 0.01;         //Learning rate
    double D            = 0.5;          //Diffusion coeff, to be used in Hastings met.algo
    int T               = 1000;          //Number of iterations (alphas)

    double alpha_old[5];
    for(int i=0;i<5;i++){
        alpha_old[i] = 0;
    }

    for(int iter=0; iter<T; iter++){

        //averages and energies
        double E_tot       = 0;          //sum of energies of all states
        double E_tot_sqrd  = 0;          //sum of energies of all states squared
        double E           = 0;          //energy after change in position
        double psi_E_tot   = 0;          //used for calculating derv of local energy
        double psi_tot     = 0;          //used for calculating derv of local energy

        double psi_ratio;                //ratio of new and old wave function
        int    N_rand;                   //randomly chosen N
        int    dim_rand;                 //randomly chosen dimension

        //Initialize position matrix for N particles in dim dimentions
        vector<vector<double>> pos_mat;
        vector<vector<double>> pos_mat_new;
        pos_mat.resize(N);
        pos_mat_new.resize(N);
        for(int i = 0; i < N; i++){
            pos_mat[i].resize(dim);
            pos_mat_new[i].resize(dim);
        }

        for(auto& i : pos_mat)
            i.resize(dim);
        for(auto& i : pos_mat_new)
            i.resize(dim);

        for(auto& particle : pos_mat)
            for(auto& coord : particle)
                coord = Random_position();

        //Initialize wave function
        WaveFunction Psi;
        Psi.setTrialWF(dim, N, a, num_or_an, HO);

        //Add initial energies to averages
        if(num_or_an == 0) {
            E = Psi.E_L_ana(pos_mat, alpha, beta);
        }
        else if(num_or_an == 1) {
            E = Psi.E_L_num(pos_mat, alpha, beta, h);

        }
        else {
            cout << "num_or_an is out of range" << endl;
        }

        E_tot += E;
        E_tot_sqrd += E*E;
        psi_tot += Psi.Psi_der(pos_mat, beta);
        psi_E_tot += Psi.Psi_der(pos_mat, beta)*E;

        double accept = 0;

        //Start Monte Carlo iterations
        for(int i=0; i<M; i++){
            //Draw random position, for one particle and one dimention
            N_rand   = nrand(seed);
            dim_rand = dimrand(seed);

            //Set new meatrix equal old one
            pos_mat_new = pos_mat;

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
                pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + D*QForce(pos_mat[N_rand][dim_rand], \
                                                alpha, beta, dim_rand)*timestep + eps_gauss(seed)*sqrt(timestep);
                psi_ratio = GreenFuncSum(pos_mat, pos_mat_new, D, timestep, N, alpha, beta, dim)*Psi.Psi_value_sqrd(pos_mat_new, \
                                                alpha, beta)/(Psi.Psi_value_sqrd(pos_mat, alpha, beta));
            }

            if(psi_ratio >= Random_position()){
                //accept and update pos_mat
                pos_mat = pos_mat_new;
                accept += 1;

                if(num_or_an == 0) {
                    E = Psi.E_L_ana(pos_mat, alpha, beta);
                }
                else if(num_or_an == 1) {
                    E = Psi.E_L_num(pos_mat, alpha, beta, h);
                }
            }

            E_tot += E;
            E_tot_sqrd += E*E;
            psi_tot += Psi.Psi_der(pos_mat, beta);
            psi_E_tot += Psi.Psi_der(pos_mat, beta)*E;
        }

        //Calculate <E_L> and <E_L**2>
        double E_L_avg = E_tot/M;
        double E_L_avg_sqrd = E_tot_sqrd/M;
        double accept_ratio = accept/M; //too small for wrong alpha; can say 0, but is not exactly 0
        double variance = E_L_avg_sqrd - E_L_avg*E_L_avg;

        double psi_E_avg = psi_E_tot/M;
        double psi_avg = psi_tot/M;
        double E_L_der = 2*(psi_E_avg - psi_avg*E_L_avg);

        cout << "--- ALPHA: " << alpha << " ---" << endl;
        cout << "E_L_der: " << E_L_der << endl;
        cout << "E_L_avg: " << E_L_avg << endl;
        cout << "Acceptance ratio: " << accept_ratio << endl;
        cout << "Variance: " << variance << "\n" << endl;

        for(int j=0; j<4; j++) {
            alpha_old[j] = alpha_old[j+1];
        }
        alpha_old[4] = alpha;



        if(abs(E_L_der)<eps||abs(alpha-alpha_old[1])/alpha_old[1] < eps) {
            cout <<"FINAL VALUES" << endl;
            if(abs(E_L_der)<eps){
                cout << "small E_L" << endl;
            }
            cout << "Alpha: " << alpha << endl;
            cout << "E_L_avg: " << E_L_avg << endl;
            cout << "E_L_der: " << E_L_der << endl;
            cout << "iteration alpha nr: " << iter << endl;
            break;
        }

        //Update alpha
        alpha = alpha - eta0 * E_L_der;//sqrt(iter + 1.0);

    }
}
