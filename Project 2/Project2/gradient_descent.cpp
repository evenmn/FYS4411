#include <iostream>
#include "wavefunction.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>
#include <fstream>
#include "hastings_tools.h"
#include "gibbs_tools.h"

using namespace Eigen;
using namespace std;

//Mersenne Twister RNG
random_device rd;                   //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0, 1);
uniform_int_distribution<> hrand(0, 1);

double random_position(){
    return dis(gen);
}

void GradientDescent(int P, double Diff, int D, int N, int MC, int iterations, int sampling, double sigma, double omega, double steplength, double timestep, double eta, bool interaction) {

    //Constants
    double psi_ratio = 0;               //ratio of new and old wave function
    int M = P*D;
    int M_rand = 0;
    int N_rand = 0;

    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);       //Gaussian distr random number generator
    normal_distribution<double> eps_gauss_small(0,0.001);       //Gaussian distr random number generator
    uniform_int_distribution<> mrand(0, M-1);         //Random number between 0 and M
    uniform_int_distribution<> nrand(0, N-1);         //Random number between 0 and N

    MatrixXd W = MatrixXd::Random(M, N);
    VectorXd a = VectorXd::Random(M);
    VectorXd b = VectorXd::Random(N);
    VectorXd X = VectorXd::Random(M);
    VectorXd X_new = VectorXd::Zero(M);
    VectorXd h = VectorXd::Zero(N);
    VectorXd Xa = X - a;
    VectorXd v = b + (W.transpose() * X)/(sigma * sigma);

    for(int i=0; i<N; i++) {
        h(i) = hrand(gen);
    }

    WaveFunction Psi;
    Psi.setTrialWF(N, M, sigma, omega);

    //Open file for writing
    ofstream myfile;
    myfile.open("../data/energy.txt");

    ofstream myfile1;
    myfile1.open("../data/local_energies_interaction_hastings.txt");

    for(int iter=0; iter<iterations; iter++) {
        //averages and energies
        double EL_tot      = 0;          //sum of energies of all states
        double EL_tot_sqrd = 0;          //sum of energies of all states squared
        double E = Psi.EL_calc(X, a, b, W, D, interaction);

        VectorXd da_tot           = VectorXd::Zero(M);
        VectorXd daE_tot          = VectorXd::Zero(M);
        VectorXd db_tot           = VectorXd::Zero(N);
        VectorXd dbE_tot          = VectorXd::Zero(N);
        MatrixXd dW_tot           = MatrixXd::Zero(M,N);
        MatrixXd dWE_tot          = MatrixXd::Zero(M,N);

        double accept = 0;
        for(int i=0; i<MC; i++) {
            X_new = X;              //Setting new matrix equal to old one

            M_rand = mrand(gen);    //Random particle and dimension


            if(sampling == 0) {
                //Standard Metropolis
                X_new(M_rand) = X(M_rand) + (2*random_position() - 1.0)*steplength;
                psi_ratio = Psi.Psi_value_sqrd(a, b, X_new, W)/Psi.Psi_value_sqrd(a, b, X, W);
            }
            else if(sampling == 1) {
                //Metropolis-Hastings
                X_new(M_rand) = X(M_rand) + Diff*QForce(X, a, b, W, N, sigma, M_rand)*timestep + eps_gauss(gen)*sqrt(timestep);
                psi_ratio = GreenFuncSum(X, X_new, a, b, W, N, sigma, timestep, D, Diff)*(Psi.Psi_value_sqrd(a, b, X_new, W)/Psi.Psi_value_sqrd(a, b, X, W));
            }
            else if(sampling == 2) {
                //Gibbs' sampling
                X(M_rand) = x_sampling(a, h, W, sigma, M_rand);
                //Continue writing sampling function for h

            }

            if(psi_ratio >= random_position()&&sampling!=2) {
                //accept and update
                X = X_new;
                accept += 1;
                cout << "accept" << "\n" << endl;
                E = Psi.EL_calc(X, a, b, W, D, interaction);
                VectorXd Xa = X - a;
                VectorXd v = b + (W.transpose() * X)/(sigma * sigma);
            }
            //cout << E << endl;
            if(iter==iterations-1) myfile1 << E << endl;


            VectorXd da = VectorXd::Zero(M);
            VectorXd db = VectorXd::Zero(N);
            MatrixXd dW = MatrixXd::Zero(M,N);
            Psi.Gradient_a(X, a, da);
            Psi.Gradient_b(b, X, W, db);
            Psi.Gradient_W(X, b, W, dW);

            da_tot   += da;
            daE_tot  += E*da;
            db_tot   += db;
            dbE_tot  += E*db;
            dW_tot   += dW;
            dWE_tot  += E*dW;

            EL_tot      += E;
            EL_tot_sqrd += E*E;

        }

        //Calculate <EL> and <EL^2>
        double EL_avg = EL_tot/MC;
        double EL_avg_sqrd = EL_tot_sqrd/MC;
        double variance = (EL_avg_sqrd - EL_avg*EL_avg)/MC;

        cout << "\n--- Iteration " << iter << " ---" << endl;
        cout << "E_L_avg: " << EL_avg << endl;
        cout << "Acceptance ratio: " << accept/MC << endl;
        cout << "Variance " << variance << endl;

        //Gradient descent
        a -= 2*eta*(daE_tot - EL_avg*da_tot)/MC;
        b -= 2*eta*(dbE_tot - EL_avg*db_tot)/MC;
        W -= 2*eta*(dWE_tot - EL_avg*dW_tot)/MC;

        //Write to file
        myfile << EL_avg << "\n";

    }
    //Close myfile
    if(myfile.is_open())  myfile.close();
    if(myfile1.is_open()) myfile1.close();
}
