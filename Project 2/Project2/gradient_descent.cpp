#include <iostream>
#include <wavefunction.h>
#include <random>
#include <eigen3/Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace std;

//Mersenne Twister RNG
random_device rd;                   //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0, 1);

double random_position(){
    return dis(gen);
}

void GradientDescent(int P, int D, int N, int MC, double sigma, double omega, double steplength) {

    //Constants
    double psi_ratio;               //ratio of new and old wave function
    double eta = 0.01;
    int M = P*D;
    int M_rand;

    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);       //Gaussian distr random number generator
    uniform_int_distribution<> mrand(0, M-1);         //Random number between 0 and N

    for(int iter=0; iter<10; iter++) {
        //averages and energies
        double EL_tot      = 0;          //sum of energies of all states
        double EL_tot_sqrd = 0;          //sum of energies of all states squared
        double EL          = 0;          //energy after change in position
        double EL_sqrd     = 0;

        MatrixXd W = MatrixXd::Random(M,N);     //Weights
        VectorXd a = VectorXd::Random(M);       //Visible biases
        VectorXd b = VectorXd::Random(N);       //Hidden biases
        VectorXd X = VectorXd::Random(M);       //Visible nodes (position coordinates)
        VectorXd X_new;

        VectorXd da_tot           = VectorXd::Zero(M);
        VectorXd daE_tot          = VectorXd::Zero(M);
        VectorXd db_tot           = VectorXd::Zero(N);
        VectorXd dbE_tot          = VectorXd::Zero(N);
        MatrixXd dW_tot           = MatrixXd::Zero(M,N);
        MatrixXd dWE_tot          = MatrixXd::Zero(M,N);

        WaveFunction Psi;
        Psi.setTrialWF(N, M, sigma, omega);

        double E = Psi.EL_calc(X, a, b, W);
        EL += E;
        EL_sqrd += E*E;

        double accept = 0;
        for(int i=0; i<MC; i++) {
            X_new = X;              //Setting new matrix equal to old one

            M_rand = mrand(gen);    //Random particle and dimension

            int sampling = 0;
            if(sampling == 0) {
                X_new(M_rand) = X(M_rand) + (2*random_position() - 1.0)*steplength;
                psi_ratio = Psi.Psi_value_sqrd(a, b, X_new, W)/Psi.Psi_value_sqrd(a, b, X, W);
            }
            else if(sampling == 1) {
                //Hastings
            }
            double r = random_position();

            //cout << psi_ratio << " " << r << endl;

            if(psi_ratio >= r) {
                //accept and update
                X = X_new;
                accept += 1;

                E = Psi.EL_calc(X, a, b, W);
            }
            //cout << "E after " << E << endl;

            VectorXd da;
            VectorXd db;
            MatrixXd dW;
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
        double EL_avg = EL_tot/(MC + 1);
        double EL_avg_sqrd = EL_tot_sqrd/(MC + 1);

        //Gradient descent
        a -= 2*eta*(daE_tot - EL_avg*da_tot)/MC;
        b -= 2*eta*(dbE_tot - EL_avg*db_tot)/MC;
        W -= 2*eta*(dWE_tot - EL_avg*dW_tot)/MC;
        //cout << "a after " << a << '\n' << endl;

        cout << "\n--- Iteration " << iter << " ---" << endl;
        cout << "E_L_avg: " << EL_avg << endl;
        cout << "Acceptance ratio: " << accept/MC << endl;

    }
}
