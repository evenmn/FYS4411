#include <iostream>
#include <cmath>

using namespace std;

double psi_T(double alpha, double beta, double x);

double psi_T_squared(double alpha, double beta, double x);

double random_position(double steplength){
    return ((double)rand() / (double)RAND_MAX)*steplength;
}

double E_L_calc(double alpha, double beta, double x, double omega_HO);


int main()
{
    double alpha = 1;
    double beta = 1;
    double omega_HO = 1;
    //loop over several alphas


    //number of MC cycles
    int M = 10;
    //steplength
    double steplength = 1;
    //initialize averages
    double E_tot = 0;
    double E_tot_sqrd= 0;
    double E = 0;
    double E_prev = 0;
    //
    double x_new = 0;


    //initialize position
    double x = random_position(steplength);

    //add initial energies to averages

    E = E_L_calc(alpha, beta, x, omega_HO);

    E_tot += E;
    E_tot_sqrd += E*E;

    for(int i=0;i<M;i++){

        x_new = random_position(steplength);

        double psi_ratio = (psi_T_squared(alpha, beta, x_new))/(psi_T_squared(alpha, beta, x));

        double E_prev = E_L_calc(alpha, beta, x, omega_HO);

        double r = ((double)rand() / (double)RAND_MAX);

        if(psi_ratio >= r){
            //accept
            x = x_new;

        }

        double E = E_L_calc(alpha, beta, x, omega_HO);
        double delta_EL = E - E_prev;

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

double psi_T(double alpha, double beta, double x){
    return exp(-alpha*x*x);
}

double psi_T_squared(double alpha, double beta, double x){
    return exp(-alpha*x*x*2);
}


double E_L_calc(double alpha, double beta, double x, double omega_HO){
    return -2*x*x*alpha*alpha + alpha + 0.5*omega_HO*omega_HO*x*x;
}
