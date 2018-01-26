#include <iostream>
#include <cmath>
#include <vec3.h>
#include <wavefunction.h>

using namespace std;

double random_position(double steplength){
    return ((double)rand() / (double)RAND_MAX)*steplength;
}

int main()
{
    //loop over several alphas
    double alpha = 1;


    double beta = 1;
    double omega_HO = 1;        //frequency
    int M = 10;                 //number of MC cycles
    double steplength = 1;      //steplength when changing position

    //averages and energies
    double E_tot = 0;
    double E_tot_sqrd= 0;
    double E = 0;
    double E_prev = 0;
    double delta_EL;

    //
    double psi_ratio;
    double r;

    WaveFunction Psi;
    Psi.setTrialWF(1,1);

    //initialize start position
    vec3 r1(random_position(steplength), 0, 0);

    //add initial energies to averages
    E = Psi.E_L(r1, alpha, omega_HO, beta);

    E_tot += E;
    E_tot_sqrd += E*E;


    for(int i=0;i<M;i++){

        vec3 r1_new(random_position(steplength), 0, 0);

        psi_ratio = Psi.Psi_value_sqrd(r1_new, alpha, beta)/(Psi.Psi_value_sqrd(r1, alpha, beta));
        E_prev = Psi.E_L(r1, alpha, omega_HO, beta);

        r = ((double)rand() / (double)RAND_MAX);

        if(psi_ratio >= r){
            //accept
            r1 = r1_new;
        }

        E = Psi.E_L(r1, alpha, omega_HO, beta);
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

