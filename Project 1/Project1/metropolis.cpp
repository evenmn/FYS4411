#include <iostream>
#include <wavefunction.h>
#include <vector>
#include <ctime>
#include <random>
#include <fstream>
#include <tools.h>
#include <test.h>

using namespace std;

//Mersenne Twister RNG
random_device rd;                   //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0, 1);

double random_position(){
    return dis(gen);
}


void Met_algo(int N, int dim, int M, double a, double steplength, double alpha[], \
              int length_alpha_1, double beta, double h, int num_or_an, int BF_H, double timestep, int one_body)
{

    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);       //Gaussian distr random number generator
    uniform_int_distribution<> nrand(0, N-1);         //Random number between 0 and N
    uniform_int_distribution<> dimrand(0, dim-1);     //Random number between 0 and dim

    //Open file for writing (will write for a specific alpha)
    //ofstream energy_file;
    //energy_file.open ("../data/local_energy.dat");
    //ofstream local_energy_var;
    //local_energy_var.open ("../data/energy.txt");

    for(int k=0; k<length_alpha_1; k++){

        double D = 0.5;                 //Diffusion coeff, to be used in Hastings met.algo

        //averages and energies
        double E_tot      = 0;          //sum of energies of all states
        double E_tot_sqrd = 0;          //sum of energies of all states squared
        double E          = 0;          //energy after change in position

        double psi_ratio;               //ratio of new and old wave function
        int    N_rand;                  //randomly chosen N
        int    dim_rand;                //randomly chosen dimension

        //Initialize position matrix for N particles in dim dimentions
        vector<vector<double>> pos_mat;
        vector<vector<double>> pos_mat_new;
        pos_mat.resize(N);
        pos_mat_new.resize(N);
        for(int i = 0; i < pos_mat.size(); i++){
            pos_mat[i].resize(dim);
            pos_mat_new[i].resize(dim);
        }

        for(auto& i : pos_mat)
            i.resize(dim);
        for(auto& i : pos_mat_new)
            i.resize(dim);

        // Initialize the position matrix
        for(int i=0; i<N; i++) {
            LOOP:
            for(int j=0; j<dim; j++) {
                pos_mat[i][j] = random_position();
            }
            // Ensure that the particles are separated by a
            for(int j=0; j<i; j++) {
                double distij = 0;
                for(int d=0; d<dim; d++) {
                    distij += (pos_mat[i][d] - pos_mat[j][d])*(pos_mat[i][d] - pos_mat[j][d]);
                }
                if(sqrt(distij) < a) {
                    goto LOOP;
                }
            }
        }


        //Initialize wave function
        WaveFunction Psi;
        Psi.setTrialWF(dim, N, a, num_or_an);

        //Add initial energies to averages
        if(num_or_an == 0) {
            E = Psi.E_L_ana(pos_mat, alpha[k], beta);
        }
        else if(num_or_an == 1) {
            E = Psi.E_L_num(pos_mat, alpha[k], beta, h);
        }
        else {
            cout << "num_or_an is out of range" << endl;
        }

        E_tot += E;
        E_tot_sqrd += E*E;

        //Define bins for the one body density measure
        int number_of_bins = 500;
        double max_radius = 3;
        double radius_step = max_radius/number_of_bins;
        double bin_array[number_of_bins];
        double bin_dist[number_of_bins];

        ofstream ob_file;

        if(one_body == 1) {
            for(int i=0; i<number_of_bins; i++){
                bin_array[i] = i * radius_step;
                bin_dist[i] = 0;
            }

            //Open file for writing (will write for a specific alpha)
            ob_file.open ("../data/ob_density.dat");
        }

        double accept = 0;

        clock_t start_time = clock();
        //Start Monte Carlo iterations
        for(int i=0; i<M; i++){

            //if(i%(M/100)==0){
            //    cout << "Iteration: " << i << " of " << M << endl;
            //}

            //Draw random position, for one particle and one dimention
            N_rand   = nrand(gen);
            dim_rand = dimrand(gen);

            //Set new meatrix equal old one
            pos_mat_new = pos_mat;

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
                pos_mat_new[N_rand][dim_rand] = pos_mat[N_rand][dim_rand] + D*QForce(pos_mat[N_rand][dim_rand], alpha[k], beta, \
                                                dim_rand)*timestep + eps_gauss(gen)*sqrt(timestep);
                psi_ratio = GreenFuncSum(pos_mat, pos_mat_new, D, timestep, N, alpha[k], beta, dim)*Psi.Psi_value_sqrd(pos_mat_new, \
                            alpha[k], beta)/(Psi.Psi_value_sqrd(pos_mat, alpha[k], beta));
            }

            if(psi_ratio >= random_position()){
                //accept and update pos_mat
                pos_mat = pos_mat_new;
                accept += 1;

                if(num_or_an == 0) {
                    E = Psi.E_L_ana(pos_mat, alpha[k], beta);
                }
                else if(num_or_an == 1) {
                    E = Psi.E_L_num(pos_mat, alpha[k], beta, h);
                }
            }
            if(one_body == 1) {
                for(int l=0; l<N; l++){
                    double r = sqrt(pos_mat[l][0]*pos_mat[l][0] + pos_mat[l][1]*pos_mat[l][1] + pos_mat[l][2]*pos_mat[l][2]);
                    double err = 1000000;
                    int bin_nr = 0;
                    for(int j=0; j<number_of_bins; j++) {
                        double e = fabs(bin_array[j] - r);
                        if(e < err) {
                            err = e;
                            bin_nr = j;
                        }
                    }
                    bin_dist[bin_nr] += 1;
                }
            }

            E_tot += E;
            E_tot_sqrd += E*E;

            //If blocking to analyze error
            //local_energy_var << E << endl;

        }
        clock_t end_time = clock();

        if(one_body == 1){
            //Write to file
            for(int j=0; j<number_of_bins; j++) {
               ob_file << bin_dist[j]/(bin_array[j]*bin_array[j]*M) << "\n";
            }
            //Close myfile
            ob_file.close();
        }

        //local_energy_var.close();


        //Calculate <E_l> and <E_L**2>
        double E_L_avg = E_tot/M;
        double E_L_avg_sqrd = E_tot_sqrd/M;
        double accept_ratio = accept*1.0/M;
        double CPU_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
        double variance = E_L_avg_sqrd - E_L_avg*E_L_avg;

        //Tests
        if(a==0) test_EL(E_L_avg, N, alpha[k], beta);


        cout << "--- ALPHA: " << alpha[k] << " ---" << endl;
        cout << "E_L_avg: " << E_L_avg << endl;
        cout << "Acceptance ratio: " << accept_ratio << endl;
        cout << "Variance: " << variance << endl;
        cout << "CPU time: " << CPU_time << "\n" << endl;


        //Write to file
        //energy_file << alpha << " " << E_L_avg << " " << variance << "\n";
    }

    //Close myfile
    //energy_file.close();
}
