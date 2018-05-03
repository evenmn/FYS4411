#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

void test_energy_convergence(double energy, double omega, int M, bool interaction) {

    double error  = 0;
    double tolerance_w_int = 0.1;
    double tolerance_wo_int = 0.05;

    if(interaction) {
        if(M==4) {
            if(fabs(energy - omega * 3) > tolerance_w_int) {
                cout << "Warning: Energy error is larger than the tolerance" << endl;
            }
        }
        else{
            //No analytical values
        }
    }
    else {
        if(fabs(energy - omega * 0.5 * M) > tolerance_wo_int) {
            cout << "Warning: Energy error is larger than the tolerance" << endl;
        }
    }
}
