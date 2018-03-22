#include <iostream>
#include <cmath>

void test_EL(double EL, int N, double alpha, double beta) {
    double tol = 0.1;
    double exact = N*(1 + beta/2)*(1/(4*alpha) + alpha);
    //std::cout << EL << " " << exact << std::endl;
    if(fabs(EL - exact) > tol) {
        std::cout << "Error: EL differs too much from the exact energy" << std::endl;
        std::exit(0);
    }
}
