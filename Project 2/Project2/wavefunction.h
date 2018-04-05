#pragma once
#include <vector>

class WaveFunction
{
private:
public:
    WaveFunction() {}
    int setTrialWF              (int X);
    double Psi_value();
    //double E_L(std::vector<std::vector<double>> &pos_mat, double alpha, double beta);
};
