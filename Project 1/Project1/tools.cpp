#include <tools.h>
#include <iostream>
#include <cmath>

double QForce(double pos, double alpha, double beta, int dim_rand){
    double QF = -4*alpha*pos;
    if (dim_rand==2){
        QF = QF*beta;
    }
    return QF;
}

double GreenFuncSum(std::vector<std::vector<double>> &pos_mat, std::vector<std::vector<double>> &pos_mat_new, \
                    double D, double timestep, int N, double alpha, double beta, int dim){
    double GreenSum = 0; //Sum of ratios for all particles

    for(int i=0; i<N; i++) {
        double GreenOld = 0;
        double GreenNew = 0;
        for(int j=0; j<dim; j++) {
            GreenOld += (pos_mat_new[i][j] - pos_mat[i][j] - D*timestep*QForce(pos_mat[i][j], alpha, beta, j))* \
                        (pos_mat_new[i][j] - pos_mat[i][j] - D*timestep*QForce(pos_mat[i][j], alpha, beta, j));
            GreenNew += (pos_mat[i][j] - pos_mat_new[i][j] - D*timestep*QForce(pos_mat_new[i][j], alpha, beta, j))* \
                        (pos_mat[i][j] - pos_mat_new[i][j] - D*timestep*QForce(pos_mat_new[i][j], alpha, beta, j));
        }

        GreenSum += exp(GreenOld/GreenNew);
    }
    return GreenSum;
}
