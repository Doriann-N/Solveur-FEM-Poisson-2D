#include "solver.hpp"
#include "linalgebra.hpp"
#include <iostream>
// Pivot de gauss
std::vector<double> solve_with_gauss(std::vector<std::vector<double>> K, std::vector<double> F){
    int n_ligne = K.size();
    int n_colonne = K[0].size();

    std::vector<double> U(n_colonne, 0.);

    // Parcourt les colones sauf la dernière
    for (int k = 0; k < n_colonne-1; k++){
        // Parcourt les lignes sous la diagonale
        for (int i = k + 1; i < n_ligne; i++){
            double f = K[i][k] / K[k][k];
            
            // On modifie toute la ligne i (on Ligne_i recoint - f * Ligne_k)
            // (donc on parcourt les colonnes)
            for (int j = 0; j < n_colonne; j++){
                K[i][j] -= f * K[k][j];
            }
            F[i] -= f * F[k];
        }
    }

    U[n_colonne - 1] = F[n_colonne-1] / K[n_colonne-1][n_colonne-1];
    for (int i = n_colonne - 2; i >= 0; i--){
        // Calcul de S
        double S = 0.;
        for (int j = i + 1; j < n_colonne; j++){
            S += K[i][j] * U[j];
        }

        U[i] = (F[i] - S) / K[i][i];
    }


    //display_mat(K);
    //display_vec(F);
    //display_vec(U);
    return U;
}



// Gradient conjugué
std::vector<double> solve_with_gc(const CSRMatrix &A, const std::vector<double> &b, int max_iter){
    std::vector<double> x(b.size(), 0.);
    std::vector<double> r = b;
    std::vector<double> p = r;
    const double eps = 1.e-14;
    for (int k = 0; k < max_iter; k++){
        auto Ap = A.mult_with_v(p);
        double r_squared = dot(r, r);
        double p_dot_Ap = dot(p, Ap);
        if (std::abs(p_dot_Ap) < eps) 
            break;
        double alpha = r_squared / p_dot_Ap;
        // x = x + alpha * p
        axpy(alpha, p, x);
        // r = r - alpha * Ap
        axpy(-alpha, Ap, r); // modifie r_k pour r_next
        double r_next_squared = dot(r, r);
        if (r_next_squared < eps) {
            std::cout << "Iteration : " << k << " Residu : " << r_next_squared << "\n";
            break;
        }  
        double beta = r_next_squared / r_squared;
        // p = beta * p 
        ax(beta, p);
        // p = p + r
        axpy(1., r, p);
    }

    return x;
}