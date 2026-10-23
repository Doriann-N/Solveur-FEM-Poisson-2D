#include "solver.hpp"
#include "linalgebra.hpp"
#include <iostream>

/* Pivot de gauss */
std::vector<double> solve_with_gauss(std::vector<std::vector<double>> K, std::vector<double> F){
    int n_ligne = K.size();
    int n_colonne = K[0].size();

    std::vector<double> U(n_colonne, 0.);

    /* Parcourt les colonnes sauf la dernière */
    for (int k = 0; k < n_colonne-1; k++){
        /* Parcourt les lignes sous la diagonale */
        for (int i = k + 1; i < n_ligne; i++){
            double f = K[i][k] / K[k][k];
            
            /* On modifie toute la ligne i (ligne_i reçoit ligne_i - f * ligne_k)
            (donc on parcourt les colonnes) */
            for (int j = 0; j < n_colonne; j++){
                K[i][j] -= f * K[k][j];
            }
            F[i] -= f * F[k];
        }
    }

    U[n_colonne - 1] = F[n_colonne-1] / K[n_colonne-1][n_colonne-1];
    for (int i = n_colonne - 2; i >= 0; i--){
        /* Calcul de S */
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


/* Gradient conjugué */
std::vector<double> solve_with_gc(const CSRMatrix &A, const std::vector<double> &b, int max_iter, const bool verbose = false){
    /* Premier point choisie est zéro */
    std::vector<double> x(b.size(), 0.);
    /* Le résidu est donc de b */
    std::vector<double> r = b;
    std::vector<double> p = r;
    const double b_squared = dot(b, b);
    const double eps = 1.e-8;

    /* Si le second membre est nul, alors la solution est le vecteur nul */
    if (b_squared == 0.) return x;

    for (int k = 1; k <= max_iter; k++){
        auto Ap = A.mult_with_v(p);
        double r_squared = dot(r, r);
        double p_dot_Ap = dot(p, Ap);
        /* Si (p,Ap) est trop petit */
        if (std::abs(p_dot_Ap) < 1.e-16){
            std::cout << "Iteration " << k << " : Problem : (p,Ap) = 0\n";
            break;
        }
        double alpha = r_squared / p_dot_Ap;
        /* x = x + alpha * p */
        axpy(alpha, p, x);
        /* r = r - alpha * Ap, modification de r_k pour r_next */
        axpy(-alpha, Ap, r);
        double r_next_squared = dot(r, r);

        auto Ax = A.mult_with_v(x);
        std::vector<double> r_exact = b;
        axpy(-1.0, Ax, r_exact);
        double residu = std::sqrt(dot(r_exact, r_exact));
        double residu_relatif = residu/std::sqrt(b_squared);
        
        if (verbose){
            if (k == max_iter) std::cout << "(END)";
            std::cout << "Iteration " << k << " : ||Ax-b|| = " 
                      << residu << " : ||Ax-b||/||b|| = " << residu_relatif << '\n';
        }
        
        if (r_next_squared < eps * eps * b_squared) {
            if (verbose) std::cout << "Iteration " << k 
                                   << " Sorti de boucle : r² / ||b||² < epsilon² (epsilon = " 
                                   << eps << ")\n";
            break;
        }

        double beta = r_next_squared / r_squared;
        // p = r + beta * p;
        std::vector<double> p_new = r;
        axpy(beta, p, p_new);
        p = p_new;
    }
    return x;
}