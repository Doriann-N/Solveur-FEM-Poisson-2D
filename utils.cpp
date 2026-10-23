#include "utils.hpp"


/* Ensemble de fonction qui m'aide pour l'affichage et autre */

void display_mat(const std::vector<std::vector<double>> &K){
    int n_i = K.size();
    int n_j = K[0].size();

    for (int i = 0; i < n_i; i++){
        for (int j = 0; j <n_j; j++){
            std::cout << K[i][j] << "  ";
        }
        std::cout << "\n";
    }
}

void display_vec(const std::vector<double> &v){
    int n = v.size();
    for (int i = 0; i < n; i++){
        std::cout << v[i];
        if (i < n-1) std::cout << ", ";
        else std::cout << "\n";
    }
}

void display_vector_int(const std::vector<int>& v){
    int size = v.size();
    for (int i = 0; i < size; i++){
        std::cout << v[i];
        if (i < size-1) std::cout << ", ";
        else std::cout << "\n";
    }
}

int ask_for_discretisation(){
    int N;
    std::cout << "Discrétisation N ? ";
    std::cin >> N;
    return N;
}