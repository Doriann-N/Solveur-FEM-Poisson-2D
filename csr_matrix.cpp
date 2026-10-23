#include "csr_matrix.hpp"

void CSRMatrix::print(){
    int size_row = row.size();
    int size_col = col.size();
    int size_values = values.size();
    printf("\nRow : \n");
    for (int i = 0; i < size_row; i++)
    {
        std::cout << row[i];
        if (i < size_row-1) std::cout << ", ";
        else std::cout << "\n";
    }
    printf("\nCol : \n");
    for (int i = 0; i < size_col; i++){
        std::cout << col[i];
        if (i < size_col-1) std::cout << ", ";
        else std::cout << "\n";
    }
    printf("\nValues : \n");
    for (int i = 0; i < size_values; i++){
        std::cout << values[i];
        if (i < size_values-1) std::cout << ", ";
        else std::cout << "\n";
    }
}

void CSRMatrix::print_row(int i) const{
    std::cout << "Ligne " << i << " :\n";
    for (int j = row[i]; j < row[i + 1]; ++j) {
        std::cout << "  A[" << i << "][" << col[j]
                  << "] = " << values[j] << '\n';
    }
}

CSRMatrix::CSRMatrix(const std::vector<std::vector<double>> &all_values){
    n_ligne = all_values.size();
    n_colonne = all_values[0].size();
    row.reserve(n_ligne + 1);
    
    for (int i = 0; i < n_ligne; i++){
        row.push_back(values.size());
        for(int j = 0; j < n_colonne; j++){
            const double value = all_values[i][j];
            // Non nul
            if (std::abs(value) > epsilon){
                values.push_back(value);
                col.push_back(j);
            }
        }
    }
    row.push_back(values.size());
}

std::vector<double> CSRMatrix::mult_with_v(const std::vector<double> &v) const{
    assert(v.size() == n_colonne);
    std::vector<double> res(n_ligne, 0.);

    #pragma omp parallel for if(n_ligne > 100)
    for (int i = 0; i < n_ligne; i++){
        int begin_index = row[i];
        int next_index = row[i+1];
        double val = 0.;
        for (int j = begin_index; j < next_index; j++){
            val += values[j] * v[col[j]];
        }
        res[i] = val;
    }
    return res;
}

/* void CSRMatrix::create_row_col_value(const std::vector<Triplet>& fusions){
    const double epsilon = 1e-15;
    int current_line = fusions[0].i;
    int size_fusions = fusions.size();
    row.push_back(0);
    for(const auto& t : fusions){
        // Ligne différente
        if (current_line != t.i){
            for (int k = 0; k < t.i - current_line; k++){
                row.push_back(values.size());
            }
            current_line = t.i;
        }
        if (std::abs(t.value) > epsilon){
            values.push_back(t.value);
            col.push_back(t.j);
        }
    }
    row.push_back(values.size());
} */

// Triplets contient les coordonnées et les valeurs de la matrice 
// (de façon trier par ligne puis colonne et de façon unique)
CSRMatrix::CSRMatrix(const int& p_lignes, const int& p_colonnes, std::vector<Triplet>& fusions){
    n_ligne = p_lignes;
    n_colonne = p_colonnes;

    // Créer row, col, value
    int current_line = fusions[0].i;
    row.push_back(0);
    for(const auto& t : fusions){
        // Ligne différente
        if (current_line != t.i){
            for (int k = 0; k < t.i - current_line; k++){
                row.push_back(values.size());
            }
            current_line = t.i;
        }
        if (std::abs(t.value) > epsilon){
            values.push_back(t.value);
            col.push_back(t.j);
        }
    }
    
    // Si jamais y a des lignes vide à la fin, on remplit
    while (row.size() < n_ligne){
        row.push_back(values.size());
    }
    // Puis on ajoute le nombre total de coef non nul
    row.push_back(values.size());
}

/* CSRMatrix::CSRMatrix(const int &lignes, const int &colonnes){
    n_ligne = lignes;
    n_colonne = colonnes;
    row.resize(lignes + 1, 0.);
}
void CSRMatrix::add_element_at(const double &element, const int &i, const int &j){
    if (have_element_at(i, j)){
            
    }
}
bool CSRMatrix::have_element_at(const int& i, const int& j){
    bool have_element_on_line_i = count_zero_on_line(i);
    if (!have_element_on_line_i) return false;
    for (int start = row[i]; start < row[i+1]; start++)
    {
        if (col[start] == j) return true;
    }
    return false;
}
int CSRMatrix::count_zero_on_line(const int &i){
    return row[i+1] - row[i];
}
int CSRMatrix::count_all_zero(){
    return row[n_ligne];
} */