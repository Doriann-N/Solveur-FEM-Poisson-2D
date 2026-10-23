#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <omp.h>
#include <functional>
#include "utils.hpp"

/* Ce format de matrice permet de garder l'information sur les zéros non nuls.
Pour décoder ce format, on regarde ligne par ligne. row[0] indique l'index où regarder col, et row[0 + 1]
nous dit le début de la ligne suivante. 

Exemple : 
row = [0, 2, 4, 7, 8] (Dernier terme = nombre de terme non nul) (taille = N + 1 pour une matrice de taille N * N)
col = [0, 1, 1, 2, 0, 1, 2, 1] (Indique la colonne où placer les termes val)
val = [2, 6, 1, 7, -1, 13, 5, 3] (Valeur non nul à placer)

On regarde row[0] = 0 et row[1] = 2, donc de l'index 0 à 1, on sait que la ligne 0 possède 2 et 6
On regarde row[1] = 2 et row[2] = 4, donc de l'index 2 à 3, on sait que la ligne 1 possède 1 et 7
On regarde row[2] = 4 et row[4] = 7, donc de l'index 4 à 6, on sait que la ligne 2 possède -1, 13 et 5
On regarde row[4] = 7 qui indique la dernière ligne, donc de l'index 7 à la fin, donc la dernière ligne possède 3

Après on les place correctement en regarde leur colonne. On obtient :

[[2, 6, 0, 0],
[0, 1, 7, 0],
[-1, 13, 5, 0],
[0, 3, 0, 0]]
*/

// Explication : https://phtournier.pages.math.cnrs.fr/5mm30/csr/
class CSRMatrix {
    private:
        double epsilon = 1e-15;
    public:
        int n_ligne = 0;
        int n_colonne = 0;
        std::vector<int> row;
        std::vector<int> col;
        std::vector<double> values;

        void print();
        void print_row(int i) const;
        //bool have_element_at(const int& i, const int& j);
        //int count_zero_on_line(const int &i);
        //int count_all_zero();
        //void create_row_col_value(const std::vector<Triplet>& fusions);
        CSRMatrix() = default;
        /* Triplets contient les coordonnées et les valeurs ajouté à 
        la matrice, on doit donc trier les triplets puis fusionner 
        les triplets et ensuite fabriquer row, col et values. */
        CSRMatrix(const int& p_lignes, const int& p_colonnes, std::vector<Triplet>& triplets);
        CSRMatrix(const std::vector<std::vector<double>> &all_values);
        std::vector<double> mult_with_v(const std::vector<double> &v) const;
};