#pragma once
#include <vector>
#include <cmath>
#include <cassert>
#include <omp.h>

// Explication : https://phtournier.pages.math.cnrs.fr/5mm30/csr/
struct CSRMatrix {
    int n_ligne;
    int n_colonne;
    std::vector<int> row;
    std::vector<int> col;
    std::vector<double> values;

    CSRMatrix(const std::vector<std::vector<double>> &all_values)
    {
        const double eps = 1e-14;
        n_ligne = all_values.size();
        n_colonne = all_values[0].size();
        row.reserve(n_ligne + 1);

        for (int i = 0; i < n_ligne; i++){
            row.push_back(values.size());
            for(int j = 0; j < n_colonne; j++){
                const double value = all_values[i][j];
                // Non nul
                if (std::abs(value) > eps){
                    values.push_back(value);
                    col.push_back(j);
                }
            }
        }
        row.push_back(values.size());
    }

    std::vector<double> mult_with_v(const std::vector<double> &v) const{
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
};