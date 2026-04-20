#pragma once
#include <vector>
#include "csr_matrix.hpp"

// Pivot de gauss
std::vector<double> solve_with_gauss(std::vector<std::vector<double>> K, std::vector<double> F);

// Gradient conjugué
std::vector<double> solve_with_gc(const CSRMatrix &A, const std::vector<double> &b, int max_iter);