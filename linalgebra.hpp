#pragma once
#include <vector>
double dot(
    const std::vector<double>& u, 
    const std::vector<double>& v
);
// axpy = alpha * x + y, on modifie y
void axpy(
    double alpha, 
    const std::vector<double>& x, 
    std::vector<double> &y
);

// ax = alpha * x, on modifie x
void ax(
    double alpha, 
    std::vector<double>& x
);