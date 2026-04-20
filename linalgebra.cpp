#include "linalgebra.hpp"
#include <cassert>

double dot(const std::vector<double>& u, const std::vector<double>& v){
    assert (u.size() == v.size());
    double res = 0;
    for (int i = 0; i < u.size(); i++){
        res += u[i] * v[i];
    }
    return res;
}
// axpy = alpha * x + y, on modifie y
void axpy(double alpha, const std::vector<double>& x, std::vector<double> &y){
    assert (x.size() == y.size());

    for (int i = 0; i< x.size(); i++){
        y[i] += alpha * x[i];
    }
}

// ax = alpha * x, on modifie x
void ax(double alpha, std::vector<double>& x){
    for (int i = 0; i < x.size(); i++){
        x[i] *= alpha;
    }
}