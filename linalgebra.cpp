#include "linalgebra.hpp"
#include <cassert>
#include <omp.h>

double dot(const std::vector<double>& u, const std::vector<double>& v){
    int n = u.size();
    assert (n == v.size());
    double res = 0;
    #pragma omp parallel for reduction(+:res) if(n > 1000)
    for (int i = 0; i < n; i++){
        res += u[i] * v[i];
    }
    return res;
}
// axpy = alpha * x + y, on modifie y
void axpy(double alpha, const std::vector<double>& x, std::vector<double> &y){
    int n = x.size();
    assert (n == y.size());

    #pragma omp parallel for if(n > 1000)
    for (int i = 0; i< n; i++){
        y[i] += alpha * x[i];
    }
}

// ax = alpha * x, on modifie x
void ax(double alpha, std::vector<double>& x){
    int n = x.size();
    #pragma omp parallel for if(n > 1000)
    for (int i = 0; i < n; i++){
        x[i] *= alpha;
    }
}