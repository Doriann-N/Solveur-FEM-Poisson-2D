#include "error.hpp"

double error_l2(const std::vector<Node> &nodes, const std::vector<Element> &triangles,
    const std::vector<double> &U_h, 
    const std::function<double(double, double)> &sol_exa){
    double error = 0.;
    int n = triangles.size();
    #pragma omp parallel for reduction(+:error) if(n > 1000)
    for (const auto& t : triangles){
        const Node &A = nodes[t.n1];
        const Node &B = nodes[t.n2];
        const Node &C = nodes[t.n3];
        const double x_c = (A.x + B.x + C.x) / 3.;
        const double y_c = (A.y + B.y + C.y) / 3.;
        const double u_h_c = (U_h[t.n1] + U_h[t.n2] + U_h[t.n3]) / 3.;
        const double e = u_h_c - sol_exa(x_c, y_c);
        error += std::abs(triangle_area(A, B, C)) * e * e;
    }

    return std::sqrt(error);
}

std::vector<double> errors(const std::vector<double>& u_approximation, const std::vector<double>& u_solution)
{
    int n = u_approximation.size();
    std::vector<double> errs(n, 0.);
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        errs[i] = std::abs(u_approximation[i] - u_solution[i]);
    }
    return errs;
}