#pragma once
#include "mesh_data.hpp"

double error_l2(const std::vector<Node> &nodes, 
    const std::vector<Element> &triangles,
    const std::vector<double> &U_h, 
    const std::function<double(double, double)> &sol_exa
);

std::vector<double> errors(
    const std::vector<double>& u_approximation, 
    const std::vector<double>& u_solution
);