#pragma once 

#include "mesh_data.hpp"
#include "csr_matrix.hpp"

std::array<std::array<double, 3>, 3> get_elementaire_mat(
    const std::vector<Node> &nodes, 
    const Element &triangle
);
std::vector<std::vector<double>> get_global_mat(
    const std::vector<Node> &nodes, 
    const std::vector<Element> &triangles
);
std::array<double, 3> get_elementaire_f(
    const std::vector<Node> &nodes, 
    const Element &triangle, 
    const std::function<double(double, double)> &f
);
std::vector<double> get_global_f(
    const std::vector<Node> &nodes, 
    const std::vector<Element> &triangles, 
    const std::function<double(double, double)> &f
);
void set_dirichlet_condition(
    const std::vector<Node> &nodes, 
    std::vector<std::vector<double>> &K, 
    std::vector<double> &F
);
void set_dirichlet_condition(
    const std::vector<Node>& nodes, 
    std::vector<Triplet>& fusions,
    std::vector<double>& F
);
std::vector<Triplet> get_global_mat_in_list_form(
    const std::vector<Node> &nodes, 
    const std::vector<Element> &triangles
);
