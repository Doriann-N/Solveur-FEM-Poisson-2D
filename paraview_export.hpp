#include "mesh_data.hpp"

std::vector<double> sample_function_on_nodes(
    const std::vector<Node> &nodes, 
    const std::function<double(double, double)> &func
);
void export_for_paraview(
    const std::vector<Node> &nodes, 
    const std::vector<Element> &triangles, 
    std::vector<double> &U, 
    const std::function<double(double, double)> &sol_exacte, 
    const std::string &file_name
);