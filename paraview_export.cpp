#include "paraview_export.hpp"


std::vector<double> sample_function_on_nodes(const std::vector<Node> &nodes, 
    const std::function<double(double, double)> &func){
    std::vector<double> sample;
    sample.reserve(nodes.size());
    for (const auto &node : nodes){
        sample.push_back(func(node.x, node.y));
    }
    return sample;
}

void export_for_paraview(const std::vector<Node> &nodes, const std::vector<Element> &triangles, 
    std::vector<double> &U, const std::function<double(double, double)> &sol_exacte, const std::string &file_name){
    std::ofstream file("Solutions/"+ file_name + ".vtk");
    file << "# vtk DataFile Version 3.0\n";
    file << "FEM mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << nodes.size() << " float\n";
    for (auto& n : nodes)
        file << n.x << " " << n.y << " 0\n";
    file << "CELLS " << triangles.size() << " " << triangles.size()*4 << "\n";
    for (auto& e : triangles)    
        file << "3 " << e.n1 << " " << e.n2 << " " << e.n3 << "\n";
    file << "CELL_TYPES " << triangles.size() << "\n";
    for (auto& e : triangles) file << "5\n";


    file << "POINT_DATA " << U.size() << "\n";
    file << "SCALARS U_h double 1\n";
    file << "LOOKUP_TABLE default\n";
    
    for (const auto &val : U) {
        file << val << "\n";
    }

    file << "SCALARS U_exact double 1\n";
    file << "LOOKUP_TABLE default\n";


    std::vector<double> U_exact = sample_function_on_nodes(nodes, sol_exacte);
    
    for (const auto &val : U_exact){
        file << val << "\n";
    }

    file << "SCALARS error double 1\n";
    file << "LOOKUP_TABLE default\n";

    std::vector<double> error(nodes.size(), 0.);
    for (int i = 0; i < nodes.size(); i++){
        error[i] = std::abs(U[i] - U_exact[i]);
    }

    for (const auto &val : error){
        file << val << "\n";
    }

    file.close();
}