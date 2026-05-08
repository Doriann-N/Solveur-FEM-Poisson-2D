#include "vtk_writer.hpp"

#include <cassert>

VTK_Writer::VTK_Writer(const std::string& p_filename){
    filename = p_filename;
}

void VTK_Writer::set_mesh(
    const std::vector<Node>& p_nodes,
    const std::vector<Element>& p_triangles){
        nodes = p_nodes;
        triangles = p_triangles;
}

void VTK_Writer::add_scalar_field(
    const std::string& field_name,
    const std::vector<double>& values){
        scalar_fields[field_name] = values;
}

void VTK_Writer::clear_fields(){
    scalar_fields.clear();
}
void VTK_Writer::set_filename(const std::string& name){
    filename = name;
}
void VTK_Writer::write() const
{

    std::ofstream file(filename + ".vtk");
    if (!file.is_open()){
        std::cout << "Impossible d'ouvrir le fichier\n";
    }

    //std::ofstream file("Solutions/"+ filename + ".vtk");
    file << "# vtk DataFile Version 3.0\n";
    file << "FEM mesh\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    file << "POINTS " << nodes.size() << " float\n";
    for (const auto& n : nodes) file << n.x << " " << n.y << " 0\n";
    
    file << "CELLS " << triangles.size() << " " << triangles.size()*4 << "\n";
    for (const auto& e : triangles) file << "3 " << e.n1 << " " << e.n2 << " " << e.n3 << "\n";
    
    file << "CELL_TYPES " << triangles.size() << "\n";
    for (const auto& e : triangles) file << "5\n";

    file << "POINT_DATA " << nodes.size() << "\n";
    for(const auto& s_field : scalar_fields)
    {
        assert(s_field.second.size() == nodes.size());

        file << "SCALARS "<< s_field.first << " double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (const auto &val : s_field.second)
        {
            file << val << "\n";
        }
    }
    file.close();
}