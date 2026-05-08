#pragma once

#include <map>
#include "mesh_data.hpp"

class VTK_Writer{

    private:
        std::string filename;
        std::vector<Node> nodes;
        std::vector<Element> triangles;
        std::map<std::string, std::vector<double>> scalar_fields;
    
    public:
        VTK_Writer(const std::string &filename);

        void set_mesh(
            const std::vector<Node>& p_nodes,
            const std::vector<Element>& p_triangles
        );

        void add_scalar_field(
            const std::string& field_name,
            const std::vector<double>& values
        );

        void clear_fields();
        void set_filename(const std::string& name);
        void write() const;
};