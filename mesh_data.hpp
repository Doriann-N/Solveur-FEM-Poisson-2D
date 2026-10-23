#pragma once
#include <functional>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

struct Node{
    double x, y;

    // on pourrait sotcker l'info dans un booléen
    bool is_on_boundary(double x_limit, double y_limit) const{
        double eps = 1e-10;
        //printf("x_lim, y_lim, %f, %f\n", x_limit, y_limit);
        bool on_boundary = std::abs(x - x_limit) < eps || std::abs(y - y_limit) < eps || std::abs(x) < eps || std::abs(y) < eps; 
        //printf("(%f, %f) (on boundary : %d)\n", x, y, on_boundary);
        return on_boundary;
    }
};
struct Element{
    int n1, n2, n3;
};

// On génère les points sur [0, 1]²...
std::vector<Node> setup_nodes(int N);

int ij_to_index(int i, int j, int N);

std::vector<Element> setup_triangles(int N, bool alternate_triangle_mesh);

double triangle_area(const Node& A, const Node& B, const Node& C);

std::vector<double> sample_function_on_nodes(
    const std::vector<Node> &nodes, 
    const std::function<double(double, double)> &func
);

void print_info_nodes(const std::vector<Node> &nodes);