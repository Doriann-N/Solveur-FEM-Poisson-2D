#pragma once
#include <functional>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

struct Node{
    double x, y;

    bool is_on_boundary(double x_limit, double y_limit) const{
        double eps = 1e-10;
        return abs(x - x_limit) < eps or abs(y - y_limit) < eps or 
        abs(x) < eps or abs(y) < eps; 
    }
};
struct Element{
    int n1, n2, n3;
};

// On génère les points sur [0, 1]²...
std::vector<Node> setup_nodes(int N);
int ij_to_index(int i, int j, int N);
std::vector<Element> setup_triangles(int N);
double triangle_area(const Node& A, const Node& B, const Node& C);

Node oui_voila();