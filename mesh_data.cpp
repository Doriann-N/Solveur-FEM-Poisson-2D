#include "mesh_data.hpp"


// On génère les points sur [0, 1]²...
std::vector<Node> setup_nodes(int N){
    std::vector<Node> nodes;
    nodes.reserve((N+1)*(N+1));
    for (int j = 0; j <= N; j++){
        for (int i = 0; i <= N; i++){
            Node node;
            node.x = (double)i / N;
            node.y = (double)j / N;

            nodes.push_back(node);
        }
    }
    return nodes;
}

int ij_to_index(int i, int j, int N){
    return i + j * (N+1);
}

// Défini l'index des points du triangle, donc en soit le maillage
// maillage actuel : (Inversé une fois sur 2)
// +*****+*****+...
// **    **    *
// * *   * *   *
// *  *  *  *  *
// *   * *   * *
// *    **    **
// +*****+*****+...
std::vector<Element> setup_triangles(int N){
    std::vector<Element> triangles;
    triangles.reserve(2*N*N);
    for (int j = 0; j < N; j++){
        for (int i = 0; i < N; i++){

            // Triangle gauche
            Element t1;
            // Triangle droit
            Element t2;
            if ((i + j) % 2 == 0){
                t1.n1 = ij_to_index(i, j, N);
                t1.n2 = ij_to_index(i+1, j, N);
                t1.n3 = ij_to_index(i, j+1, N);

                t2.n1 = ij_to_index(i+1, j+1, N);
                t2.n2 = ij_to_index(i, j+1, N);
                t2.n3 = ij_to_index(i+1, j, N);
            } else{
                t1.n1 = ij_to_index(i, j, N);
                t1.n2 = ij_to_index(i+1, j, N);
                t1.n3 = ij_to_index(i+1, j+1, N);

                t2.n1 = ij_to_index(i, j, N);
                t2.n2 = ij_to_index(i+1, j+1, N);
                t2.n3 = ij_to_index(i, j+1, N);
            }
            
            triangles.push_back(t1);
            triangles.push_back(t2); 
        }
    }
    return triangles;
}
//#include <iostream>
// Air orienté, moitié du produit vectoriel 
double triangle_area(const Node& A, const Node& B, const Node& C){
    double double_area = ((B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y));
    //if (double_area < 0) std::cout << "Air négative, attention...\n";
    return 0.5 * double_area;
}

Node oui_voila(){
    return Node{1, 1};
}