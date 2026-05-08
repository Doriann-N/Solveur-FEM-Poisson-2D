#include "fem.hpp"
// Matrice élémentaire, pour un triangle
std::array<std::array<double, 3>, 3> get_elementaire_mat(const std::vector<Node> &nodes, 
                                    const Element &triangle){
    std::array<std::array<double, 3>, 3> Ke = {};
    
    const Node& A = nodes[triangle.n1];
    const Node& B = nodes[triangle.n2];
    const Node& C = nodes[triangle.n3];
    double area = triangle_area(A, B, C);

    std::array<std::array<double, 2>, 3> gradients;
    double a = 1./(2 * area);
    gradients[0][0] = (B.y - C.y) * a;
    gradients[0][1] = (C.x - B.x) * a;

    gradients[1][0] = (C.y - A.y) * a;
    gradients[1][1] = (A.x - C.x) * a;

    gradients[2][0] = (A.y - B.y) * a;
    gradients[2][1] = (B.x - A.x) * a;


    // Calcul de Ke(i,j) = air * grad_i . grad_j
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            double dot = gradients[i][0] * gradients[j][0] + gradients[i][1] * gradients[j][1];
            Ke[i][j] = area * dot;
        }
    }
    return Ke;
}

// Get the K in the KU = F system
std::vector<std::vector<double>> 
get_global_mat(const std::vector<Node> &nodes, const std::vector<Element> &triangles){
    std::vector<std::vector<double>> K;
    int n = nodes.size();
    K.resize(n, std::vector<double>(n, 0.));

    for (const auto& t : triangles){
        auto Ke = get_elementaire_mat(nodes, t);
        const int indices[3] = {t.n1, t.n2, t.n3};
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                K[indices[i]][indices[j]] += Ke[i][j]; 
            }
        }
    }
    return K;
}

std::array<double, 3> get_elementaire_f(const std::vector<Node> &nodes, 
                        const Element &triangle, const std::function<double(double, double)> &f){
    std::array<double, 3> F_e = {};
    
    const Node &A = nodes[triangle.n1];
    const Node &B = nodes[triangle.n2];
    const Node &C = nodes[triangle.n3];
    

    // Améliorable avec une meilleure quadrature, à voir plus tard
    // On évalue f au centre du triangle
    double f_val = f((A.x + B.x + C.x) / 3., (A.y + B.y + C.y) / 3.);
    double area = triangle_area(A, B, C);
    F_e[0] = f_val * area / 3.;
    F_e[1] = f_val * area / 3.;
    F_e[2] = f_val * area / 3.;
    return F_e;
}

// Get the F in the KU = F system
std::vector<double> get_global_f(const std::vector<Node> &nodes, 
                const std::vector<Element> &triangles, const std::function<double(double, double)> &f){
    std::vector<double> F(nodes.size(), 0.);
    for (const auto& t : triangles){
        int indices[3] = {t.n1, t.n2, t.n3};
        
        auto F_e = get_elementaire_f(nodes, t, f);
        for (int i = 0; i < 3; i++){
            F[indices[i]] += F_e[i];
        } 
    }
    return F;
}

// Modifie K et F pour prendre en compte les conditions au bord
void set_dirichlet_condition(const std::vector<Node> &nodes, 
            std::vector<std::vector<double>> &K, std::vector<double> &F){
    double x_limit = 1.;
    double y_limit = 1.;
    int n = nodes.size();
    for (int i = 0; i < n; i++){
        //printf("On boundary %d: %d\n", i, nodes[i].is_on_boundary(x_limit, y_limit));
        bool on_boundary = nodes[i].is_on_boundary(x_limit, y_limit);
        if (on_boundary){
            for (int j = 0; j < n; j++){
                K[i][j] = 0.;
                //K[j][i] = 0.;
            }
            K[i][i] = 1.;
            F[i] = 0.;
        }
    }
}