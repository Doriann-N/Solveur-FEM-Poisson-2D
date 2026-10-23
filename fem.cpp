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

// Renvoit une liste de triplet contenant (i, j, value) de la matrice
// de façon ordonnée par ligne puis colonne
std::vector<Triplet> get_global_mat_in_list_form(const std::vector<Node> &nodes, const std::vector<Element> &triangles){
    const int n_ligne = nodes.size();
    std::vector<Triplet> triplets;
    for (const auto& t : triangles){
        auto Ke = get_elementaire_mat(nodes, t);
        const int indices[3] = {t.n1, t.n2, t.n3};
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                Triplet triplet(indices[i], indices[j], Ke[i][j]);
                triplets.push_back(triplet);
            }
        }
    }

    // Trier par ligne, puis par colonne
    std::sort(triplets.begin(), triplets.end());
    //std::printf("TRILLAGE :\n");
    //for(const auto &t: triplets) t.print();

    // Fusionner les triplets
    const int size_triplets = triplets.size();
    std::vector<Triplet> fusions;
    Triplet current = triplets[0];
    for(int j = 1; j < size_triplets; j++){
        const Triplet& t = triplets[j];
        if (current.i == t.i && current.j == t.j){
            current.value += t.value;
        }
        else{
            fusions.push_back(current);
            current = t;
        }
    }
    fusions.push_back(current);
    return fusions;
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
    //double f_AB = f((A.x + B.x)/2., (A.y + B.y)/2.);
    //double f_AC = f((A.x + C.x)/2., (A.y + C.y)/2.);
    //double f_BC = f((B.x + C.x)/2., (B.y + C.y)/2.);
    //f_val = (f_AB + f_AC + f_BC) / 3.;
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

// fusions est donc déjà trier par ligne puis colonne
void set_dirichlet_condition(const std::vector<Node>& nodes, 
            std::vector<Triplet>& fusions, std::vector<double>& F){
    double x_limit = 1.;
    double y_limit = 1.;
    
    int n = nodes.size();

    std::vector<bool> boundary(n);
    for (int i = 0; i < n; i++){
        boundary[i] = nodes[i].is_on_boundary(x_limit, y_limit);
    }
    // modifier matrice 
    for (auto& t : fusions){
        bool boundary_i = boundary[t.i];
        bool boundary_j = boundary[t.j];
        if (boundary_i || boundary_j){
            if (t.i == t.j && boundary_i) {
                t.value = 1.;
            }
            else {
                t.value = 0.;
            }
        }
    }

    // modifier le second membre
    for(int i=0; i<n; i++)
    {
        if(boundary[i])
            F[i] = 0.;
    }
}