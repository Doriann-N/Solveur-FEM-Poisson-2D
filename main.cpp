#include <iostream>
#include <vector>
#include <utility>
#include <array>
#include <cmath>
#include <functional>
#include <fstream>
#include <cassert>

struct Node{
    double x, y;

    bool is_on_boundary(double x_limit, double y_limit) const{
        double eps = 1e-10;
        return std::abs(x - x_limit) < eps or std::abs(y - y_limit) < eps or 
        std::abs(x) < eps or std::abs(y) < eps; 
    }
};
struct Element{
    int n1, n2, n3;
};

// implémentation basé sur : https://phtournier.pages.math.cnrs.fr/5mm30/csr/
struct CSRMatrix {
    int n_ligne;
    int n_colonne;
    std::vector<int> row;
    std::vector<int> col;
    std::vector<double> values;

    CSRMatrix(const std::vector<std::vector<double>> &all_values)
    {
        const double eps = 1e-14;
        n_ligne = all_values.size();
        n_colonne = all_values[0].size();
        row.reserve(n_ligne + 1);

        for (int i = 0; i < n_ligne; i++){
            row.push_back(values.size());
            for(int j = 0; j < n_colonne; j++){
                const double value = all_values[i][j];
                // Non nul
                if (std::abs(value) > eps){
                    values.push_back(value);
                    col.push_back(j);
                }
            }
        }
        row.push_back(values.size());
    }

    std::vector<double> mult_with_v(const std::vector<double> &v) const{
        assert(v.size() == n_colonne);

        std::vector<double> res(n_ligne, 0.);

        for (int i = 0; i < n_ligne; i++){
            
            int begin_index = row[i];
            int next_index = row[i+1];
            double val = 0.;
            for (int j = begin_index; j < next_index; j++){
                val += values[j] * v[col[j]];
            }
            res[i] = val;
        }

        return res;
    }
};

// On génère les points sur [0, 1]²... 
std::vector<Node> setup_nodes(int N){
    std::vector<Node> nodes;
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
// maillage actuel : 
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

// Air orienté, moitié du produit vectoriel 
double triangle_area(const Node& A, const Node& B, const Node& C){
    double double_area = ((B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y));
    if (double_area < 0) std::cout << "Air négative, attention...\n";
    return 0.5 * double_area;
}

double error_l2(const std::vector<Node> &nodes, const std::vector<Element> &triangles,
    const std::vector<double> &U_h, 
    const std::function<double(double, double)> &sol_exa){
    double error = 0.;

    for (const auto& t : triangles){
        const Node &A = nodes[t.n1];
        const Node &B = nodes[t.n2];
        const Node &C = nodes[t.n3];
        const double x_c = (A.x + B.x + C.x) / 3.;
        const double y_c = (A.y + B.y + C.y) / 3.;
        const double u_h_c = (U_h[t.n1] + U_h[t.n2] + U_h[t.n3]) / 3.;
        const double e = u_h_c - sol_exa(x_c, y_c);
        error += std::abs(triangle_area(A, B, C)) * e * e;
    }

    return std::sqrt(error);
}

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
    std::ofstream file(file_name + ".vtk");
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

// Test pour voir coord des nodes, index triangle et air triangle
void test1() {
    int N = 2;
    std::vector<Node> nodes = setup_nodes(N);
    std::vector<Element> triangles = setup_triangles(N);

    std::cout << "Nodes:\n";
    for (int i = 0; i < nodes.size(); i++) {
        std::cout << i << " : (" << nodes[i].x << ", " << nodes[i].y << ")\n";
    }

    std::cout << "\nTriangles:\n";
    for (int i = 0; i < triangles.size(); i++) {
        std::cout << i << " : " << triangles[i].n1 << " " << triangles[i].n2 << " " << triangles[i].n3 << "\n";
    }

    std::cout << "\nAir Triangle\n";
    for(auto& t : triangles){
        double area = triangle_area(nodes[t.n1], nodes[t.n2], nodes[t.n3]);

        if (area <= 0){
            std::cout << "Triangle mal orienté !\n";
        } else{
            std::cout << "Air : " << area << "\n";
        }
    }
}


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
        if (nodes[i].is_on_boundary(x_limit, y_limit)){
            
            for (int j = 0; j < nodes.size(); j++){
                K[i][j] = 0.;
                K[j][i] = 0.;
            }
            K[i][i] = 1.;
            F[i] = 0.;
        }
    }
}

void display_mat(const std::vector<std::vector<double>> &K){
    int n_i = K.size();
    int n_j = K[0].size();

    for (int i = 0; i < n_i; i++){
        for (int j = 0; j <n_j; j++){
            std::cout << K[i][j] << "  ";
        }
        std::cout << "\n";
    }
}
void display_vec(const std::vector<double> &v){
    int n = v.size();

    for (int i = 0; i < n; i++){
        std::cout << v[i] << "\n";
    }
}

// Pivot de gauss
std::vector<double> solve_with_gauss(std::vector<std::vector<double>> K, std::vector<double> F){
    int n_ligne = K.size();
    int n_colonne = K[0].size();

    std::vector<double> U(n_colonne, 0.);

    // Parcourt les colones sauf la dernière
    for (int k = 0; k < n_colonne-1; k++){
        // Parcourt les lignes sous la diagonale
        for (int i = k + 1; i < n_ligne; i++){
            double f = K[i][k] / K[k][k];
            
            // On modifie toute la ligne i (on Ligne_i recoint - f * Ligne_k)
            // (donc on parcourt les colonnes)
            for (int j = 0; j < n_colonne; j++){
                K[i][j] -= f * K[k][j];
            }
            F[i] -= f * F[k];
        }
    }

    U[n_colonne - 1] = F[n_colonne-1] / K[n_colonne-1][n_colonne-1];
    for (int i = n_colonne - 2; i >= 0; i--){
        // Calcul de S
        double S = 0.;
        for (int j = i + 1; j < n_colonne; j++){
            S += K[i][j] * U[j];
        }

        U[i] = (F[i] - S) / K[i][i];
    }


    //display_mat(K);
    //display_vec(F);
    //display_vec(U);
    return U;
}

double dot(const std::vector<double>& u, const std::vector<double>& v){
    assert (u.size() == v.size());
    double res = 0;
    for (int i = 0; i < u.size(); i++){
        res += u[i] * v[i];
    }
    return res;
}
// axpy = alpha * x + y, on modifie y
void axpy(double alpha, const std::vector<double>& x, std::vector<double> &y){
    assert (x.size() == y.size());

    for (int i = 0; i< x.size(); i++){
        y[i] += alpha * x[i];
    }
}

// ax = alpha * x, on modifie x
void ax(double alpha, std::vector<double>& x){
    for (int i = 0; i < x.size(); i++){
        x[i] *= alpha;
    }
}


// Gradient conjugué
std::vector<double> solve_with_gc(const CSRMatrix &A, const std::vector<double> &b, int max_iter = 5){
    std::vector<double> x(b.size(), 0.);
    std::vector<double> r = b;
    std::vector<double> p = r;
    const double eps = 1.e-14;
    for (int k = 0; k < max_iter; k++){
        auto Ap = A.mult_with_v(p);
        double r_squared = dot(r, r);
        double p_dot_Ap = dot(p, Ap);
        if (std::abs(p_dot_Ap) < eps) 
            break;
        double alpha = r_squared / p_dot_Ap;
        // x = x + alpha * p
        axpy(alpha, p, x);
        // r = r - alpha * Ap
        axpy(-alpha, Ap, r); // modifie r_k pour r_next
        double r_next_squared = dot(r, r);
        if (r_next_squared < eps) {
            std::cout << "Iteration : " << k << " Residu : " << r_next_squared << "\n";
            break;
        }  
        double beta = r_next_squared / r_squared;
        // p = beta * p 
        ax(beta, p);
        // p = p + r
        axpy(1., r, p);
    }

    return x;
}

const double PI = 3.141592653589793;


void test2(){
    // On résoud : -laplacien(u) = f + dirichlet avec une méthode d'élément fini P1
    int N;
    std::cout << "Discrétisation N ?\n";
    std::cin >> N;
    std::vector<Node> nodes = setup_nodes(N);
    std::vector<Element> triangles = setup_triangles(N);
    std::function<double(double, double)> f = [](double x, double y)->double{ 
        return 2 * PI * PI * sin(PI * x) * sin(PI * y);};
    std::function<double(double, double)> solution_exacte = [](double x, double y)->double{
        return sin(PI * x) * sin(PI*y);
    };
    auto K = get_global_mat(nodes, triangles);
    auto F = get_global_f(nodes, triangles, f);
    set_dirichlet_condition(nodes, K, F);
    //display_mat(K);
    //display_vec(F);
    std::vector<double> U = solve_with_gauss(K, F);
    //display_mat(K);
    //display_vec(F);
    std::string file_name = "solution_comparer_N" + std::to_string(N);
    export_for_paraview(nodes, triangles, U, solution_exacte, file_name);
    
    double err = error_l2(nodes, triangles, U, solution_exacte);
    std::cout << "Error L2 (N=" + std::to_string(N) + ") : " << err << "\n";
}

void test3(){
    CSRMatrix mat = CSRMatrix({{1, 0, 2}, {0, 3, 0}, {4, 0, 5}});
    std::vector<double> vect = {1, 2, 3};
    std::cout << "Info Matrice\n";
    for (const auto& v : mat.row){
        std::cout << v << ", ";
    }
    std::cout << "\n";
    for (const auto& v : mat.col){
        std::cout << v << ", ";
    }
    std::cout << "\n";
    for (const auto& v : mat.values){
        std::cout << v << ", ";
    }
    std::cout << "\n";
    display_vec(mat.mult_with_v(vect));
    //test2();
}

void test4(){
    CSRMatrix mat = CSRMatrix({{4, 1}, {1, 3}});
    std::vector<double> b = {{1}, {2}};
    std::vector<double> x = solve_with_gc(mat, b);
    display_vec(x);
}

int main(){
    test4();
    return 0;
}
