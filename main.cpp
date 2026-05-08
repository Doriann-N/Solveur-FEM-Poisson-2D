#include "mesh_data.hpp"
#include "csr_matrix.hpp"
#include "fem.hpp"
#include "error.hpp"
#include "linalgebra.hpp"
#include "solver.hpp"
#include "vtk_writer.hpp"

const double PI = 3.141592653589793;

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
    for (int i = 0; i < n; i++){ std::cout << v[i] << "\n";}
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

/* // Test pour voir coord des nodes, index triangle et air triangle
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
    std::vector<double> x = solve_with_gc(mat, b, 5);
    display_vec(x);
} */

void check_l2_convergence(){
    // Membre de droite de l'équation
    std::function<double(double, double)> f = [](double x, double y)->double{ 
        return 2 * PI * PI * sin(PI * x) * sin(PI * y);
    };
    // Solution exacte
    std::function<double(double, double)> solution_exacte = [](double x, double y)->double{
        return sin(PI * x) * std::sin(PI*y);
    };
    
    std::vector<int> Ns = {5, 10, 20, 40, 80, 160};
    std::vector<double> errs;
    for(const int N : Ns){
        std::vector<Node> nodes = setup_nodes(N);
        std::vector<Element> triangles = setup_triangles(N);
        auto K = get_global_mat(nodes, triangles);
        auto F = get_global_f(nodes, triangles, f);
        set_dirichlet_condition(nodes, K, F);
        CSRMatrix csr_mat_K = CSRMatrix(K);
        std::vector<double> U = solve_with_gc(csr_mat_K, F, 5, false);
        double err = error_l2(nodes, triangles, U, solution_exacte);
        //std::cout << "Error L2 (N=" + std::to_string(N) + ") : " << err << "\n";
        errs.push_back(err);
    }

    int last = errs.size()-1;
    std::cout<< "N | L2 Error | Errors[i] / Errors[i+1]\n";
    for (int i = 0; i <= last; i++){
        std::cout << Ns[i] << " | " << errs[i];
        if (i != last){
            std::cout << " | " << errs[i] / errs[i+1];
        }
        std::cout << "\n";
    }

}


/* int main(){
    int N;
    std::cout << "Discrétisation N ?\n";
    std::cin >> N;
    std::vector<Node> nodes = setup_nodes(N);
    std::vector<Element> triangles = setup_triangles(N);


    std::function<double(double, double)> f = [](double x, double y)->double{ 
        return 2 * PI * PI * sin(PI * x) * sin(PI * y);
    };
    std::function<double(double, double)> solution_exacte = [](double x, double y)->double{
        return sin(PI * x) * std::sin(PI*y);
    };
    
    
    auto K = get_global_mat(nodes, triangles);
    auto F = get_global_f(nodes, triangles, f);
    bool verbose = false;
    if (verbose){
        display_mat(K);
        printf("---\n");
        display_vec(F);
    }   
    set_dirichlet_condition(nodes, K, F);
    if (verbose){
        printf("Condition de dirichlet appliquée\n");
        display_mat(K);
        printf("---\n");
        display_vec(F);
    } 
    CSRMatrix csr_mat_K = CSRMatrix(K);
    //std::vector<double> U = solve_with_gauss(K, F);
    std::vector<double> U = solve_with_gc(csr_mat_K, F, 5, true);
    //export_for_paraview(nodes, triangles, U, solution_exacte, file_name);
    
    double err = error_l2(nodes, triangles, U, solution_exacte);
    std::cout << "Error L2 (N=" + std::to_string(N) + ") : " << err << "\n";
    //if (verbose) display_vec(U);

    std::vector<double> U_sol = sample_function_on_nodes(nodes, solution_exacte);
    std::vector<double> errs = errors(U, U_sol);
    
    std::string file_name = "Solutions/solution_comparer_N" + std::to_string(N) + " new";

    VTK_Writer writer(file_name);
    writer.set_mesh(nodes, triangles);
    writer.add_scalar_field("U_h", U);
    writer.add_scalar_field("U_exa", U_sol);
    writer.add_scalar_field("Errors", errs);
    writer.write();

    return 0;
} */

int main(){
    check_l2_convergence();
}