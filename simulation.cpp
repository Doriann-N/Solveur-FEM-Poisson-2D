#include "simulation.hpp"
#include "fem.hpp"
#include "solver.hpp"
#include "error.hpp"
#include "vtk_writer.hpp"

Simulation::Simulation(const int p_N, std::function<double(double, double)> p_rhs_function, bool alternate_triangle_mesh){
    N = p_N;
    rhs_function = p_rhs_function;
    nodes = setup_nodes(N);
    triangles = setup_triangles(N, alternate_triangle_mesh);
}

void Simulation::solve(int max_iter, bool verbose){
    /* 
    Calcul de matrice d'assemblage pour le 
    système Ax = b (avant A = K, x = U, et b = F)
    Ici j'adopte une autre convention d'écriture 
    */
    auto A_list = get_global_mat_in_list_form(nodes, triangles);
    b = get_global_f(nodes, triangles, rhs_function);
    
    /* Condition de Dirichlet */
    set_dirichlet_condition(nodes, A_list, b);
    
    /* Matrice CSR du système à résoudre */
    A_csr = CSRMatrix(nodes.size(), nodes.size(), A_list);

    solution = solve_with_gc(A_csr, b, max_iter, verbose);
}

/* Pour sauvegarder la solution pour être lue dans ParaView. */
void Simulation::save_vtk(const std::string& p_filename){
    /* Afin de s'assurer qu'on a calculer une solution, on pourrait faire mieux avec un booléen */
    assert(!solution.empty());
    std::string file_name = "Solutions/" + p_filename;

    VTK_Writer writer(file_name);
    writer.set_mesh(nodes, triangles);
    writer.add_scalar_field("U_h", solution);

    if (solution_function){
        std::vector<double> solution_exacte = sample_function_on_nodes(nodes, solution_function);
        std::vector<double> errs = errors(solution, solution_exacte);
        writer.add_scalar_field("U_exa", solution_exacte);
        writer.add_scalar_field("Errors", errs);
    }
    writer.write();
}

/* Permet de définir une solution exacte afin de pouvoir calculer l'erreur en norme L2 par exemple */
void Simulation::set_solution_function(std::function<double(double, double)> solution_exacte){
    solution_function = solution_exacte;
}

/* Pour obtenir la norme L2 de l'erreur, il faut la solution exacte. */
const float Simulation::get_l2_error(bool print_error) const{
    /* On vérifie que la solution a bien été calculer */
    assert(!solution.empty());
    /* On vérifie qu'on a bien la solution exacte */
    assert(solution_function);

    double err = error_l2(nodes, triangles, solution, solution_function);
    if (print_error) std::cout << "Erreur L2 : " << err << "\n";
    return err;
}

/* Getters */

const CSRMatrix& Simulation::get_matrix() const{
    return A_csr;
}
const std::vector<Node>& Simulation::get_nodes() const{
    return nodes;
}
const std::vector<Element>& Simulation::get_triangles() const{
    return triangles;
}
const std::vector<double>& Simulation::get_solution() const{
    return solution;
}
const std::vector<double>& Simulation::get_rhs_vector() const{
    return b;
}