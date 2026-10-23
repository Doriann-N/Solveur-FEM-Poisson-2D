#pragma once

#include <functional>
#include "mesh_data.hpp"
#include "csr_matrix.hpp"

/* 
'Solveur' de -Δ(u) = f + Dirichlet sur le bord de [0, 1]² par élément 
fini P1 sur un maillage triangulaire alterné
*/

class Simulation{
    private:
        /* Discrétisation du maillage sur un axe */
        int N;
        /* Matrice du système Ax = b */
        CSRMatrix A_csr;
        /* Membre de droite du système Ax = b */
        std::vector<double> b;
        std::vector<Node> nodes;
        std::vector<Element> triangles;
        /* Solution calculée numériquement */
        std::vector<double> solution;
        /* Fonction du membre de droite dans l'équation */
        std::function<double(double, double)> rhs_function;
        /* Fonction solution du problème si on l'a connait */
        std::function<double(double, double)> solution_function;
    public:
        void solve(int max_iter, bool verbose);
        void save_vtk(const std::string& p_filename);
        /* Permet de définir la fonction solution du problème */
        void set_solution_function(std::function<double(double, double)> solution_exacte);
        Simulation(const int N, std::function<double(double, double)> rhs_function, bool alternate_triangle_mesh);
        const std::vector<Element>& get_triangles() const;
        const std::vector<double>& get_solution() const;
        const std::vector<Node>& get_nodes() const;
        const CSRMatrix& get_matrix() const;
        const std::vector<double>& get_rhs_vector() const;
        const float get_l2_error() const;
};