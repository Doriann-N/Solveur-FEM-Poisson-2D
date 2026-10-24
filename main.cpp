#include "mesh_data.hpp"
#include "csr_matrix.hpp"
#include "fem.hpp"
#include "error.hpp"
#include "linalgebra.hpp"
#include "solver.hpp"
#include "vtk_writer.hpp"
#include "simulation.hpp"
#include <iomanip>

const double PI = 3.141592653589793;


void check_l2_convergence(std::function<double(double, double)> f, std::function<double(double, double)> solution_exacte){
    std::vector<int> Ns = {10, 20, 40, 80, 160, 320, 640, 1280};
    std::vector<double> errs;
    std::vector<double> times;
    for(const int N : Ns){
        bool alternate_triangle = true;
        Simulation sim(N, f, alternate_triangle);
        sim.set_solution_function(solution_exacte);
        double start_time = omp_get_wtime();
        sim.solve(50, false);
        double time = omp_get_wtime() - start_time;
        double err = sim.get_l2_error();
        errs.push_back(err);
        times.push_back(time);
    }
    int last = errs.size()-1;
    std::cout << "----------------------------------------------------------------------------------\n";
    std::cout << "N | Erreur L2 | Temps (s) | Max Threads | Rapport d'erreurs | Ordre expérimental\n";
    for (int i = 0; i <= last; i++){
        std::cout << Ns[i] << " | ";
        std::cout << std::scientific << errs[i] << " | ";
        std::cout << std::fixed << times[i] << " | " << omp_get_max_threads();
        if (i != last){
            std::cout << " | " << errs[i] / errs[i+1];
            double p = std::log(errs[i] / errs[i+1]) / std::log(2.0);
            std::cout << " | " << p;
        }
        std::cout << "\n";
    }
    std::cout << "----------------------------------------------------------------------------------\n";
}


int main(){
    /* Membre de droite dans -Δ(u) = f */ 
    std::function<double(double, double)> f = [](double x, double y)->double{ 
        return 2. * PI * PI * sin(PI * x) * sin(PI * y);
    };
    
    /* Solution exacte */
    std::function<double(double, double)> solution_exacte = [](double x, double y)->double{
        return sin(PI * x) * sin(PI * y);
    };
    
    /* Discrétisation du maillage */
    int N = ask_for_discretisation();
    /* Alterner les triangles du maillages ? Un maillage 
    alterné résulte en une convergence en 2 étapes du gradient 
    conjugué, je ne sais pas si c'est lié très exactement à ça */
    bool alternate_triangle = ask_for_boolean("Alterner les triangles pour le mesh ?");
    /* Demande de nom pour sauvegarder le fichier, vide 
    impliquera aucune sauvegarde */
    std::string filename = ask_for_name_file();
    /* Prépare le maillage */
    Simulation sim(N, f, alternate_triangle);
    /* Indique si on va print des info lié à la 
    convergence avec le gradient conjugué */
    bool verbose = true;
    sim.solve(50, verbose);
    /* On définit la solution exacte pour calculer l'erreur */
    sim.set_solution_function(solution_exacte);
    float err = sim.get_l2_error();
    
    if (!filename.empty()) sim.save_vtk(filename);
    //print_info_nodes(sim.get_nodes());
    //const auto& A = sim.get_matrix();
    //const auto& b = sim.get_rhs_vector();

    /* Pour voir le tableau comparatif des convergences */
    bool show_convergence = ask_for_boolean("Montrer l'étude de convergence ?");
    if (show_convergence){
        std::cout << "Étude de la convergence :\n";
        check_l2_convergence(f, solution_exacte);
    }

    return 0;
}