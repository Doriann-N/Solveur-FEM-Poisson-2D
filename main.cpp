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
        Simulation sim(N, f, false);
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
    

    int N = ask_for_discretisation();
    Simulation sim(N, f, true);
    sim.solve(50, true);
    sim.set_solution_function(solution_exacte);
    float err = sim.get_l2_error();
    
    //sim.save_vtk("Other_End_solution_N_" + std::to_string(N));
    //print_info_nodes(sim.get_nodes());
    //const auto& A = sim.get_matrix();
    //const auto& b = sim.get_rhs_vector();

    //check_l2_convergence(f, solution_exacte);

    return 0;
}