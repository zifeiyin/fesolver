#include "../solver/solver.hpp"

namespace IncompNS
{
using namespace dealii ;


/**
 * compute initial guess for the solution
 */
template <int dim>
void StationaryNavierStokes<dim>::compute_initial_guess
(
    double step_size
)
{
    const double target_Re = 1.0 / viscosity;
    bool is_initial_step = true;
    for ( double Re = 1000.0; Re < target_Re; Re = std::min(Re + step_size, target_Re) )
    {
        viscosity = 1.0 / Re;
        std::cout << "Searching for initial guess with Re = " << Re
                  << std::endl;
        newton_iteration(1e-12, 50, 0, is_initial_step, false);
        is_initial_step = false;
    }
}
template void StationaryNavierStokes<ELEMENT_DIM>::compute_initial_guess(double step_size) ;

}