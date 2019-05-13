#include "solver.hpp"
#include "../boundary/boundary.hpp"
#include "../preconditioner/preconditioner.hpp"

namespace IncompNS
{
using namespace dealii;

/**
 *  constructor of the Stationary Navier Stokes solver 
 */
template <int dim>
StationaryNavierStokes<dim>::StationaryNavierStokes
(
    const unsigned int degree
)
    : viscosity(1.0 / 7500.0)
    , gamma(1.0)
    , degree(degree)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1)
    , dof_handler(triangulation)
{}
  
/**
 * run the simulation
 */
template <int dim>
void StationaryNavierStokes<dim>::run
(
    const unsigned int refinement
)
{
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(5);
    const double Re = 1.0 / viscosity;
    if (Re > 1000.0)
    {
        std::cout << "Searching for initial guess ..." << std::endl;
        const double step_size = 2000.0;
        compute_initial_guess(step_size);
        std::cout << "Found initial guess." << std::endl;
        std::cout << "Computing solution with target Re = " << Re << std::endl;
        viscosity = 1.0 / Re;
        newton_iteration(1e-12, 50, refinement, false, true);
    }
    else
    {
        newton_iteration(1e-12, 50, refinement, true, true);
    }
}

//to avoid linker issue
template class StationaryNavierStokes<ELEMENT_DIM> ;

} // namespace IncompNS
