#include "solver.hpp"

namespace NavierStokes
{

/**
 * Run the incompressible navier stokes solver
 */
template <int dim>
void IncompressibleNavierStokes<dim>::run()
{
    printf("Run the simulation \n") ;
}
template void IncompressibleNavierStokes<ELEMENT_DIM>::run() ;


/**
 * Constructor
 */
template <int dim>
IncompressibleNavierStokes<dim>::IncompressibleNavierStokes(
    const RunTimeParameters::DataStorage&   data
)
:
meshStorage( data )
{
    cout << "finishing constructor" << endl ;
}
template IncompressibleNavierStokes<ELEMENT_DIM>::IncompressibleNavierStokes( 
   const RunTimeParameters::DataStorage&   data 
) ;


}