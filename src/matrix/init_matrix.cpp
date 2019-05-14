#include "../solver/solver.hpp"

namespace NavierStokes
{

template <int dim> 
void NavierStokesProjection<dim>::initialize_velocity_matrices()
{
    {
      DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(), dof_handler_velocity.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler_velocity, dsp);
      sparsity_pattern_velocity.copy_from (dsp);
    }
    vel_Laplace_plus_Mass.reinit (sparsity_pattern_velocity);
    for (unsigned int d=0; d<dim; ++d)
      vel_it_matrix[d].reinit (sparsity_pattern_velocity);
    vel_Mass.reinit (sparsity_pattern_velocity);
    vel_Laplace.reinit (sparsity_pattern_velocity);
    vel_Advection.reinit (sparsity_pattern_velocity);
    MatrixCreator::create_mass_matrix (dof_handler_velocity,
                                       quadrature_velocity,
                                       vel_Mass);
    MatrixCreator::create_laplace_matrix (dof_handler_velocity,
                                          quadrature_velocity,
                                          vel_Laplace);
}
template void NavierStokesProjection<ELEMENT_DIM>::initialize_velocity_matrices() ;


template <int dim>
void NavierStokesProjection<dim>::initialize_pressure_matrices()
{
    {
      DynamicSparsityPattern dsp(dof_handler_pressure.n_dofs(), dof_handler_pressure.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler_pressure, dsp);
      sparsity_pattern_pressure.copy_from (dsp);
    }
    pres_Laplace.reinit (sparsity_pattern_pressure);
    pres_iterative.reinit (sparsity_pattern_pressure);
    pres_Mass.reinit (sparsity_pattern_pressure);
    MatrixCreator::create_laplace_matrix (dof_handler_pressure,
                                          quadrature_pressure,
                                          pres_Laplace);
    MatrixCreator::create_mass_matrix (dof_handler_pressure,
                                       quadrature_pressure,
                                       pres_Mass);
}
template void NavierStokesProjection<ELEMENT_DIM>::initialize_pressure_matrices() ;


}
