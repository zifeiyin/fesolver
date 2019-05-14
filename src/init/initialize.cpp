#include "../solver/solver.hpp"

namespace NavierStokes
{

template <int dim>
void NavierStokesProjection<dim>::initialize()
{
    vel_Laplace_plus_Mass = 0.;
    vel_Laplace_plus_Mass.add (1./Re, vel_Laplace);
    vel_Laplace_plus_Mass.add (1.5/dt, vel_Mass);
    EquationData::Pressure<dim> pres (t_0);
    VectorTools::interpolate (dof_handler_pressure, pres, pres_n_minus_1);
    pres.advance_time (dt);
    VectorTools::interpolate (dof_handler_pressure, pres, pres_n);
    phi_n = 0.;
    phi_n_minus_1 = 0.;
    for (unsigned int d=0; d<dim; ++d)
      {
        vel_exact.set_time (t_0);
        vel_exact.set_component(d);
        VectorTools::interpolate (dof_handler_velocity, Functions::ZeroFunction<dim>(), u_n_minus_1[d]);
        vel_exact.advance_time (dt);
        VectorTools::interpolate (dof_handler_velocity, Functions::ZeroFunction<dim>(), u_n[d]);
      }
}
template void NavierStokesProjection<ELEMENT_DIM>::initialize() ;

template <int dim>
void NavierStokesProjection<dim>::initialize_gradient_operator()
{
    {
      DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(), dof_handler_pressure.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler_velocity, dof_handler_pressure, dsp);
      sparsity_pattern_pres_vel.copy_from (dsp);
    }
    InitGradPerTaskData per_task_data (0, fe_velocity.dofs_per_cell,
                                       fe_pressure.dofs_per_cell);
    InitGradScratchData scratch_data (fe_velocity,
                                      fe_pressure,
                                      quadrature_velocity,
                                      update_gradients | update_JxW_values,
                                      update_values);
    for (unsigned int d=0; d<dim; ++d)
      {
        pres_Diff[d].reinit (sparsity_pattern_pres_vel);
        per_task_data.d = d;
        WorkStream::run (IteratorPair (IteratorTuple (dof_handler_velocity.begin_active(),
                                                      dof_handler_pressure.begin_active()
                                                     )
                                      ),
                         IteratorPair (IteratorTuple (dof_handler_velocity.end(),
                                                      dof_handler_pressure.end()
                                                     )
                                      ),
                         *this,
                         &NavierStokesProjection<dim>::assemble_one_cell_of_gradient,
                         &NavierStokesProjection<dim>::copy_gradient_local_to_global,
                         scratch_data,
                         per_task_data
                        );
      }
}
template void NavierStokesProjection<ELEMENT_DIM>::initialize_gradient_operator() ;


}