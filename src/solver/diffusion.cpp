#include "../solver/solver.hpp"

namespace NavierStokes
{

template <int dim>
void NavierStokesProjection<dim>::diffusion_step (const bool reinit_prec)
{
    pres_tmp.equ (-1., pres_n);
    pres_tmp.add (-4./3., phi_n, 1./3., phi_n_minus_1);
    assemble_advection_term();
    for (unsigned int d=0; d<dim; ++d)
      {
        force[d] = 0.;
        v_tmp.equ (2./dt,u_n[d]);
        v_tmp.add (-.5/dt,u_n_minus_1[d]);
        vel_Mass.vmult_add (force[d], v_tmp);
        pres_Diff[d].vmult_add (force[d], pres_tmp);
        u_n_minus_1[d] = u_n[d];
        vel_it_matrix[d].copy_from (vel_Laplace_plus_Mass);
        vel_it_matrix[d].add (1., vel_Advection);
        vel_exact.set_component(d);
        boundary_values.clear();
        for (std::vector<types::boundary_id>::const_iterator
             boundaries = boundary_ids.begin();
             boundaries != boundary_ids.end();
             ++boundaries)
          {
            switch (*boundaries)
              {
              case 1:
                VectorTools::
                interpolate_boundary_values (dof_handler_velocity,
                                             *boundaries,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
                break;
              case 2:
                VectorTools::
                interpolate_boundary_values (dof_handler_velocity,
                                             *boundaries,
                                             vel_exact,
                                             boundary_values);
                break;
              case 3:
                if (d != 0)
                  VectorTools::
                  interpolate_boundary_values (dof_handler_velocity,
                                               *boundaries,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values);
                break;
              case 4:
                VectorTools::
                interpolate_boundary_values (dof_handler_velocity,
                                             *boundaries,
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        MatrixTools::apply_boundary_values (boundary_values,
                                            vel_it_matrix[d],
                                            u_n[d],
                                            force[d]);
      }
    Threads::TaskGroup<void> tasks;
    for (unsigned int d=0; d<dim; ++d)
      {
        if (reinit_prec)
          prec_velocity[d].initialize (vel_it_matrix[d],
                                       SparseILU<double>::
                                       AdditionalData (vel_diag_strength,
                                                       vel_off_diagonals));
        tasks += Threads::new_task (&NavierStokesProjection<dim>::
                                    diffusion_component_solve,
                                    *this, d);
      }
    tasks.join_all();
}
template void NavierStokesProjection<ELEMENT_DIM>::diffusion_step (const bool reinit_prec) ;




template <int dim>
void NavierStokesProjection<dim>::diffusion_component_solve (const unsigned int d)
{
    SolverControl solver_control (vel_max_its, vel_eps*force[d].l2_norm());
    SolverGMRES<> gmres (solver_control,
                         SolverGMRES<>::AdditionalData (vel_Krylov_size));
    gmres.solve (vel_it_matrix[d], u_n[d], force[d], prec_velocity[d]);
}
template void NavierStokesProjection<ELEMENT_DIM>::diffusion_component_solve (const unsigned int d) ;



}