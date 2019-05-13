#include "../solver/solver.hpp"

namespace IncompNS
{
    using namespace dealii ;

template <int dim>
void StationaryNavierStokes<dim>::assemble
( 
    const bool initial_step,
    const bool assemble_matrix 
){
    if (assemble_matrix)    system_matrix = 0;
    
    system_rhs = 0;
    
    QGauss<dim>     quadrature_formula(degree + 2);
    FEValues<dim>   fe_values(  fe,
                                quadrature_formula,
                                update_values | update_quadrature_points |
                                update_JxW_values | update_gradients);
    
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
    std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
    std::vector<double>         present_pressure_values(n_q_points);
    std::vector<double>         div_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double>         phi_p(dofs_per_cell);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        local_matrix = 0;
        local_rhs    = 0;
        fe_values[velocities].get_function_values(evaluation_point,
                                                  present_velocity_values);
        fe_values[velocities].get_function_gradients(
          evaluation_point, present_velocity_gradients);
        fe_values[pressure].get_function_values(evaluation_point,
                                                present_pressure_values);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                phi_u[k]      = fe_values[velocities].value(k, q);
                phi_p[k]      = fe_values[pressure].value(k, q);
              }
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                if (assemble_matrix)
                  {
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        local_matrix(i, j) +=
                          (viscosity *
                             scalar_product(grad_phi_u[j], grad_phi_u[i]) +
                           present_velocity_gradients[q] * phi_u[j] * phi_u[i] +
                           grad_phi_u[j] * present_velocity_values[q] *
                             phi_u[i] -
                           div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                           gamma * div_phi_u[j] * div_phi_u[i] +
                           phi_p[i] * phi_p[j]) *
                          fe_values.JxW(q);
                      }
                  }
                double present_velocity_divergence =
                  trace(present_velocity_gradients[q]);
                local_rhs(i) +=
                  (-viscosity * scalar_product(present_velocity_gradients[q],
                                               grad_phi_u[i]) -
                   present_velocity_gradients[q] * present_velocity_values[q] *
                     phi_u[i] +
                   present_pressure_values[q] * div_phi_u[i] +
                   present_velocity_divergence * phi_p[i] -
                   gamma * present_velocity_divergence * div_phi_u[i]) *
                  fe_values.JxW(q);
              }
          }
        cell->get_dof_indices(local_dof_indices);
        const AffineConstraints<double> &constraints_used =
          initial_step ? nonzero_constraints : zero_constraints;
        if (assemble_matrix)
          {
            constraints_used.distribute_local_to_global(local_matrix,
                                                        local_rhs,
                                                        local_dof_indices,
                                                        system_matrix,
                                                        system_rhs);
          }
        else
          {
            constraints_used.distribute_local_to_global(local_rhs,
                                                        local_dof_indices,
                                                        system_rhs);
          }
      }
    if (assemble_matrix)
      {
        pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
        pressure_mass_matrix.copy_from(system_matrix.block(1, 1));
        system_matrix.block(1, 1) = 0;
      }
  }

    template void StationaryNavierStokes<2>::assemble(  const bool initial_step,
                                                        const bool assemble_matrix  ) ;


    template <int dim>
    void StationaryNavierStokes<dim>::assemble_system(const bool initial_step)
    {
        assemble(initial_step, true);
    }
    template void StationaryNavierStokes<ELEMENT_DIM>::assemble_system(const bool initial_step) ;
    

    template <int dim>
    void StationaryNavierStokes<dim>::assemble_rhs(const bool initial_step)
    {
        assemble(initial_step, false);
    }

    template void StationaryNavierStokes<ELEMENT_DIM>::assemble_rhs(const bool initial_step) ;

}