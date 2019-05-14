#include "../solver/solver.hpp"

namespace NavierStokes
{

template <int dim>
void NavierStokesProjection<dim>::output_results (const unsigned int step)
{
    assemble_vorticity ( (step == 1));
    const FESystem<dim> joint_fe (fe_velocity, dim,
                                  fe_pressure, 1,
                                  fe_velocity, 1);
    DoFHandler<dim> joint_dof_handler (triangulation);
    joint_dof_handler.distribute_dofs (joint_fe);
    Assert (joint_dof_handler.n_dofs() ==
            ((dim + 1)*dof_handler_velocity.n_dofs() +
             dof_handler_pressure.n_dofs()),
            ExcInternalError());
    Vector<double> joint_solution (joint_dof_handler.n_dofs());
    std::vector<types::global_dof_index> loc_joint_dof_indices (joint_fe.dofs_per_cell),
        loc_vel_dof_indices (fe_velocity.dofs_per_cell),
        loc_pres_dof_indices (fe_pressure.dofs_per_cell);
    typename DoFHandler<dim>::active_cell_iterator
    joint_cell = joint_dof_handler.begin_active(),
    joint_endc = joint_dof_handler.end(),
    vel_cell   = dof_handler_velocity.begin_active(),
    pres_cell  = dof_handler_pressure.begin_active();
    for (; joint_cell != joint_endc; ++joint_cell, ++vel_cell, ++pres_cell)
      {
        joint_cell->get_dof_indices (loc_joint_dof_indices);
        vel_cell->get_dof_indices (loc_vel_dof_indices);
        pres_cell->get_dof_indices (loc_pres_dof_indices);
        for (unsigned int i=0; i<joint_fe.dofs_per_cell; ++i)
          switch (joint_fe.system_to_base_index(i).first.first)
            {
            case 0:
              Assert (joint_fe.system_to_base_index(i).first.second < dim,
                      ExcInternalError());
              joint_solution (loc_joint_dof_indices[i]) =
                u_n[ joint_fe.system_to_base_index(i).first.second ]
                (loc_vel_dof_indices[ joint_fe.system_to_base_index(i).second ]);
              break;
            case 1:
              Assert (joint_fe.system_to_base_index(i).first.second == 0,
                      ExcInternalError());
              joint_solution (loc_joint_dof_indices[i]) =
                pres_n (loc_pres_dof_indices[ joint_fe.system_to_base_index(i).second ]);
              break;
            case 2:
              Assert (joint_fe.system_to_base_index(i).first.second == 0,
                      ExcInternalError());
              joint_solution (loc_joint_dof_indices[i]) =
                rot_u (loc_vel_dof_indices[ joint_fe.system_to_base_index(i).second ]);
              break;
            default:
              Assert (false, ExcInternalError());
            }
      }
    std::vector<std::string> joint_solution_names (dim, "v");
    joint_solution_names.emplace_back("p");
    joint_solution_names.emplace_back("rot_u");
    DataOut<dim> data_out;
    data_out.attach_dof_handler (joint_dof_handler);
    std::vector< DataComponentInterpretation::DataComponentInterpretation >
    component_interpretation (dim+2,
                              DataComponentInterpretation::component_is_part_of_vector);
    component_interpretation[dim]
      = DataComponentInterpretation::component_is_scalar;
    component_interpretation[dim+1]
      = DataComponentInterpretation::component_is_scalar;
    data_out.add_data_vector (joint_solution,
                              joint_solution_names,
                              DataOut<dim>::type_dof_data,
                              component_interpretation);
    data_out.build_patches (deg + 1);
    std::ofstream output (("solution-" +
                           Utilities::int_to_string (step, 5) +
                           ".vtk").c_str());
    data_out.write_vtk (output);
}
template void NavierStokesProjection<ELEMENT_DIM>::output_results (const unsigned int step) ;


}