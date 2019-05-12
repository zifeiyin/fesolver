#include "solver.hpp"
#include "../boundary/boundary.hpp"
#include "../preconditioner/preconditioner.hpp"

namespace IncompNS
{
	using namespace dealii;


  	template <int dim>
  	StationaryNavierStokes<dim>::StationaryNavierStokes(const unsigned int degree)
    	: viscosity(1.0 / 7500.0)
    	, gamma(1.0)
    	, degree(degree)
    	, triangulation(Triangulation<dim>::maximum_smoothing)
    	, fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1)
    	, dof_handler(triangulation)
  	{}
  
  	template <int dim>
  	void StationaryNavierStokes<dim>::setup_dofs()
  	{
    	system_matrix.clear();
    	pressure_mass_matrix.clear();
    	dof_handler.distribute_dofs(fe);
    	std::vector<unsigned int> block_component(dim + 1, 0);
    	block_component[dim] = 1;
    	DoFRenumbering::component_wise(dof_handler, block_component);
    	dofs_per_block.resize(2);
    	DoFTools::count_dofs_per_block(dof_handler,
    	                               dofs_per_block,
    	                               block_component);
    	unsigned int dof_u = dofs_per_block[0];
    	unsigned int dof_p = dofs_per_block[1];
    	FEValuesExtractors::Vector velocities(0);
    	{
    		nonzero_constraints.clear();
      		DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
      		VectorTools::interpolate_boundary_values(dof_handler,
     		                                          0,
     		                                          BoundaryValues<dim>(),
     		                                          nonzero_constraints,
     		                                          fe.component_mask(velocities));
    	}
    	nonzero_constraints.close();
    	{
      		zero_constraints.clear();
      		DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
      		VectorTools::interpolate_boundary_values(dof_handler,
      		                                         0,
      		                                         Functions::ZeroFunction<dim>(
      		                                           dim + 1),
      		                                         zero_constraints,
      		                                         fe.component_mask(velocities));
    	}
    	zero_constraints.close();
    	std::cout << "Number of active cells: " << triangulation.n_active_cells()
    	          << std::endl
    	          << "Number of degrees of freedom: " << dof_handler.n_dofs()
    	          << " (" << dof_u << " + " << dof_p << ')' << std::endl;
  		}
  		
	template <int dim>
  	void StationaryNavierStokes<dim>::initialize_system()
  	{
    	{
      	BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
      	DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
      	sparsity_pattern.copy_from(dsp);
    }
    system_matrix.reinit(sparsity_pattern);
    present_solution.reinit(dofs_per_block);
    newton_update.reinit(dofs_per_block);
    system_rhs.reinit(dofs_per_block);
  }
  template <int dim>
  void StationaryNavierStokes<dim>::assemble(const bool initial_step,
                                             const bool assemble_matrix)
  {
    if (assemble_matrix)
      system_matrix = 0;
    system_rhs = 0;
    QGauss<dim> quadrature_formula(degree + 2);
    FEValues<dim> fe_values(fe,
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
  template <int dim>
  void StationaryNavierStokes<dim>::assemble_system(const bool initial_step)
  {
    assemble(initial_step, true);
  }
  template <int dim>
  void StationaryNavierStokes<dim>::assemble_rhs(const bool initial_step)
  {
    assemble(initial_step, false);
  }
  template <int dim>
  void StationaryNavierStokes<dim>::solve(const bool initial_step)
  {
    const AffineConstraints<double> &constraints_used =
      initial_step ? nonzero_constraints : zero_constraints;
    SolverControl solver_control(system_matrix.m(),
                                 1e-4 * system_rhs.l2_norm(),
                                 true);
    SolverFGMRES<BlockVector<double>> gmres(solver_control);
    SparseILU<double>                 pmass_preconditioner;
    pmass_preconditioner.initialize(pressure_mass_matrix,
                                    SparseILU<double>::AdditionalData());
    const BlockSchurPreconditioner<SparseILU<double>> preconditioner(
      gamma,
      viscosity,
      system_matrix,
      pressure_mass_matrix,
      pmass_preconditioner);
    gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
    std::cout << "FGMRES steps: " << solver_control.last_step() << std::endl;
    constraints_used.distribute(newton_update);
  }
  template <int dim>
  void StationaryNavierStokes<dim>::refine_mesh()
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    FEValuesExtractors::Vector velocity(0);
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      QGauss<dim - 1>(degree + 1),
      std::map<types::boundary_id, const Function<dim> *>(),
      present_solution,
      estimated_error_per_cell,
      fe.component_mask(velocity));
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.0);
    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
    triangulation.execute_coarsening_and_refinement();
    setup_dofs();
    BlockVector<double> tmp(dofs_per_block);
    solution_transfer.interpolate(present_solution, tmp);
    nonzero_constraints.distribute(tmp);
    initialize_system();
    present_solution = tmp;
  }
 


  template <int dim>
  void StationaryNavierStokes<dim>::compute_initial_guess(double step_size)
  {
    const double target_Re = 1.0 / viscosity;
    bool is_initial_step = true;
    for (double Re = 1000.0; Re < target_Re;
         Re        = std::min(Re + step_size, target_Re))
      {
        viscosity = 1.0 / Re;
        std::cout << "Searching for initial guess with Re = " << Re
                  << std::endl;
        newton_iteration(1e-12, 50, 0, is_initial_step, false);
        is_initial_step = false;
      }
  }


  template <int dim>
  void StationaryNavierStokes<dim>::run(const unsigned int refinement)
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
  	
  template class StationaryNavierStokes<ELEMENT_DIM> ;

} // namespace IncompNS
