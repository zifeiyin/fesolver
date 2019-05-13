#include "../solver/solver.hpp"
#include "../boundary/boundary.hpp"

namespace IncompNS
{
	using namespace dealii;

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

    template void StationaryNavierStokes<ELEMENT_DIM>::setup_dofs() ;

  		
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

    template void StationaryNavierStokes<ELEMENT_DIM>::initialize_system() ;


}