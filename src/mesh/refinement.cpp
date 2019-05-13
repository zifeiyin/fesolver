#include "../solver/solver.hpp"

namespace IncompNS
{
	using namespace dealii;

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
 
    template void StationaryNavierStokes<ELEMENT_DIM>::refine_mesh() ;


}