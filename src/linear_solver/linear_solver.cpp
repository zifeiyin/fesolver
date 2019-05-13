#include "../solver/solver.hpp"
#include "../preconditioner/preconditioner.hpp"

namespace IncompNS
{
using namespace dealii;

/**
 * solve linear equation
 */
template <int dim>
void StationaryNavierStokes<dim>::solve
(
    const bool initial_step
)
{
    const AffineConstraints<double>& 
        constraints_used = initial_step ? nonzero_constraints : zero_constraints;
    
    SolverControl solver_control(
        system_matrix.m(),          1e-4 * system_rhs.l2_norm(),        true    );
    
    SolverFGMRES<BlockVector<double>> gmres(    solver_control  );
    
    SparseILU<double> pmass_preconditioner;
    
    pmass_preconditioner.initialize(    
        pressure_mass_matrix,       SparseILU<double>::AdditionalData()         );
    
    const BlockSchurPreconditioner<SparseILU<double>> preconditioner(
        gamma,                      viscosity,                  system_matrix,
        pressure_mass_matrix,       pmass_preconditioner                        );

    gmres.solve(
        system_matrix,              newton_update,              system_rhs, 
        preconditioner                                                          ) ;

    std::cout << "FGMRES steps: " << solver_control.last_step() << std::endl;

    constraints_used.distribute( newton_update );
  }
  template void StationaryNavierStokes<ELEMENT_DIM>::solve(const bool initial_step) ;


/**
 * newton iteration
 */
template <int dim>
void StationaryNavierStokes<dim>::newton_iteration
(
    const double        tolerance,
    const unsigned int  max_n_line_searches,
    const unsigned int  max_n_refinements,
    const bool          is_initial_step,
    const bool          output_result
)
{
    bool first_step = is_initial_step;
    for ( unsigned int refinement_n = 0; refinement_n < max_n_refinements + 1; 
          ++refinement_n )
    {
        unsigned int line_search_n = 0;
        double       last_res      = 1.0;
        double       current_res   = 1.0;
        std::cout << "grid refinements: " << refinement_n << std::endl
                  << "viscosity: " << viscosity << std::endl;
        while ((first_step || (current_res > tolerance)) &&
               line_search_n < max_n_line_searches)
        {
            if (first_step)
            {
                setup_dofs();
                initialize_system();
                evaluation_point = present_solution;
                assemble_system(first_step);
                solve(first_step);
                present_solution = newton_update;
                nonzero_constraints.distribute(present_solution);
                first_step       = false;
                evaluation_point = present_solution;
                assemble_rhs(first_step);
                current_res = system_rhs.l2_norm();
                std::cout << "The residual of initial guess is " << current_res
                          << std::endl;
                last_res = current_res;
            }
            else
            {
                evaluation_point = present_solution;
                assemble_system(first_step);
                solve(first_step);
                for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5)
                {
                    evaluation_point = present_solution;
                    evaluation_point.add(alpha, newton_update);
                    nonzero_constraints.distribute(evaluation_point);
                    assemble_rhs(first_step);
                    current_res = system_rhs.l2_norm();
                    std::cout << "  alpha: " << std::setw(10) << alpha
                              << std::setw(0) << "  residual: " << current_res
                              << std::endl;
                    if (current_res < last_res)
                      break;
                }

                {
                  present_solution = evaluation_point;
                  std::cout << "  number of line searches: " << line_search_n
                            << "  residual: " << current_res << std::endl;
                  last_res = current_res;
                }

                ++line_search_n;
            }
            if (output_result)
            {
                output_results( max_n_line_searches * refinement_n + line_search_n );
                if (current_res <= tolerance)
                {
                  process_solution(refinement_n);
                }
            }
        }
        if (refinement_n < max_n_refinements)
        {
            refine_mesh();
        }
    }
}

template void StationaryNavierStokes<ELEMENT_DIM>::newton_iteration
(
    const double        tolerance,
    const unsigned int  max_n_line_searches,
    const unsigned int  max_n_refinements,
    const bool          is_initial_step,
    const bool          output_result
) ;


}