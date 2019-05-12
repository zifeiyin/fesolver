#ifndef SOLVER_H
#define SOLVER_H
#include "../headers/top.hpp"

namespace IncompNS
{
	using namespace dealii;

  	template <int dim>
  	class StationaryNavierStokes
  	{
  	public:
    	StationaryNavierStokes(     const unsigned int  degree              );
    	void run(                   const unsigned int  refinement          );
  	private:
    	void setup_dofs();
    	void initialize_system();
    	void assemble(              const bool          initial_step, 
                                    const bool          assemble_matrix     );
    	void assemble_system(       const bool          initial_step        );
    	void assemble_rhs(          const bool          initial_step        );
    	void solve(                 const bool          initial_step        );
    	void refine_mesh();
    	void process_solution(      unsigned int        refinement          );
    	void output_results(        const unsigned int  refinement_cycle    ) const;
    	void newton_iteration(	    const double        tolerance,
                          		    const unsigned int  max_n_line_searches,
                          		    const unsigned int  max_n_refinements,
                          		    const bool          is_initial_step,
                          		    const bool          output_result       );
    	void compute_initial_guess( double              step_size           );

    	double                                          viscosity;
    	double                                          gamma;
    	const unsigned int                              degree;
    	std::vector<types::global_dof_index>            dofs_per_block;
    	Triangulation<dim>                              triangulation;
    	FESystem<dim>                                   fe;
    	DoFHandler<dim>                                 dof_handler;
    	AffineConstraints<double>                       zero_constraints;
    	AffineConstraints<double>                       nonzero_constraints;
    	BlockSparsityPattern                            sparsity_pattern;
    	BlockSparseMatrix<double>                       system_matrix;
    	SparseMatrix<double>                            pressure_mass_matrix;
    	BlockVector<double>                             present_solution;
    	BlockVector<double>                             newton_update;
    	BlockVector<double>                             system_rhs;
    	BlockVector<double>                             evaluation_point;
  	};

}

#endif