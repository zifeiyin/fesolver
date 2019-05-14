#ifndef SOLVER_H
#define SOLVER_H
#include "../headers/top.hpp"
#include "../prm/prm.hpp"
#include "../eqn/eqnData.hpp"

namespace NavierStokes
{

  template <int dim>
  class NavierStokesProjection
  {
  public:
    NavierStokesProjection (const RunTimeParameters::Data_Storage &data);
    void run (const bool         verbose    = false,
              const unsigned int n_plots = 10);
  protected:
    RunTimeParameters::MethodFormulation        type;
    const unsigned int                          deg;
    const double                                dt;
    const double                                t_0, T, Re;
    EquationData::Velocity<dim>                 vel_exact;
    std::map<types::global_dof_index, double>   boundary_values;
    std::vector<types::boundary_id>             boundary_ids;
    Triangulation<dim>                          triangulation;
    FE_Q<dim>                                   fe_velocity;
    FE_Q<dim>                                   fe_pressure;
    DoFHandler<dim>                             dof_handler_velocity;
    DoFHandler<dim>                             dof_handler_pressure;
    QGauss<dim>                                 quadrature_pressure;
    QGauss<dim>                                 quadrature_velocity;
    SparsityPattern                             sparsity_pattern_velocity;
    SparsityPattern                             sparsity_pattern_pressure;
    SparsityPattern                             sparsity_pattern_pres_vel;
    SparseMatrix<double>                        vel_Laplace_plus_Mass;
    SparseMatrix<double>                        vel_it_matrix[dim];
    SparseMatrix<double>                        vel_Mass;
    SparseMatrix<double>                        vel_Laplace;
    SparseMatrix<double>                        vel_Advection;
    SparseMatrix<double>                        pres_Laplace;
    SparseMatrix<double>                        pres_Mass;
    SparseMatrix<double>                        pres_Diff[dim];
    SparseMatrix<double>                        pres_iterative;

    Vector<double>                              pres_n;
    Vector<double>                              pres_n_minus_1;
    Vector<double>                              phi_n;
    Vector<double>                              phi_n_minus_1;
    Vector<double>                              u_n[dim];
    Vector<double>                              u_n_minus_1[dim];
    Vector<double>                              u_star[dim];
    Vector<double>                              force[dim];
    Vector<double>                              v_tmp;
    Vector<double>                              pres_tmp;
    Vector<double>                              rot_u;
    SparseILU<double>                           prec_velocity[dim];
    SparseILU<double>                           prec_pres_Laplace;
    SparseDirectUMFPACK                         prec_mass;
    SparseDirectUMFPACK                         prec_vel_mass;

    DeclException2 (ExcInvalidTimeStep,
                    double, double,
                    << " The time step " << arg1 << " is out of range."
                    << std::endl
                    << " The permitted range is (0," << arg2 << "]");

    void create_triangulation_and_dofs (const unsigned int n_refines);
    void initialize();
    void interpolate_velocity ();
    void diffusion_step (const bool reinit_prec);
    void projection_step (const bool reinit_prec);
    void update_pressure (const bool reinit_prec);


  private:
    unsigned int                                vel_max_its;
    unsigned int                                vel_Krylov_size;
    unsigned int                                vel_off_diagonals;
    unsigned int                                vel_update_prec;
    double                                      vel_eps;
    double                                      vel_diag_strength;
    void initialize_velocity_matrices();
    void initialize_pressure_matrices();
    typedef std::tuple< typename DoFHandler<dim>::active_cell_iterator,
            typename DoFHandler<dim>::active_cell_iterator
            > IteratorTuple;
    typedef SynchronousIterators<IteratorTuple> IteratorPair;
    void initialize_gradient_operator();


    #include "../headers/scratch_structs.hpp"

    void assemble_one_cell_of_gradient (const IteratorPair  &SI,
                                        InitGradScratchData &scratch,
                                        InitGradPerTaskData &data);
    void copy_gradient_local_to_global (const InitGradPerTaskData &data);
    void assemble_advection_term();



    void assemble_one_cell_of_advection (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                         AdvectionScratchData &scratch,
                                         AdvectionPerTaskData &data);
    void copy_advection_local_to_global (const AdvectionPerTaskData &data);
    void diffusion_component_solve (const unsigned int d);
    void output_results (const unsigned int step);
    void assemble_vorticity (const bool reinit_prec);
  };


}

#endif