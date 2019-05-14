#include "solver.hpp"

namespace NavierStokes
{

using namespace dealii;

template <int dim>
void NavierStokesProjection<dim>::run (
    const bool verbose,
    const unsigned int output_interval
)
{
    ConditionalOStream verbose_cout (std::cout, verbose);
    const unsigned int n_steps =  static_cast<unsigned int>((T - t_0)/dt);
    vel_exact.set_time (2.*dt);
    output_results(1);
    for (unsigned int n = 2; n<=n_steps; ++n)
      {
        if (n % output_interval == 0)
          {
            verbose_cout << "Plotting Solution" << std::endl;
            output_results(n);
          }
        std::cout << "Step = " << n << " Time = " << (n*dt) << std::endl;
        verbose_cout << "  Interpolating the velocity " << std::endl;
        interpolate_velocity();
        verbose_cout << "  Diffusion Step" << std::endl;
        if (n % vel_update_prec == 0)
          verbose_cout << "    With reinitialization of the preconditioner"
                       << std::endl;
        diffusion_step ((n%vel_update_prec == 0) || (n == 2));
        verbose_cout << "  Projection Step" << std::endl;
        projection_step ( (n == 2));
        verbose_cout << "  Updating the Pressure" << std::endl;
        update_pressure ( (n == 2));
        vel_exact.advance_time(dt);
      }
    output_results (n_steps);
}
template void NavierStokesProjection<ELEMENT_DIM>::run ( 
    const bool          verbose,
    const unsigned int  output_interval 
) ;


template <int dim>
NavierStokesProjection<dim>::NavierStokesProjection(const RunTimeParameters::Data_Storage &data)
:
    type (data.form),
    deg (data.pressure_degree),
    dt (data.dt),
    t_0 (data.initial_time),
    T (data.final_time),
    Re (data.Reynolds),
    vel_exact (data.initial_time),
    fe_velocity (deg+1),
    fe_pressure (deg),
    dof_handler_velocity (triangulation),
    dof_handler_pressure (triangulation),
    quadrature_pressure (deg+1),
    quadrature_velocity (deg+2),
    vel_max_its (data.vel_max_iterations),
    vel_Krylov_size (data.vel_Krylov_size),
    vel_off_diagonals (data.vel_off_diagonals),
    vel_update_prec (data.vel_update_prec),
    vel_eps (data.vel_eps),
    vel_diag_strength (data.vel_diag_strength)
{
    if (deg < 1)
      std::cout << " WARNING: The chosen pair of finite element spaces is not stable."
                << std::endl
                << " The obtained results will be nonsense"
                << std::endl;
    AssertThrow (!  ( (dt <= 0.) || (dt > .5*T)), ExcInvalidTimeStep (dt, .5*T));
    create_triangulation_and_dofs (data.n_global_refines);
    initialize();
}
template NavierStokesProjection<ELEMENT_DIM>::NavierStokesProjection(const RunTimeParameters::Data_Storage &data) ;




template <int dim>
void NavierStokesProjection<dim>::create_triangulation_and_dofs (
    const unsigned int n_refines
)
{
    GridIn<dim> grid_in;
    grid_in.attach_triangulation (triangulation);
    {
      std::string filename = "nsbench2.inp";
      std::ifstream file (filename.c_str());
      Assert (file, ExcFileNotOpen (filename.c_str()));
      grid_in.read_ucd (file);
    }
    std::cout << "Number of refines = " << n_refines
              << std::endl;
    triangulation.refine_global (n_refines);
    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;
    boundary_ids = triangulation.get_boundary_ids();
    dof_handler_velocity.distribute_dofs (fe_velocity);
    DoFRenumbering::boost::Cuthill_McKee (dof_handler_velocity);
    dof_handler_pressure.distribute_dofs (fe_pressure);
    DoFRenumbering::boost::Cuthill_McKee (dof_handler_pressure);
    initialize_velocity_matrices();
    initialize_pressure_matrices();
    initialize_gradient_operator();
    pres_n.reinit (dof_handler_pressure.n_dofs());
    pres_n_minus_1.reinit (dof_handler_pressure.n_dofs());
    phi_n.reinit (dof_handler_pressure.n_dofs());
    phi_n_minus_1.reinit (dof_handler_pressure.n_dofs());
    pres_tmp.reinit (dof_handler_pressure.n_dofs());
    for (unsigned int d=0; d<dim; ++d)
      {
        u_n[d].reinit (dof_handler_velocity.n_dofs());
        u_n_minus_1[d].reinit (dof_handler_velocity.n_dofs());
        u_star[d].reinit (dof_handler_velocity.n_dofs());
        force[d].reinit (dof_handler_velocity.n_dofs());
      }
    v_tmp.reinit (dof_handler_velocity.n_dofs());
    rot_u.reinit (dof_handler_velocity.n_dofs());
    std::cout << "dim (X_h) = " << (dof_handler_velocity.n_dofs()*dim)
              << std::endl
              << "dim (M_h) = " << dof_handler_pressure.n_dofs()
              << std::endl
              << "Re        = " << Re
              << std::endl
              << std::endl;
}
template void NavierStokesProjection<ELEMENT_DIM>::create_triangulation_and_dofs ( const unsigned int n_refines) ;



template <int dim>
void NavierStokesProjection<dim>::interpolate_velocity()
{
    for (unsigned int d=0; d<dim; ++d)
    {
        u_star[d].equ (2., u_n[d]);
        u_star[d] -=  u_n_minus_1[d];
    }
}
template void NavierStokesProjection<ELEMENT_DIM>::interpolate_velocity() ;



template <int dim>
void NavierStokesProjection<dim>::update_pressure (const bool reinit_prec)
{
    pres_n_minus_1 = pres_n;
    switch (type)
      {
      case RunTimeParameters::METHOD_STANDARD:
        pres_n += phi_n;
        break;
      case RunTimeParameters::METHOD_ROTATIONAL:
        if (reinit_prec)
          prec_mass.initialize (pres_Mass);
        pres_n = pres_tmp;
        prec_mass.solve (pres_n);
        pres_n.sadd(1./Re, 1., pres_n_minus_1);
        pres_n += phi_n;
        break;
      default:
        Assert (false, ExcNotImplemented());
      };
}
template void NavierStokesProjection<ELEMENT_DIM>::update_pressure (const bool reinit_prec);



template <int dim>
void NavierStokesProjection<dim>::assemble_vorticity (const bool reinit_prec)
{
    Assert (dim == 2, ExcNotImplemented());
    if (reinit_prec)
      prec_vel_mass.initialize (vel_Mass);
    FEValues<dim> fe_val_vel (fe_velocity, quadrature_velocity,
                              update_gradients |
                              update_JxW_values |
                              update_values);
    const unsigned int dpc = fe_velocity.dofs_per_cell,
                       nqp = quadrature_velocity.size();
    std::vector<types::global_dof_index> ldi (dpc);
    Vector<double> loc_rot (dpc);
    std::vector< Tensor<1,dim> > grad_u1 (nqp), grad_u2 (nqp);
    rot_u = 0.;
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler_velocity.begin_active(),
    end  = dof_handler_velocity.end();
    for (; cell != end; ++cell)
      {
        fe_val_vel.reinit (cell);
        cell->get_dof_indices (ldi);
        fe_val_vel.get_function_gradients (u_n[0], grad_u1);
        fe_val_vel.get_function_gradients (u_n[1], grad_u2);
        loc_rot = 0.;
        for (unsigned int q=0; q<nqp; ++q)
          for (unsigned int i=0; i<dpc; ++i)
            loc_rot(i) += (grad_u2[q][0] - grad_u1[q][1]) *
                          fe_val_vel.shape_value (i, q) *
                          fe_val_vel.JxW(q);
        for (unsigned int i=0; i<dpc; ++i)
          rot_u (ldi[i]) += loc_rot(i);
      }
    prec_vel_mass.solve (rot_u);
}
template void NavierStokesProjection<ELEMENT_DIM>::assemble_vorticity (const bool reinit_prec) ;



}