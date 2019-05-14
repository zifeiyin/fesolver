#include "../solver/solver.hpp"

namespace NavierStokes
{

template <int dim>
void NavierStokesProjection<dim>::projection_step (const bool reinit_prec)
{
    pres_iterative.copy_from (pres_Laplace);
    pres_tmp = 0.;
    for ( unsigned int d = 0 ; d < dim ; ++ d ){
        pres_Diff[d].Tvmult_add (pres_tmp, u_n[d]);
    }
    phi_n_minus_1 = phi_n;
    
    static std::map<types::global_dof_index, double> bval;
    if (reinit_prec){
      VectorTools::interpolate_boundary_values (dof_handler_pressure, 3,
                                                Functions::ZeroFunction<dim>(), bval);
    }
    MatrixTools::apply_boundary_values (bval, pres_iterative, phi_n, pres_tmp);
    
    if (reinit_prec){
        prec_pres_Laplace.initialize(
            pres_iterative,
            SparseILU<double>::AdditionalData (vel_diag_strength, vel_off_diagonals) 
        );
    }

    SolverControl solvercontrol (vel_max_its, vel_eps*pres_tmp.l2_norm());
    
    SolverCG<> cg (solvercontrol);
    
    cg.solve (pres_iterative, phi_n, pres_tmp, prec_pres_Laplace);
    
    phi_n *= 1.5/dt;
}
template void NavierStokesProjection<ELEMENT_DIM>::projection_step (const bool reinit_prec) ;



}