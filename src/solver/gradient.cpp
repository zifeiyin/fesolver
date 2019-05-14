#include "../solver/solver.hpp"

namespace NavierStokes
{

template <int dim>
void NavierStokesProjection<dim>::assemble_one_cell_of_gradient (
    const IteratorPair  &SI,
    InitGradScratchData &scratch,
    InitGradPerTaskData &data
)
{
    scratch.fe_val_vel.reinit (std::get<0> (*SI));
    scratch.fe_val_pres.reinit (std::get<1> (*SI));
    std::get<0> (*SI)->get_dof_indices (data.vel_local_dof_indices);
    std::get<1> (*SI)->get_dof_indices (data.pres_local_dof_indices);
    data.local_grad = 0.;
    for (unsigned int q=0; q<scratch.nqp; ++q)
      {
        for (unsigned int i=0; i<data.vel_dpc; ++i)
          for (unsigned int j=0; j<data.pres_dpc; ++j)
            data.local_grad (i, j) += -scratch.fe_val_vel.JxW(q) *
                                      scratch.fe_val_vel.shape_grad (i, q)[data.d] *
                                      scratch.fe_val_pres.shape_value (j, q);
      }
}
template void NavierStokesProjection<ELEMENT_DIM>::assemble_one_cell_of_gradient (
    const IteratorPair  &SI,
    InitGradScratchData &scratch,
    InitGradPerTaskData &data
);




template <int dim>
void NavierStokesProjection<dim>::copy_gradient_local_to_global(
    const InitGradPerTaskData &data
)
{
    for (unsigned int i=0; i<data.vel_dpc; ++i)
      for (unsigned int j=0; j<data.pres_dpc; ++j)
        pres_Diff[data.d].add (data.vel_local_dof_indices[i], data.pres_local_dof_indices[j],
                               data.local_grad (i, j) );
}
template void NavierStokesProjection<ELEMENT_DIM>::copy_gradient_local_to_global(
    const InitGradPerTaskData &data
) ;



}