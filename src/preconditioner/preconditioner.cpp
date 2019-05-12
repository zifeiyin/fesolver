#include "preconditioner.hpp"

namespace IncompNS
{
	using namespace dealii;
  
  	template <>
  	BlockSchurPreconditioner<SparseILU<double>>::BlockSchurPreconditioner(
    	double                           gamma,
    	double                           viscosity,
    	const BlockSparseMatrix<double> &S,
    	const SparseMatrix<double> &     P,
    	const SparseILU<double> &         Mppreconditioner)
    	: gamma(gamma)
    	, viscosity(viscosity)
    	, stokes_matrix(S)
    	, pressure_mass_matrix(P)
    	, mp_preconditioner(Mppreconditioner)
  	{	
    	A_inverse.initialize(stokes_matrix.block(0, 0));
  	}



  	template <>
  	void BlockSchurPreconditioner<SparseILU<double>>::vmult(
    	BlockVector<double> &      dst,
    	const BlockVector<double> &src) const
  	{
    	Vector<double> utmp(src.block(0));
    	{
      	SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
      	SolverCG<>    cg(solver_control);
      	dst.block(1) = 0.0;
      	cg.solve(pressure_mass_matrix,
      	         dst.block(1),
      	         src.block(1),
      	         mp_preconditioner);
      	dst.block(1) *= -(viscosity + gamma);
    	}
    	{
      		stokes_matrix.block(0, 1).vmult(utmp, dst.block(1));
      		utmp *= -1.0;
      		utmp += src.block(0);
    	}
    	A_inverse.vmult(dst.block(0), utmp);
  	}



}