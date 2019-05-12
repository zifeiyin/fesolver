#include "../headers/top.hpp"

namespace IncompNS
{
    using namespace dealii ;

    template <class PreconditionerMp>
  	class BlockSchurPreconditioner : public Subscriptor
  	{
  	public:
    	BlockSchurPreconditioner(double                           gamma,
    	                         double                           viscosity,
    	                         const BlockSparseMatrix<double> &S,
    	                         const SparseMatrix<double> &     P,
    	                         const PreconditionerMp &         Mppreconditioner);
    	void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const;
  	private:
    	const double                     gamma;
    	const double                     viscosity;
    	const BlockSparseMatrix<double> &stokes_matrix;
    	const SparseMatrix<double> &     pressure_mass_matrix;
    	const PreconditionerMp &         mp_preconditioner;
    	SparseDirectUMFPACK              A_inverse;
  	};
  

      
}