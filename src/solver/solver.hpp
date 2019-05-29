#ifndef SOLVER_H
#define SOLVER_H
#include "../headers/top.hpp"
#include "../parameters/parameters.hpp"

namespace NavierStokes
{

/**
 * Incompressible Navier-Stokes solver using Stablized Finite Element method,
 * the Steamline-upwind Petrov-Galerkin (SUPG) method or the Galerkin/Least-squares 
 * method. Since it is Continuous Galerkin method, it only applys to incompressible
 * flow or weakly compressible flow
 */
template <int dim>
class IncompressibleNavierStokes
{
public:
  	IncompressibleNavierStokes( const RunTimeParameters::DataStorage& data );
   	void run() ;
protected:

private:


};

/**
 * Compressible Navier-Stokes solver which uses Discontinuous Galerkin method.
 * Awaiting development
 */





}

#endif