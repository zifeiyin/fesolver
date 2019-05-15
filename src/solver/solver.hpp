#ifndef SOLVER_H
#define SOLVER_H
#include "../headers/top.hpp"
#include "../parameters/parameters.hpp"

namespace NavierStokes
{

template <int dim>
class IncompressibleNavierStokes
{
public:
  	IncompressibleNavierStokes( const RunTimeParameters::DataStorage& data );
   	void run() ;
protected:

private:


};


}

#endif