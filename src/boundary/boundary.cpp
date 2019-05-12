#include "boundary.hpp"

namespace IncompNS
{
	using namespace dealii;

	template <>
  	double BoundaryValues<ELEMENT_DIM>::value	
	(
		const Point<ELEMENT_DIM> & p,
        const unsigned int component	
	) const
  	{
    	Assert(component < this->n_components,
           		ExcIndexRange(component, 0, this->n_components));

    	if (component == 0 && std::abs(p[ELEMENT_DIM - 1] - 1.0) < 1e-10) {
      		return 1.0;
		} else {
    		return 0;
		}
  	}


}