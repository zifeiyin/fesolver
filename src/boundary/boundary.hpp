#include "../headers/top.hpp"

namespace IncompNS
{
	using namespace dealii;

    template <int dim>
  	class BoundaryValues : public Function<dim>
  	{
  	public:
    	BoundaryValues()
    	  : Function<dim>(dim + 1)
    	{}
    	virtual double value(const Point<dim> & p,
    	                     const unsigned int component) const override;
  	};
  	
}