#include "eqnData.hpp"

namespace NavierStokes
{
namespace EquationData
{

using namespace dealii ;

template <int dim>
MultiComponentFunction<dim>::MultiComponentFunction ( const double initial_time )
:
    Function<dim> (1, initial_time), comp(0)
{}
template MultiComponentFunction<ELEMENT_DIM>::MultiComponentFunction ( const double initial_time ) ;



template <int dim>
void MultiComponentFunction<dim>::set_component( const unsigned int d )
{
    Assert (d<dim, ExcIndexRange (d, 0, dim));
    comp = d;
}
template void MultiComponentFunction<ELEMENT_DIM>::set_component( const unsigned int d ) ;



template <int dim>
Velocity<dim>::Velocity (const double initial_time)
:
MultiComponentFunction<dim> (initial_time)
{}
template Velocity<ELEMENT_DIM>::Velocity (const double initial_time) ;

    
template <int dim>
void Velocity<dim>::value_list (
    const std::vector<Point<dim> > &points,
    std::vector<double> &values,
    const unsigned int ) const
{
    const unsigned int n_points = points.size();
    Assert (    values.size() == n_points,
                ExcDimensionMismatch (values.size(), n_points));
    for (unsigned int i=0; i<n_points; ++i){
        values[i] = Velocity<dim>::value (points[i]);
    }
}
template void Velocity<ELEMENT_DIM>::value_list (
    const std::vector<Point<ELEMENT_DIM> > &points,
    std::vector<double> &values,
    const unsigned int 
) const ;



template <int dim>
double Velocity<dim>::value (   
    const Point<dim> &p,
    const unsigned int 
) const
{
    if (this->comp == 0)
    {
        const double Um = 1.5;
        const double H  = 4.1;
        return 4.*Um*p(1)*(H - p(1))/(H*H);
    }
    else{
        return 0.;
    }
}
template double Velocity<ELEMENT_DIM>::value (   
    const Point<ELEMENT_DIM> &p,
    const unsigned int 
) const ;



template <int dim>
Pressure<dim>::Pressure (   const double initial_time )
:
Function<dim> (1, initial_time)
{}
template Pressure<ELEMENT_DIM>::Pressure (   const double initial_time ) ;


template <int dim>
double Pressure<dim>::value (   
    const Point<dim> &p,
    const unsigned int 
) const
{
    return 25.-p(0);
}
template double Pressure<ELEMENT_DIM>::value (   
    const Point<ELEMENT_DIM> &p,
    const unsigned int 
) const ;



template <int dim>
void Pressure<dim>::value_list (    
    const std::vector<Point<dim> > &points,
    std::vector<double> &values,
    const unsigned int
) const
{
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points, ExcDimensionMismatch (values.size(), n_points));
    for (unsigned int i=0; i<n_points; ++i){
        values[i] = Pressure<dim>::value (points[i]);
    }
}
template void Pressure<ELEMENT_DIM>::value_list (    
    const std::vector<Point<ELEMENT_DIM> > &points,
    std::vector<double> &values,
    const unsigned int
) const ;


}
}