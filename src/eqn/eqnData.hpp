#ifndef EQUATION_DATA_H
#define EQUATION_DATA_H

#include "../headers/top.hpp"

namespace NavierStokes
{
namespace EquationData
{

using namespace dealii ;

template <int dim>
class MultiComponentFunction: public Function<dim>
{
public:
    MultiComponentFunction (const double initial_time = 0.);
    void set_component (const unsigned int d);
protected:
    unsigned int comp;
};


template <int dim>
class Velocity : public MultiComponentFunction<dim>
{
public:
    Velocity (const double initial_time = 0.0);
    virtual double value (      const Point<dim> &p,
                                const unsigned int component = 0) const;
    virtual void value_list(    const std::vector< Point<dim> > &points,
                                std::vector<double> &values,
                                const unsigned int component = 0) const;
};


template <int dim>
class Pressure: public Function<dim>
{
public:
    Pressure (const double initial_time = 0.0);
    virtual double value (      const Point<dim> &p,
                                const unsigned int component = 0) const;
    virtual void value_list (   const std::vector< Point<dim> > &points,
                                std::vector<double> &values,
                                const unsigned int component = 0) const;
};



}
}

#endif