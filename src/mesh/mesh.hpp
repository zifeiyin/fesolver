#ifndef MESH_H
#define MESH_H
#include "../headers/top.hpp"
#include "../parameters/parameters.hpp"

namespace NavierStokes
{
/**
 * storage of mesh data structure from reading 
 */
template<int dim>
class MeshStorage
{
public:
    MeshStorage( const RunTimeParameters::DataStorage&   data ) ;
    GridIn<dim>         grid_in ;
    Triangulation<dim>   triangulation;
protected:
private:

} ;




}

#endif