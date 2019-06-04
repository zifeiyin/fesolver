#include "mesh.hpp"

namespace NavierStokes
{

/**
 * read mesh based on directory provided by input file
 */
template<int dim>
MeshStorage<dim>::MeshStorage(
    const RunTimeParameters::DataStorage&   input_data
){

    grid_in.attach_triangulation (triangulation);
    
    string read_name = input_data.mesh_file_name() ;
    
    istringstream str_mesh( read_name ) ;
    
    grid_in.read_ucd( str_mesh );

}
template MeshStorage<ELEMENT_DIM>::MeshStorage(
   const RunTimeParameters::DataStorage&   input_data 
);


}