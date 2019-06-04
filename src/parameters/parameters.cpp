#include "parameters.hpp"

namespace NavierStokes
{
namespace RunTimeParameters
{

/**
 * return file name as constant
 */
const string DataStorage::mesh_file_name() const {
    return this->mesh_file ;
}





}
}