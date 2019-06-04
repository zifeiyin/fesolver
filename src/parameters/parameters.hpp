#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "../headers/top.hpp"

namespace NavierStokes
{
namespace RunTimeParameters
{
/**
 * storage of data structure from input file 
 * prmHd: parameter Handler from dealii library
 */
class DataStorage
{
public:
    void                read_data(  const string file_name ) ;
    DataStorage() ;

    const string        mesh_file_name() const ;

protected:
    ParameterHandler    prmHd ;

private:
    string              mesh_file ;

    double              dt ;
    double              initial_time ;
    double              final_time ;

    double              nu ;

    double              convergence_limit ;
    
    void                print_settings() ;

} ;




}
}

#endif