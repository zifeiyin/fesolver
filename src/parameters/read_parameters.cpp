#include "parameters.hpp"

namespace NavierStokes
{
namespace RunTimeParameters
{
/**
 * read input file which is the problem.input file in current folder
 */
void DataStorage::read_data(
    const string input_file_name
){
    std::ifstream file ( input_file_name ) ;

    AssertThrow( file, ExcFileNotOpen( input_file_name ) ) ;

    prmHd.parse_input( file ) ;

    prmHd.enter_subsection( "Mesh" ) ;
    {
        mesh_file       = prmHd.get( "mesh file"        ) ;
    }
    prmHd.leave_subsection() ;

    prmHd.enter_subsection( "Physical data" ) ;
    {
        dt              = prmHd.get_double( "dt"            ) ;
        initial_time    = prmHd.get_double( "initial time"  ) ;
        final_time      = prmHd.get_double( "final time"    ) ;
        nu              = prmHd.get_double( "viscosity"     ) ;
    }
    prmHd.leave_subsection() ;

    prmHd.enter_subsection( "Solver control" ) ;
    {
        convergence_limit = prmHd.get_double( "convergence limit" ) ;
    }
    prmHd.leave_subsection() ;

    print_settings() ;
}

/**
 * print out data that read from input file
 */
void DataStorage::print_settings()
{
    cout << "printing input data"                           << endl ;
    cout << "viscosity = "          << nu                   << endl ;
    cout << "dt = "                 << dt                   << endl ;
    cout << "initial time = "       << initial_time         << endl ;
    cout << "final time = "         << final_time             << endl ;
    cout << "convergence limit = "  << convergence_limit    << endl ;
}


}
}