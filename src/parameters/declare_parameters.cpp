#include "parameters.hpp"

namespace NavierStokes
{
namespace RunTimeParameters
{
/**
 * Default constructor of Data storage
 */
DataStorage::DataStorage() :
    dt                  (   0.01    ),
    initial_time        (   0.0     ),
    final_time          (   100.0   ),
    nu                  (   1.5e-5  ),
    convergence_limit   (   1.e-4   )
{
    prmHd.enter_subsection(     "Mesh"                                          ) ;
    {
        prmHd.declare_entry(    "mesh file",
                                "",
                                Patterns::Anything(),
                                "The mesh file name "                           ) ;
    }
    prmHd.leave_subsection();

    prmHd.enter_subsection(     "Physical data"                                 ) ;
    {
        prmHd.declare_entry (   "dt",   
                                "0.0",
                                Patterns::Double(0.0) ,
                                "The time step of the transient simulation"     ) ;

        prmHd.declare_entry (   "initial time",   
                                "0.0",
                                Patterns::Double(0.0) ,
                                "The initial time of the simulation"            ) ;

        prmHd.declare_entry (   "final time",   
                                "0.0",
                                Patterns::Double(0.0) ,
                                "The final time of the simulation"              ) ;
                            
        prmHd.declare_entry (   "viscosity",   
                                "0.0",
                                Patterns::Double(0.0) ,
                                "The kinematic viscosity of fluid"              ) ;
    }
    prmHd.leave_subsection() ;

    prmHd.enter_subsection(     "Solver control"                                ) ;
    {
        prmHd.declare_entry (   "convergence limit",   
                                "0.0001",
                                Patterns::Double(0.0) ,
                                "The convergence limit"                         ) ;
    }
    prmHd.leave_subsection() ;


}



}
}