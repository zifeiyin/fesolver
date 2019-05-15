#include "../solver/solver.hpp"

/**
 * main function to solve incompressible Navier-Stokes equation
 * Author: Dr. Zifei Yin at SJTU
 * Email:  zifeiyin@outlook.com
 */
int main()
{
  	try
    {
      	using namespace NavierStokes;

		RunTimeParameters::DataStorage data ;

		data.read_data( "problem.inp" ) ;

      	IncompressibleNavierStokes<ELEMENT_DIM> problem( data ) ;
      	problem.run ();
    }
  	catch (std::exception &exc)
    {
      	std::cerr 	<< std::endl << std::endl
                	<< "----------------------------------------------------"
                	<< std::endl;
      	std::cerr 	<< "Exception on processing: " << std::endl
                	<< exc.what() << std::endl
                	<< "Aborting!" << std::endl
                	<< "----------------------------------------------------"
                	<< std::endl;
      	return 1;
    }
  	catch (...)
    {
      	std::cerr 	<< std::endl << std::endl
                	<< "----------------------------------------------------"
                	<< std::endl;
      	std::cerr 	<< "Unknown exception!" << std::endl
                	<< "Aborting!" << std::endl
                	<< "----------------------------------------------------"
                	<< std::endl;
      	return 1;
    }
  	std::cout 	<< "----------------------------------------------------"
            	<< std::endl
            	<< "Apparently everything went fine!"
            	<< std::endl
            	<< "Don't forget to brush your teeth :-)"
            	<< std::endl << std::endl;
  	return 0;
}
