#include "solver.hpp"

int main()
{
  	try
    {
      	using namespace dealii;
      	using namespace IncompNS;
      	StationaryNavierStokes<2> flow(1); //degree
      	flow.run(4);
    }
  	
	catch (std::exception &exc)
    {
      	std::cerr 	<< std::endl
                	<< std::endl
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
      	std::cerr 	<< std::endl
                	<< std::endl
                	<< "----------------------------------------------------"
                	<< std::endl;
      	std::cerr 	<< "Unknown exception!" << std::endl
                	<< "Aborting!" << std::endl
                	<< "----------------------------------------------------"
                	<< std::endl;
      	return 1;
    }
  	return 0;
}