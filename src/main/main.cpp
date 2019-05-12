#include "../solver/solver.hpp"

int main()
{
  	try
    {
      	using namespace dealii;
      	using namespace IncompNS;
      	StationaryNavierStokes<ELEMENT_DIM> flow(2); //degree
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