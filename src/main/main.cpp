#include "../solver/solver.hpp"

int main()
{
  try
    {
      using namespace dealii;
      using namespace NavierStokes;
      RunTimeParameters::Data_Storage data;
      data.read_data ("parameter-file.prm");
      deallog.depth_console (data.verbose ? 2 : 0);
      NavierStokesProjection<2> test (data);
      test.run (data.verbose, data.output_interval);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  std::cout << "----------------------------------------------------"
            << std::endl
            << "Apparently everything went fine!"
            << std::endl
            << "Don't forget to brush your teeth :-)"
            << std::endl << std::endl;
  return 0;
}
