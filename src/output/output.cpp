#include "../solver/solver.hpp"

namespace IncompNS
{
	using namespace dealii;

  template <>
  void StationaryNavierStokes<ELEMENT_DIM>::process_solution(unsigned int refinement)
  {
    std::ostringstream filename;
    filename << (1.0/viscosity) << "-line-" << refinement << ".txt";
    std::ofstream f (filename.str().c_str());
    f << "# y u_x u_y" << std::endl;
    Point<ELEMENT_DIM> p;
    p(0)= 0.5;
    p(1)= 0.5;
    f << std::scientific;
    for (unsigned int i=0; i<=100; ++i)
      {
        p(ELEMENT_DIM-1) = i/100.0;
        Vector<double> tmp_vector(ELEMENT_DIM+1);
        VectorTools::point_value(dof_handler, present_solution, p, tmp_vector);
        f << p(ELEMENT_DIM-1);
        for (int j=0; j<ELEMENT_DIM; j++)
          f << " " << tmp_vector(j);
        f << std::endl;
      }
  }

  template <>
  void StationaryNavierStokes<ELEMENT_DIM>::output_results(
    const unsigned int output_index) const
  {
    std::vector<std::string> solution_names(ELEMENT_DIM, "velocity");
    solution_names.emplace_back("pressure");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        ELEMENT_DIM, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    DataOut<ELEMENT_DIM> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(present_solution,
                             solution_names,
                             DataOut<ELEMENT_DIM>::type_dof_data,
                             data_component_interpretation);
    data_out.build_patches();
    std::ofstream output(std::to_string(1.0 / viscosity) + "-solution-" +
                         Utilities::int_to_string(output_index, 4) + ".vtk");
    data_out.write_vtk(output);
  }


}