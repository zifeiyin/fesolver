/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

void first_grid()
{
    const int dim = 2 ;
    Triangulation<2> triangulation;
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(5);

    typename Triangulation<dim>::cell_iterator 
    cell = triangulation.begin() , 
    endc = triangulation.end()   ;
    for ( ; cell != endc ; ++cell ){
        for ( unsigned int face_number = 0 ; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number ) 
        {
            if ((std::fabs(cell->face(face_number)->center()(1) - (1.0)) < 1e-12) ) {
                cell->face(face_number)->set_boundary_id (1);
            } else if ( (std::fabs(cell->face(face_number)->center()(1) - (0.0)) < 1e-12) ){
                cell->face(face_number)->set_boundary_id (2);
            } else if ( (std::fabs(cell->face(face_number)->center()(0) - (0.0)) < 1e-12) ){
                cell->face(face_number)->set_boundary_id (2);
            } else if ( (std::fabs(cell->face(face_number)->center()(0) - (1.0)) < 1e-12) ){
                cell->face(face_number)->set_boundary_id (2);
            }
        }
    }

    std::ofstream out("block.ucd");
    GridOut       grid_out;
    GridOutFlags::Ucd ucd_flags(false, true, false) ;

    grid_out.set_flags( ucd_flags ) ;
    grid_out.write_ucd(triangulation, out);
    std::cout << "Grid written to block.ucd" << std::endl;

}

int main()
{
  first_grid();
}
