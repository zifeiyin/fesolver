#ifndef PROBLEM_H
#define PROBLEM_H

#include "../headers/top.hpp"

namespace NavierStokes
{
namespace RunTimeParameters
{
    
enum MethodFormulation
{
    METHOD_STANDARD,
    METHOD_ROTATIONAL
};

class Data_Storage
{
public:
    Data_Storage();
    void read_data (const char *filename);
    MethodFormulation   form;
    double              initial_time,
                        final_time,
                        Reynolds;
    double              dt;
    unsigned int        n_global_refines,
                        pressure_degree;
    unsigned int        vel_max_iterations,
                        vel_Krylov_size,
                        vel_off_diagonals,
                        vel_update_prec;
    double              vel_eps,
                        vel_diag_strength;
    bool                verbose;
    unsigned int        output_interval;
protected:
    ParameterHandler    prm;
};


}

}


#endif