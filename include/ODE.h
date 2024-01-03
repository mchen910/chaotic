#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <functional>
#include <vector>

#include "solver.h"


namespace DES 
{


/**
 * @brief Solve a 
 * 
 */
template <typename T, typename V>
class ODE : public DiffEq<T, V>
{

public:
    ODE(function_t<T, V> func, timeBound_t t, T initialCondition) 
    : DiffEq(func, initialCondition, t) 
    {

    }

};






}


#endif