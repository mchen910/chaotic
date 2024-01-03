#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <functional>
#include <vector>

#include "solver.h"


namespace DES 
{

/**
 * @brief Solve a singular first order ordinary differential equation of the form 
 * dy/dt = f(t, y, ...). If the ODE is second order or above, use the ODESystem class 
 * and decompose the differential equation into multiple first order equations before 
 * proceeding.
 * 
 */
template <typename T, typename V>
class ODE : public DiffEq<T, V>
{

public:
    ODE(function_t<T, V>& func, timeBound_t& bounds, T initialCondition, double timeStep) 
    : DiffEq<T, V>(func, bounds, initialCondition, timeStep) { }

};


/**
 * @brief Solve a system of first order ordinary differential equations, all of the 
 * form dy_j/dt = f_j(t, y_1, y_2, ...). 
 * 
 * @tparam T 
 * @tparam V 
 */
template <typename T, typename V>
class ODESystem : public System<T, V>
{

private:
    timeBound_t _timeBound;
    double _timeStep;
    size_t _equations;

public:
    ODESystem(iv_t& iValues, timeBound_t& bounds, double timeStep, std::initializer_list<function_t<T, V>> funcs)
    : System(iValues, funcs) 
    { 
        _timeBound = bounds;
        _timeStep = timeStep;
        _equations = _functions.size();
    };


};


template <typename T, typename V>
struct Algorithms
{
    DataFrame<T>  _euler  (ODESystem<T, V>& ode);
    DataFrame<T>  _RK4    (ODESystem<T, V>& ode);
    DataFrame<T>  _RK45   (ODESystem<T, V>& ode);
    DataFrame<T>  _TSIT5  (ODESystem<T, V>& ode);
};

}


#endif