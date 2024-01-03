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
class ODESystem : public DiffEqSystem<T, V>
{

private:
    timeBound_t _timeBound;
    double _timeStep;
    size_t _equations;

public:
    ODESystem(iv_t& iValues, std::initializer_list<function_t<T, V>> funcs, timeBound_t& bounds, double timeStep)
    : DiffEqSystem<T, V>(iValues, funcs)
    { 
        _timeBound = bounds;
        _timeStep = timeStep;
        _equations = funcs.size();
    };

    std::vector<T> _eval(std::vector<V>& inputs); 
};


// Forward declare the algorithms

template <typename T, typename V>
DataFrame<T>  _euler  (ODESystem<T, V> const& ode);


template <typename T, typename V>
DataFrame<T>  _RK4    (ODESystem<T, V> const& ode);


template <typename T, typename V>
DataFrame<T>  _RK45   (ODESystem<T, V> const& ode);


template <typename T, typename V>
DataFrame<T>  _TSIT5  (ODESystem<T, V> const& ode);





/**
 * @brief Evaluate an entire system of functions with one set of 
 * inputs. Useful for coupled Runge-Kutta schemes, which use common sets  
 * of inputs for every iteration.
 * 
 * @tparam T 
 * @tparam V 
 * @param inputs 
 * @return std::vector<T> 
 */
template <typename T, typename V>
std::vector<T> ODESystem<T, V>::_eval(std::vector<V>& inputs) 
{
    std::vector<T> res;
    for (function_t<T, V> &func : ODESystem<T, V>::_functions)
        res.push_back(func(inputs));

    return res;
}





/**
 * @brief Solve an ordinary differential equation numerically.
 * 
 * @tparam T 
 * @tparam V 
 * @param eq 
 * @param alg 
 * @return DataFrame<T> 
 */
template <typename T, typename V>
DataFrame<T> solve(ODE<T, V>& eq, algorithm_t alg) 
{

    std::cout << alg << std::endl;
    switch (alg)
    {
    case ALGORITHIM_EULER:
        return _euler(eq);
        break;
    
    case ALGORITHM_RK4:
        return _RK4(eq);

    default:
        break;
    }
}



template <typename T, typename V>
DataFrame<T> solve(ODESystem<T, V> const& eq, algorithm_t alg) 
{
    switch (alg)
    {
    case ALGORITHIM_EULER:
        return _euler(eq);

    case ALGORITHM_RK4:
        return _RK4(eq);
    
    default:
        break;
    }
}



}

#endif