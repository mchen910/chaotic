#ifndef DIFFEQ_ODE_H
#define DIFFEQ_ODE_H

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
template <typename T>
class ODE : public DiffEq<T>
{

public:
    ODE(function_t<T>& func, timeBound_t<T>& bounds, T initialCondition, T timeStep) 
    : DiffEq<T>(func, bounds, initialCondition, timeStep) { }

};


/**
 * @brief Solve a system of first order ordinary differential equations, all of the 
 * form dy_j/dt = f_j(t, y_1, y_2, ...). 
 * 
 * @tparam T 
 */
template <typename T>
class ODESystem : public DiffEqSystem<T>
{

private:
    timeBound_t<T> _timeBound;
    T _timeStep;
    size_t _equations;

public:
    ODESystem(iv_t<T>& iValues, std::initializer_list<function_t<T>> funcs, timeBound_t<T>& bounds, T timeStep)
        : DiffEqSystem<T>(iValues, funcs)
        , _timeBound(bounds)    // must be done to avoid explicit initialization errors
    { 
        _timeStep = timeStep;
        _equations = funcs.size();
    };

    std::vector<T> _eval(std::vector<T>& inputs); 
    
    size_t          getNumEquations();      
    timeBound_t<T>  getTimeBound();
    T               getTimeStep();
};



// Forward declare the algorithms

template <typename T>   DataFrame<T>  _EULER    (ODESystem<T>& ode);
template <typename T>   DataFrame<T>  _RK4      (ODESystem<T>& ode);
template <typename T>   DataFrame<T>  _RKF45    (ODESystem<T>& ode);
template <typename T>   DataFrame<T>  _TSIT5    (ODESystem<T>& ode);

template <typename T>   DataFrame<T>  _RKF45    (ODESystem<T>& ode, T maxError);


/**
 * @brief Evaluate an entire system of functions with one set of 
 * inputs. Useful for coupled Runge-Kutta schemes, which use common sets  
 * of inputs for every iteration.
 * 
 * @tparam T 
 * @param inputs 
 * @return std::vector<T> 
 */
template <typename T>
std::vector<T> ODESystem<T>::_eval(std::vector<T>& inputs) 
{
    std::vector<T> res;
    for (function_t<T> &func : ODESystem<T>::_functions)
        res.push_back(func(inputs));

    return res;
}


/**
 * @brief Get the number of equations in an ODESystem object.
 * @tparam T 
 * @return size_t Number of equations present.
 */
template <typename T>
size_t ODESystem<T>::getNumEquations() 
{
    return _equations;
}


/**
 * @brief Get the time bounds in an ODESystem object.
 * 
 * @tparam T 
 * @return timeBound_t 
 */
template <typename T>
timeBound_t<T> ODESystem<T>::getTimeBound()
{
    return _timeBound;
}


/**
 * @brief Get the timestep that the ODESystem object propagates forward by.
 * 
 * @tparam T 
 * @tparam V 
 * @return double 
 */
template <typename T>
T ODESystem<T>::getTimeStep()
{
    return _timeStep;
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
template <typename T>
DataFrame<T> solve(ODE<T>& eq, algorithm_t alg) 
{

    std::cout << alg << std::endl;
    switch (alg)
    {
    case ALGORITHM_EULER:
        return _EULER(eq);
        break;
    
    case ALGORITHM_RK4:
        return _RK4(eq);

    default:
        break;
    }
}


template <typename T>
DataFrame<T> solve(ODESystem<T>& eq, algorithm_t alg) 
{
    switch (alg)
    {
        case ALGORITHM_EULER:   return _EULER(eq);
        case ALGORITHM_RK4:     return _RK4(eq);
        case ALGORITHM_RKF45:   return _RKF45(eq);
        
        default:                throw std::runtime_error("Invalid algorithm");
    }
}



}


#endif