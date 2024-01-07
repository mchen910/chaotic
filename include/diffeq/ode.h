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
    ODE(function_t<T>& func, timeBound_t<T>& bounds, iv_t<T>& initialCondition, T timeStep) 
    : DiffEq<T>(func, bounds, initialCondition, timeStep) { }

    T               _eval(std::vector<T>& input);
    T               getTimeStep();
    timeBound_t<T>  getTimeBound();
    iv_t<T>         getInitialCondition();

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
    size_t _equations;

public:
    ODESystem() = default;
    ODESystem(iv_t<T>& iValues, std::initializer_list<function_t<T>> funcs, timeBound_t<T>& bounds, T timeStep)
        : DiffEqSystem<T>(iValues, funcs, bounds, timeStep)
    { 
        _equations = funcs.size();
    };

    std::vector<T> _eval(std::vector<T>& inputs);  
    iv_t<T>         getInitialConditions();    
    timeBound_t<T>  getTimeBound();
    T               getTimeStep();
    size_t          getNumEquations(); 

};



template <typename T>   DataFrame<T>  _EULER    (ODE<T>& ode);
template <typename T>   DataFrame<T>  _RK4      (ODE<T>& ode);
template <typename T>   DataFrame<T>  _RKF45    (ODE<T>& ode);
template <typename T>   DataFrame<T>  _RKF45    (ODE<T>& ode, T maxError);
template <typename T>   DataFrame<T>  _TSIT5    (ODE<T>& ode);

template <typename T>   DataFrame<T>  _EULER    (ODESystem<T>& ode);
template <typename T>   DataFrame<T>  _RK4      (ODESystem<T>& ode);
template <typename T>   DataFrame<T>  _RKF45    (ODESystem<T>& ode);
template <typename T>   DataFrame<T>  _RKF45    (ODESystem<T>& ode, T maxError);
template <typename T>   DataFrame<T>  _TSIT5    (ODESystem<T>& ode);

template <typename T>   std::vector<T>  _EULER_i  (ODE<T>& ode);
template <typename T>   std::vector<T>  _RK4_i    (ODE<T>& ode);
template <typename T>   std::vector<T>  _RKF45_i  (ODE<T>& ode);
template <typename T>   std::vector<T>  _RKF45_i  (ODE<T>& ode, T maxError);
template <typename T>   std::vector<T>  _TSIT5_i  (ODE<T>& ode);

template <typename T>   std::vector<T>  _EULER_i  (ODESystem<T>& ode);
template <typename T>   std::vector<T>  _RK4_i    (ODESystem<T>& ode);
template <typename T>   std::vector<T>  _RKF45_i  (ODESystem<T>& ode);
template <typename T>   std::vector<T>  _RKF45_i  (ODESystem<T>& ode, T maxError);
template <typename T>   std::vector<T>  _TSIT5_i  (ODESystem<T>& ode);



template <typename T>
T ODE<T>::_eval(std::vector<T>& input) 
{
    return this->_func(input);
}


template <typename T>
timeBound_t<T> ODE<T>::getTimeBound() 
{
    return this->_bounds;
}


template <typename T>
iv_t<T> ODE<T>::getInitialCondition()
{
    return this->_iCondition;
}


template <typename T>
T ODE<T>::getTimeStep() 
{
    return this->_timeStep;
}


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


template <typename T>
iv_t<T> ODESystem<T>::getInitialConditions() 
{
    return this->_iValues;
}


/**
 * @brief Get the number of equations in an ODESystem object.
 * @tparam T 
 * @return size_t Number of equations present.
 */
template <typename T>
size_t ODESystem<T>::getNumEquations() 
{
    return this->_equations;
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
    return this->_timeBound;
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
    return this->_timeStep;
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
    switch (alg)
    {
        case ALGORITHM_EULER:   return _EULER(eq);
        case ALGORITHM_RK4:     return _RK4(eq);
        case ALGORITHM_RKF45:   return _RKF45(eq);

        default:                throw std::runtime_error("Invalid algorithm");
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


template <typename T>
std::vector<T> solve_i(ODE<T>& eq, algorithm_t alg) 
{
    switch (alg)
    {
        case ALGORITHM_EULER:   return _EULER_i(eq);
        case ALGORITHM_RK4:     return _RK4_i(eq);
        case ALGORITHM_RKF45:   return _RKF45_i(eq);

        default:                throw std::runtime_error("Invalid algorithm");
    }
}


template <typename T>
std::vector<T> solve_i(ODESystem<T>& eq, algorithm_t alg)
{
    switch (alg)
    {
        case ALGORITHM_EULER:   return _EULER_i(eq);
        case ALGORITHM_RK4:     return _RK4_i(eq);
        case ALGORITHM_RKF45:   return _RKF45_i(eq);

        default:                throw std::runtime_error("Invalid algorithm");
    }
}


} // namespace DES


#endif