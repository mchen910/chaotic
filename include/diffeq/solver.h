#ifndef SOLVER_H
#define SOLVER_H

#include <initializer_list>
#include <functional>
#include <vector>

#include "dataframe.h"

#include <iostream>


namespace DES
{

typedef unsigned int                algorithm_t;
typedef std::pair<float, float>     timeBound_t;
typedef std::vector<float>          iv_t;



/**
 * @brief std::function wrapper. Enforces a vector of inputs.
 * 
 * @tparam T Output type of function
 * @tparam V Input type of function - function_t is structured such that all 
 * inputs to the function must be of one type.
 */
template <typename T, typename V>
struct function_t
{
    std::function<T(std::vector<V>)> _func;

    function_t(std::function<T(std::vector<V>)> func) {
        _func = func;
    }

    T operator()(std::vector<V> args) 
    {
        return _func(args);
    }
};



/**
 * @brief Base class for a differential equation. The ODE and PDE 
 * classes both inherit from this class. The most general form of a 
 * differential equation is dy/dt = f(t, ...).
 * 
 * @tparam T Type of output of the function f.
 * @tparam V Type of inputs to the function f.
 */
template <typename T, typename V>   
class DiffEq
{

protected:
    function_t<T, V> func;
    timeBound_t bounds;
    double initialCondition;
    double timeStep;

public:
    /**
     * @brief Construct a new Diff Eq object.
     * 
     * @param func 
     * @param bounds 
     * @param initialCondition 
     * @param timeStep 
     */
    DiffEq(function_t<T, V>& func, timeBound_t& bounds, double initialCondition, double timeStep) {
        this->func = func;
        this->initialCondition = initialCondition;
        this->bounds = bounds;
        this->timeStep = timeStep;
    }

};


template <typename T, typename V>
class DiffEqSystem
{

public:
    std::vector<function_t<T, V>> _functions;
    std::vector<float> _iValues;

    DiffEqSystem(iv_t& iValues, std::initializer_list<function_t<T, V>> funcs) 
        : _functions(funcs)
    {
        _iValues = iValues;
    } 


    std::vector<T> eval(std::vector<V>& inputs); 
    std::vector<function_t<T, V>> getFunctions();

};



} // namespace DES

#endif