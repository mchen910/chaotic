#ifndef DIFFEQ_SOLVER_H
#define DIFFEQ_SOLVER_H

#include <initializer_list>
#include <functional>
#include <vector>

#include "dataframe.h"


namespace DES
{

typedef unsigned int                algorithm_t;


/**
 * @brief Time bound type, encapsulating information about when to start and 
 * stop the solver. The bounds [a, b) are inclusive at the start and exclusive 
 * at the end.
 * 
 * @tparam T Type of time bound, which can be used to specify precision.
 */
template <typename T>
struct timeBound_t
{
    T first, second;

    timeBound_t(T start, T end) : first(start), second(end) { }
    timeBound_t(std::initializer_list<T> bounds) : first(bounds.begin()[0]), second(bounds.begin()[1]) { }
};


template <typename T>
struct iv_t
{
    std::vector<T> vec;

    iv_t(std::vector<T> values)             : vec(values) { }
    iv_t(std::initializer_list<T> values)   : vec(values) { }
};




/**
 * @brief std::function wrapper. Enforces a vector of inputs.
 * 
 * @tparam T Input/output type of function - function_t is structured such that all 
 * inputs to the function must be of one type, and the output must be of the same 
 * type in order to keep things consistent.
 */
template <typename T>
struct function_t
{
    std::function<T(std::vector<T>)> _func;

    function_t(std::function<T(std::vector<T>)> func) {
        _func = func;
    }

    T operator()(std::vector<T> args) 
    {
        return _func(args);
    }
};



/**
 * @brief Base class for a differential equation. The ODE and PDE 
 * classes both inherit from this class. The most general form of a 
 * differential equation is dy/dt = f(t, ...).
 * 
 * @tparam T Type of inputs/output to the differential equation.
 */
template <typename T>
class DiffEq
{

protected:
    function_t<T> _func;
    timeBound_t<T> _bounds;
    iv_t<T> _iCondition;
    T _timeStep;

public:
    /**
     * @brief Construct a new Diff Eq object.
     * 
     * @param func 
     * @param bounds 
     * @param initialCondition 
     * @param timeStep 
     */
    DiffEq(function_t<T>& func, timeBound_t<T>& bounds, iv_t<T>& initialCondition, T timeStep)
        : _func(func)   // must initialize here to avoid explicit initialization errors
        , _iCondition(initialCondition)
        , _bounds(bounds)
        , _timeStep(timeStep)
    { }

};


template <typename T>
class DiffEqSystem
{

protected:
    std::vector<function_t<T>> _functions;
    timeBound_t<T> _timeBound;
    iv_t<T> _iValues;
    T _timeStep;

public:
    DiffEqSystem(iv_t<T>& iValues, std::initializer_list<function_t<T>> funcs, timeBound_t<T>& bounds, T timeStep) 
        : _functions(funcs)
        , _iValues(iValues)
        , _timeBound(bounds)
        , _timeStep(timeStep)
    { } 

};

} // namespace DES


#endif