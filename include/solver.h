#ifndef SOLVER_H
#define SOLVER_H

#include <functional>
#include <vector>

#include "dataframe.h"


#define ALGORITHIM_EULER    0x001
#define ALGORITHM_RK4       0x002
#define ALGORITHM_RK45      0x003



namespace DES
{

typedef unsigned int                                uint_t;
typedef unsigned int                                algorithm_t;
typedef std::pair<float, float>                     timeBound_t;


template <typename T, typename V>
struct function_t
{
    std::function<T(vector<V>)> _func;
    std::vector<V> _inputs;

    T operator(...V args) : _inputs { args } 
    {
        return _func(_inputs);
    }
};


template <typename T, typename V>   
class DiffEq
{

protected:
    function_t<T, V> func;
    timeBound_t bounds;
    double initialCondition;

public:
    DiffEq(function_t<T, V> func, double initialCondition, timeBound_t bounds) {
        this->func = func;
        this->initialCondition = initialCondition;
        this->bounds = bounds;
    }
};


template <typename T, typename V>
DataFrame<T> solve(DiffEq<T, V> eq, algorithm_t alg);



} // namespace solver



#endif