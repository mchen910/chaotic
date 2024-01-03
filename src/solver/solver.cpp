#include "solver.h"

namespace DES {


template <typename T, typename V>
std::vector<T> System<T, V>::eval(std::vector<V>& inputs) 
{
    std::vector<T> res;
    for (function_t &func : _functions)
        res.push_back(func(inputs));

    return res;
}


template <typename T, typename V>
DataFrame<T> solve(DiffEq<T, V> eq, algorithm_t alg) 
{
    switch (alg)
    {
    case ALGORITHIM_EULER:
        return _euler(eq);
        break;
    
    default:
        break;
    }
}





};