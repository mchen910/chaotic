#include "solver.h"

namespace DES {

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