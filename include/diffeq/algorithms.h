#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "dataframe.h"
#include "ode.h"

#include <iostream>

namespace DES
{


struct ButcherTableau
{
    // Butcher Tableaus have empty entries, use a macro for readability
    #define TAB_NULL    0.0
    #define TAB_NULLF   0.0f


    const double RK38_TAB[5][5] {
        {       0.0,   TAB_NULL,  TAB_NULL,  TAB_NULL,  TAB_NULL   },
        {   1.0/3.0,    1.0/3.0,  TAB_NULL,  TAB_NULL,  TAB_NULL   },
        {   2.0/3.0,   -1.0/3.0,       1.0,  TAB_NULL,  TAB_NULL   },
        {       1.0,        1.0,      -1.0,       1.0,  TAB_NULL   },
        {  TAB_NULL,    1.0/8.0,   3.0/8.0,   3.0/8.0,   3.0/8.0   }
    };


    const double RK38_TABF[5][5] {
        {       0.0f,   TAB_NULLF,  TAB_NULLF,  TAB_NULLF,  TAB_NULLF   },
        {   1.0/3.0f,    1.0/3.0f,  TAB_NULLF,  TAB_NULLF,  TAB_NULLF   },
        {   2.0/3.0f,   -1.0/3.0f,       1.0f,  TAB_NULLF,  TAB_NULLF   },
        {       1.0f,        1.0f,      -1.0f,       1.0f,  TAB_NULLF   },
        {  TAB_NULLF,    1.0/8.0f,   3.0/8.0f,   3.0/8.0f,   3.0/8.0f   }
    };


    const double RKF45_TAB[8][7] {
        {         0.0,        TAB_NULL,         TAB_NULL,         TAB_NULL,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {     1.0/4.0,         1.0/4.0,         TAB_NULL,         TAB_NULL,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {     3.0/8.0,        3.0/32.0,         9.0/32.0,         TAB_NULL,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {   12.0/13.0,   1932.0/2197.0,   -7200.0/2197.0,    7296.0/2197.0,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {         1.0,     439.0/216.0,             -8.0,     3680.0/513.0,     -845.0/4104.0,     TAB_NULL,   TAB_NULL    },
        {     1.0/2.0,       -8.0/27.0,              2.0,   -3544.0/2565.0,     1859.0/4104.0,   -11.0/40.0,   TAB_NULL    },
        {    TAB_NULL,      25.0/216.0,              0.0,    1408.0/2565.0,     2197.0/1404.0,     -1.0/5.0,        0.0    },
        {    TAB_NULL,      16.0/135.0,              0.0,   6656.0/12825.0,   28561.0/56430.0,    -9.0/50.0,   2.0/55.0    }
    };

};

template <typename T>
DataFrame<T> _EULER(ODESystem<T> & ode) 
{
    return DataFrame<T>(1, 2);
}


template <typename T>
DataFrame<T> _RK4(ODESystem<T>& ode) 
{
    // if (ode._functions.size() == 1) 
    //     return _RK4(ODE<T, V>(ode.))

    timeBound_t<T> tBound = ode.getTimeBound();
    size_t m = ode.getNumEquations();
    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1); // +1 for time column

    // add the first row (initial values)
    res.addRow(ode._iValues.vec);

    // keep track of the last set of values
    size_t row = 0;


    while (t < tBound.second)
    {
        std::vector<T> k_1j, k_2j, k_3j, k_4j;
        std::vector<T> inputs(res.getRow(row));

        k_1j = ode._eval(inputs);   // k1

        // modify eval inputs
        for (int i = 0; i < inputs.size(); i++)
            inputs[i] += (i == 0 ? 0.5 * h : 0.5 * h * k_1j[i-1]);

        k_2j = ode._eval(inputs);   // k2

        // modify again
        for (int i = 0; i < inputs.size(); i++)
            inputs[i] += (i == 0 ? 0.0 : 0.5 * h * (k_2j[i-1] - k_1j[i-1]));

        k_3j = ode._eval(inputs);   // k3

        // modify for the last time
        for (int i = 0; i < inputs.size(); i++)
            inputs[i] += (i == 0 ? 0.5 * h : 0.5 * h * (2 * k_3j[i-1] - k_2j[i-1]));

        k_4j = ode._eval(inputs);   // k4
        
        std::vector<T> result(res.getRow(row));
        result[0] = t + h;
        for (int j = 1; j <= m; j++) 
            result[j] += h * (k_1j[j-1] + 2 * k_2j[j-1] + 2 * k_3j[j-1] + k_4j[j-1]) / 6;

        res.addRow(result);
        
        row++;
        t += h;
    }

    return res;
}


template <typename T>
DataFrame<T> _RK45(ODESystem<T>& ode) 
{

}















template <typename T>
void printVector(const std::vector<T>& v)
{
    for (const auto& i: v)
        std::cout << i << ' ';
    std::cout << '\n';
}



} // namespace DES



#endif