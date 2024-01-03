#ifndef ALGORITHMS_H
#define ALGORITHMS_H


#include "dataframe.h"
#include "ode.h"


namespace DES
{

template <typename T, typename V>
DataFrame<T> _RK4(ODESystem<T, V> const& ode) 
{
    // if (ode._functions.size() == 1) 
    //     return _RK4(ODE<T, V>(ode.))

    DataFrame<T> res(0, ode._equations + 1); // +1 for time column
    float h = (ode._timeBound.second - ode._timeBound.first) / ode._timeStep;
    float t = ode._timeBound.first;
    size_t m = ode._equations;

    // add the first row (initial values)
    res.addRow(ode._iValues);

    // keep track of the last set of values
    std::vector<T> inputs = res.getRow(0);
    size_t row = 0;

    while (t < ode._timeBound.second)
    {
        T k_1j[m], k_2j[m], k_3j[m], k_4j[m];

        // k1
        for (int j = 0; j < m; j++) 
            k_1j[j] = ode._eval(inputs);

        // modify eval inputs
        for (int i = 0; i < inputs.size(); i++)
            inputs[i] += (i == 0 ? 0.5 * h : 0.5 * k_1j[i-1]);

        // k2
        for (int j = 0; j < m; j++) 
            k_2j[j] = ode._eval(inputs);

        // modify again
        for (int i = 0; i < inputs.size(); i++)
            inputs[i] += (i == 0 ? 0.5 * h : 0.5 * k_2j[i-1]);

        // k3
        for (int j = 0; j < m; j++) 
            k_3j[j] = ode._eval(inputs);

        // modify for the last time
        for (int i = 0; i < inputs.size(); i++)
            inputs[i] += (i == 0 ? h : k_3j[i-1]);
        
        // k4
        for (int j = 0; j < m; j++) 
            k_4j[j] = ode._eval(inputs);

        
        std::vector<T> result = res.getRow(row);
        for (int j = 0; j < m; j++) 
            result += h * (k_1j[j] + 2 * k_2j[j] + 2 * k_3j[j] + k_4j[j]) / 6;

        res.addRow(result);
        row++;
        t += h;
    }
    
}



} // namespace DES



#endif