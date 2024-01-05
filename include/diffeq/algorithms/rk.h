#ifndef DIFFEQ_ALGORITHMS_RK_H
#define DIFFEQ_ALGORITHMS_RK_H

#include <cmath>

#include "../ode.h"
#include "tableau.h"


#define  _LOOP_TO_M(i, n)       for (size_t i = n; i <= m; i++)
#define  _ABS(x)                (x >= 0 ? x : -(x))
#define  _MIN(a, b)             (a > b ? b : a)


namespace DES 
{


template <typename T>
DataFrame<T> _EULER(ODESystem<T> & ode) 
{
   timeBound_t<T> tBound = ode.getTimeBound();

    size_t m = ode.getNumEquations();
    size_t row = 0;

    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1);
    res.addRow(ode.getInitialConditions().vec);


#if defined(DIFFEQ_FLOAT_PRECISION)
    #define  _EULER_TIME(a, b)        EULER_TABF[b][0] - EULER_TABF[a][0]
    #define  _EULER_UNIT(a, b)        EULER_TABF[a][b] * k_##b##j[i-1]
#elif defined(DIFFEQ_DOUBLE_PRECISION)
    #define  _EULER_TIME(a, b)        EULER_TAB[b][0]  - EULER_TAB[a][0]
    #define  _EULER_UNIT(a, b)        EULER_TAB[a][b]  * k_##b##j[i-1] 
#elif defined(DIFFEQ_LONG_DOUBLE_PRECISION)
    #define  _EULER_TIME(a, b)        EULER_TABL[b][0] - EULER_TABL[a][0]
    #define  _EULER_UNIT(a, b)        EULER_TABL[a][b] * k_##b##j[i-1]
#else 
    #define  _EULER_TIME(a, b)        (T)(EULER_TAB[b][0]) - (T)(EULER_TAB[b][0])
    #define  _EULER_UNIT(a, b)        (T)(EULER_TAB[a][b] * k_##b##j[i-1])
#endif

    using namespace ButcherTableau;
    do
    {
        std::vector<T> k_1j, k_2j, k_3j, k_4j;
        std::vector<T> inputs(res.getRow(row));

        k_1j = ode._eval(inputs);
        
        std::vector<T> result(res.getRow(row));
        result[0] = t + h;

        _LOOP_TO_M(i, 1)
            result[i] += h * (_EULER_UNIT(1, 1));

        res.addRow(result);
        
        row++;
        t += h;
        
    } while (t < tBound.second);

    return res;
}



template <typename T>
DataFrame<T> _RK4(ODE<T>& ode) 
{
    timeBound_t<T> tBound = ode.getTimeBound();
    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, 2);
    res.addRow(ode.getInitialCondition().vec);

    size_t row = 0;
    
    do
    {
        T k_1, k_2, k_3, k_4;
        std::vector<T> inputs(res.getRow(row));

        k_1 = ode._eval(inputs);

        inputs[0] += h / 2;
        inputs[1] += h * k_1 / 2;

        k_2 = ode._eval(inputs);

        inputs[1] += h * (k_2 / 2 - k_1 / 2);

        k_3 = ode._eval(inputs);

        inputs[0] += h / 2;
        inputs[1] += h * (k_3 - k_2 / 2);

        k_4 = ode._eval(inputs);
        
        std::vector<T> result(res.getRow(row));
        result[0] += h;
        result[1] += h * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;

        res.addRow(result);
        
        row++;
        t += h;
        
    } while (t < tBound.second);

    return res;
}



template <typename T>
DataFrame<T> _RK4(ODESystem<T>& ode) 
{
    timeBound_t<T> tBound = ode.getTimeBound();

    size_t m = ode.getNumEquations();
    size_t row = 0;

    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1);
    res.addRow(ode.getInitialConditions().vec);


#if defined(DIFFEQ_FLOAT_PRECISION)
    #define  _RK4_TIME(a, b)        RK4_TABF[b][0] - RK4_TABF[a][0]
    #define  _RK4_UNIT(a, b)        RK4_TABF[a][b] * k_##b##j[i-1]
#elif defined(DIFFEQ_DOUBLE_PRECISION)
    #define  _RK4_TIME(a, b)        RK4_TAB[b][0]  - RK4_TAB[a][0]
    #define  _RK4_UNIT(a, b)        RK4_TAB[a][b]  * k_##b##j[i-1] 
#elif defined(DIFFEQ_LONG_DOUBLE_PRECISION)
    #define  _RK4_TIME(a, b)        RK4_TABL[b][0] - RK4_TABL[a][0]
    #define  _RK4_UNIT(a, b)        RK4_TABL[a][b] * k_##b##j[i-1]
#else 
    #define  _RK4_TIME(a, b)        (T)(RK4_TAB[b][0]) - (T)(RK4_TAB[b][0])
    #define  _RK4_UNIT(a, b)        (T)(RK4_TAB[a][b] * k_##b##j[i-1])
#endif

    using namespace ButcherTableau;
    do
    {
        std::vector<T> k_1j, k_2j, k_3j, k_4j;
        std::vector<T> inputs(res.getRow(row));

        k_1j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RK4_TIME(0, 1) : _RK4_UNIT(1, 1));

        k_2j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RK4_TIME(1, 2) : 
                ( _RK4_UNIT(2, 1) + _RK4_UNIT(2, 2)
                - _RK4_UNIT(1, 1)));

        k_3j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RK4_TIME(2, 3) :
                ( _RK4_UNIT(3, 1) + _RK4_UNIT(3, 2) + _RK4_UNIT(3, 3)
                - _RK4_UNIT(2, 1) - _RK4_UNIT(2, 2)));

        k_4j = ode._eval(inputs);
        
        std::vector<T> result(res.getRow(row));
        result[0] = t + h;

        _LOOP_TO_M(i, 1)
            result[i] += h * (
                  _RK4_UNIT(4, 1)
                + _RK4_UNIT(4, 2)
                + _RK4_UNIT(4, 3)
                + _RK4_UNIT(4, 4)
            );

        res.addRow(result);
        
        row++;
        t += h;
        
    } while (t < tBound.second);

    return res;
}



template <typename T>
std::vector<T> _RK4_i(ODE<T>& ode) 
{
    // The only things that are needed to propagate through time are the last 
    // row, which contains the time information as well
    static std::vector<T> lastRow = ode.getInitialCondition().vec;

    T h = ode.getTimeStep();

    T k_1, k_2, k_3, k_4;
    std::vector<T> inputs(lastRow);

    k_1 = ode._eval(inputs);

    inputs[0] += h / 2;
    inputs[1] += h * k_1 / 2;

    k_2 = ode._eval(inputs);

    inputs[1] += h * (k_2 / 2 - k_1 / 2);

    k_3 = ode._eval(inputs);

    inputs[0] += h / 2;
    inputs[1] += h * (k_3 - k_2 / 2);

    k_4 = ode._eval(inputs);
        
    lastRow[0] += h;
    lastRow[1] += h * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;

    return std::vector<T>(lastRow);
}



template <typename T>
std::vector<T> _RK4_i(ODESystem<T>& ode) 
{
    static std::vector<T> lastRow = ode.getInitialConditions().vec;
    T h = ode.getTimeStep();

#if defined(DIFFEQ_FLOAT_PRECISION)
    #define  _RK4_TIME(a, b)        RK4_TABF[b][0] - RK4_TABF[a][0]
    #define  _RK4_UNIT(a, b)        RK4_TABF[a][b] * k_##b##j[i-1]
#elif defined(DIFFEQ_DOUBLE_PRECISION)
    #define  _RK4_TIME(a, b)        RK4_TAB[b][0]  - RK4_TAB[a][0]
    #define  _RK4_UNIT(a, b)        RK4_TAB[a][b]  * k_##b##j[i-1] 
#elif defined(DIFFEQ_LONG_DOUBLE_PRECISION)
    #define  _RK4_TIME(a, b)        RK4_TABL[b][0] - RK4_TABL[a][0]
    #define  _RK4_UNIT(a, b)        RK4_TABL[a][b] * k_##b##j[i-1]
#else 
    #define  _RK4_TIME(a, b)        (T)(RK4_TAB[b][0]) - (T)(RK4_TAB[b][0])
    #define  _RK4_UNIT(a, b)        (T)(RK4_TAB[a][b] * k_##b##j[i-1])
#endif

    using namespace ButcherTableau;
    std::vector<T> k_1j, k_2j, k_3j, k_4j;
    std::vector<T> inputs(lastRow);

    k_1j = ode._eval(inputs);

    _LOOP_TO_M(i, 0)
        inputs[i] += h * (i == 0 ? _RK4_TIME(0, 1) : _RK4_UNIT(1, 1));

    k_2j = ode._eval(inputs);

    _LOOP_TO_M(i, 0)
        inputs[i] += h * (i == 0 ? _RK4_TIME(1, 2) : 
            ( _RK4_UNIT(2, 1) + _RK4_UNIT(2, 2)
            - _RK4_UNIT(1, 1)));

    k_3j = ode._eval(inputs);

    _LOOP_TO_M(i, 0)
        inputs[i] += h * (i == 0 ? _RK4_TIME(2, 3) :
            ( _RK4_UNIT(3, 1) + _RK4_UNIT(3, 2) + _RK4_UNIT(3, 3)
            - _RK4_UNIT(2, 1) - _RK4_UNIT(2, 2)));

    k_4j = ode._eval(inputs);
    
    lastRow[0] += h;

    _LOOP_TO_M(i, 1)
        lastRow[i] += h * (
              _RK4_UNIT(4, 1)
            + _RK4_UNIT(4, 2)
            + _RK4_UNIT(4, 3)
            + _RK4_UNIT(4, 4)
        );

    return std::vector<T>(lastRow);
}



template <typename T>
DataFrame<T> _RKF45(ODE<T>& ode, T maxError) 
{
    timeBound_t<T> tBound = ode.getTimeBound();
    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, 2);
    res.addRow(ode.getInitialCondition().vec);

    size_t row = 0;

#if defined(DIFFEQ_FLOAT_PRECISION)
    #define  _RKF45_ENTRY(a, b)         RKF45_TABF[a][b]
    #define  _RKF45_POW(a, b)           powf(a, b)
#elif defined(DIFFEQ_DOUBLE_PRECISION)
    #define  _RKF45_ENTRY(a, b)         RKF45_TAB[a][b]  
    #define  _RKF45_POW(a, b)           pow(a, b)
#elif defined(DIFFEQ_LONG_DOUBLE_PRECISION)
    #define  _RKF45_ENTRY(a, b)         RKF45_TABL[a][b]
    #define  _RKF45_POW(a, b)           powl(a, b)
#else 
    #define  _RKF45_ENTRY(a, b)         (T)(RKF45_TAB[a][b])
    #define  _RKF45_POW(a, b)           (T)pow(a, b)
#endif


    using namespace ButcherTableau;
    while (t < tBound.second)
    {
        // h = _MIN(h, tBound.second - t); // ensure that the last entry stops when h exceeds the upper bound
        T k_1, k_2, k_3, k_4, k_5, k_6;
        std::vector<T> inputs(res.getRow(row));

        k_1 = ode._eval(inputs);

        inputs[0] += h * _RKF45_ENTRY(1, 0);
        inputs[1] += h * _RKF45_ENTRY(1, 1) * k_1;

        k_2 = ode._eval(inputs);

        inputs[0] += h * (_RKF45_ENTRY(2, 0) - _RKF45_ENTRY(1, 0));
        inputs[1] += h * (_RKF45_ENTRY(2, 1) * k_1 + _RKF45_ENTRY(2, 2) * k_2 
            - _RKF45_ENTRY(1, 1) * k_1);

        k_3 = ode._eval(inputs);

        inputs[0] += h * (_RKF45_ENTRY(3, 0) - _RKF45_ENTRY(2, 0));
        inputs[1] += h * (_RKF45_ENTRY(3, 1) * k_1 + _RKF45_ENTRY(3, 2) * k_2 + _RKF45_ENTRY(3, 3) * k_3
            - _RKF45_ENTRY(2, 1) * k_1 - _RKF45_ENTRY(2, 2) * k_2);

        k_4 = ode._eval(inputs);

        inputs[0] += h * (_RKF45_ENTRY(4, 0) - _RKF45_ENTRY(3, 0));
        inputs[1] += h * (_RKF45_ENTRY(4, 1) * k_1 + _RKF45_ENTRY(4, 2) * k_2 + _RKF45_ENTRY(4, 3) * k_3 + _RKF45_ENTRY(4, 4) * k_4 
            - _RKF45_ENTRY(3, 1) * k_1 - _RKF45_ENTRY(3, 2) * k_2 - _RKF45_ENTRY(3, 3) * k_3);

        k_5 = ode._eval(inputs);

        inputs[0] += h * (_RKF45_ENTRY(5, 0) - _RKF45_ENTRY(4, 0));
        inputs[1] += h * (_RKF45_ENTRY(5, 1) * k_1 + _RKF45_ENTRY(5, 2) * k_2 + _RKF45_ENTRY(5, 3) * k_3 + _RKF45_ENTRY(5, 4) * k_4 + _RKF45_ENTRY(5, 5) * k_5
            - _RKF45_ENTRY(4, 1) * k_1 - _RKF45_ENTRY(4, 2) * k_2 - _RKF45_ENTRY(4, 3) * k_3 - _RKF45_ENTRY(4, 4) * k_4);

        k_6 = ode._eval(inputs);

        std::vector<T> result(res.getRow(row));
        
        T w1 = result[1] + h * (_RKF45_ENTRY(6, 1) * k_1 + _RKF45_ENTRY(6, 2) * k_2 + _RKF45_ENTRY(6, 3) * k_3 + _RKF45_ENTRY(6, 4) * k_4 + _RKF45_ENTRY(6, 5) * k_5 + _RKF45_ENTRY(6, 6) * k_6);
        T w2 = result[1] + h * (_RKF45_ENTRY(7, 1) * k_1 + _RKF45_ENTRY(7, 2) * k_2 + _RKF45_ENTRY(7, 3) * k_3 + _RKF45_ENTRY(7, 4) * k_4 + _RKF45_ENTRY(7, 5) * k_5 + _RKF45_ENTRY(7, 6) * k_6);
        
        T R = _ABS(w2 - w1) / h;
        T delta = 0.84 * _RKF45_POW(maxError / R, 0.25);

        if (R <= maxError) {
            result[0] += h;
            result[1] = w1;
            res.addRow(result);

            h *= delta;
            t += h;
            row++;
        
        } else {
            h *= delta;
        }

    }

    return res;

}



template <typename T>
DataFrame<T> _RKF45(ODE<T>& ode) 
{
    return _RKF45(ode, (T)DEFAULT_MAX_ERROR);
}



template <typename T>
DataFrame<T> _RKF45(ODESystem<T>& ode, T maxError) 
{
    timeBound_t<T> tBound = ode.getTimeBound();

    size_t m = ode.getNumEquations();
    size_t row = 0;

    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1);
    res.addRow(ode.getInitialConditions().vec);


#if defined(DIFFEQ_FLOAT_PRECISION)
    #define  _RKF45_TIME(a, b)        RKF45_TABF[b][0] - RKF45_TABF[a][0]
    #define  _RKF45_UNIT(a, b)        RKF45_TABF[a][b] * k_##b##j[i-1]
    #define  _RKF45_SQRT(a)           sqrtf(a)
    #define  _RKF45_POW(a, b)         powf(a, b)
#elif defined(DIFFEQ_DOUBLE_PRECISION)
    #define  _RKF45_TIME(a, b)        RKF45_TAB[b][0]  - RKF45_TAB[a][0]
    #define  _RKF45_UNIT(a, b)        RKF45_TAB[a][b]  * k_##b##j[i-1] 
    #define  _RKF45_SQRT(a)           sqrt(a)
    #define  _RKF45_POW(a, b)         pow(a, b)
#elif defined(DIFFEQ_LONG_DOUBLE_PRECISION)
    #define  _RKF45_TIME(a, b)        RKF45_TABL[b][0] - RKF45_TABL[a][0]
    #define  _RKF45_UNIT(a, b)        RKF45_TABL[a][b] * k_##b##j[i-1]
    #define  _RKF45_SQRT(a)           sqrtl(a)
    #define  _RKF45_POW(a, b)         powl(a, b)
#else 
    #define  _RKF45_TIME(a, b)        (T)(RKF45_TAB[b][0]) - (T)(RKF45_TAB[b][0])
    #define  _RKF45_UNIT(a, b)        (T)(RKF45_TAB[a][b] * k_##b##j[i-1])
    #define  _RKF45_SQRT(a)           (T)sqrt(a)
    #define  _RKF45_POW(a, b)         (T)pow(a, b)
#endif

    using namespace ButcherTableau;
    while (t < tBound.second)
    {
        std::vector<T> k_1j, k_2j, k_3j, k_4j, k_5j, k_6j;
        std::vector<T> inputs(res.getRow(row));

        k_1j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RKF45_TIME(0, 1) : _RKF45_UNIT(1, 1));
             
        k_2j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RKF45_TIME(1, 2) : 
                ( _RKF45_UNIT(2, 1) + _RKF45_UNIT(2, 2) 
                - _RKF45_UNIT(1, 1)));
        
        k_3j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RKF45_TIME(2, 3) :
                ( _RKF45_UNIT(3, 1) + _RKF45_UNIT(3, 2) + _RKF45_UNIT(3, 3)
                - _RKF45_UNIT(2, 1) - _RKF45_UNIT(2, 2)));

        k_4j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RKF45_TIME(3, 4) : 
                ( _RKF45_UNIT(4, 1) + _RKF45_UNIT(4, 2) + _RKF45_UNIT(4, 3) + _RKF45_UNIT(4, 4)
                - _RKF45_UNIT(3, 1) - _RKF45_UNIT(3, 2) - _RKF45_UNIT(3, 3)));

        k_5j = ode._eval(inputs);

        _LOOP_TO_M(i, 0)
            inputs[i] += h * (i == 0 ? _RKF45_TIME(4, 5) : 
                ( _RKF45_UNIT(5, 1) + _RKF45_UNIT(5, 2) + _RKF45_UNIT(5, 3) + _RKF45_UNIT(5, 4) + _RKF45_UNIT(5, 5)
                - _RKF45_UNIT(4, 1) - _RKF45_UNIT(4, 2) - _RKF45_UNIT(4, 3) - _RKF45_UNIT(4, 4)));
        
        k_6j = ode._eval(inputs);

        std::vector<T> result(res.getRow(row)), w1(m + 1), w2(m + 1);
        T R = 0;    // calculate the error as the magnitude of the difference between w1 and w2

        _LOOP_TO_M(i, 1)
        {
            w1[i] = result[i] + h * (
                  _RKF45_UNIT(6, 1)
                + _RKF45_UNIT(6, 2)
                + _RKF45_UNIT(6, 3)
                + _RKF45_UNIT(6, 4)
                + _RKF45_UNIT(6, 5)
                + _RKF45_UNIT(6, 6)
            );
            
            w2[i] = result[i] + h * (
                  _RKF45_UNIT(7, 1)
                + _RKF45_UNIT(7, 2)
                + _RKF45_UNIT(7, 3)
                + _RKF45_UNIT(7, 4)
                + _RKF45_UNIT(7, 5)
                + _RKF45_UNIT(7, 6)
            );

            R += _ABS(w2[i] - w1[i]) * _ABS(w2[i] - w1[i]);
        }

        R = _RKF45_SQRT(R) / h;
        T delta = 0.84 * _RKF45_POW(maxError / R, 0.25);

        if (R <= maxError) 
        {
            result[0] += h;
            _LOOP_TO_M(i, 1) 
            {
                result[i] = w1[i];
            }
            res.addRow(result);

            h *= delta;
            t += h;
            row++;

        } else {
            h *= delta;
        }
    } 

    return res;
}



template <typename T>
DataFrame<T> _RKF45(ODESystem<T>& ode)
{
    return _RKF45(ode, (T)DEFAULT_MAX_ERROR);
} 



} // namespace DES


#endif