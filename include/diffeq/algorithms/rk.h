#ifndef DIFFEQ_ALGORITHMS_RK_H
#define DIFFEQ_ALGORITHMS_RK_H

#include <cmath>

#include "../ode.h"


namespace DES 
{

/**
 * @brief The Butcher Tableau contains the weights and nodes used in 
 * Runge-Kutta numerical methods. The first column contains the weights 
 * for the time variable, the last row of the table contains the weights 
 * for the values of k computed, and the body of the table are the coefficients 
 * used for calculating k.
 * 
 * For adaptive methods (e.g. RKF45), there is a second row attached that is used 
 * to compute the error at each step, in order to modify the timestep.
 * 
 * Each tableau has a float and double version, depending on memory/precision needs.
 * 
 */
namespace ButcherTableau
{
    // Butcher Tableaus have empty entries, use a macro for readability
    #define TAB_NULL    0.0
    #define TAB_NULLF   0.0f
    #define TAB_NULLL   0.0l



    const double EULER_TAB[2][2] 
    {
        {        0.0,   0.0   },
        {   TAB_NULL,   1.0   }
    };


    const float EULER_TABF[2][2] 
    {
        {        0.0f,   0.0f   },
        {   TAB_NULLF,   1.0f   }
    };


    const long double EULER_TABL[2][2]
    {
        {        0.0l,   0.0l   },
        {   TAB_NULLL,   1.0l   }
    };



    const double RK4_TAB[5][5] 
    {
        {        0.0,   TAB_NULL,   TAB_NULL,   TAB_NULL,   TAB_NULL   },
        {    1.0/2.0,    1.0/2.0,   TAB_NULL,   TAB_NULL,   TAB_NULL   },
        {    1.0/2.0,        0.0,    1.0/2.0,   TAB_NULL,   TAB_NULL   },
        {        1.0,        0.0,        0.0,        1.0,   TAB_NULL   },
        {   TAB_NULL,    1.0/6.0,    1.0/3.0,    1.0/3.0,    1.0/6.0   }
    };


    const float RK4_TABF[5][5] 
    {
        {        0.0f,   TAB_NULLF,   TAB_NULLF,   TAB_NULLF,   TAB_NULLF   },
        {   1.0f/2.0f,   1.0f/2.0f,   TAB_NULLF,   TAB_NULLF,   TAB_NULLF   },
        {   1.0f/2.0f,        0.0f,   1.0f/2.0f,   TAB_NULLF,   TAB_NULLF   },
        {        1.0f,        0.0f,        0.0f,        1.0f,   TAB_NULLF   },
        {   TAB_NULLF,   1.0f/6.0f,   1.0f/3.0f,   1.0f/3.0f,   1.0f/6.0f   }
    };


    const long double RK4_TABL[5][5]
    {
        {        0.0l,   TAB_NULLL,   TAB_NULLL,   TAB_NULLL,   TAB_NULLL   },
        {   1.0l/2.0l,   1.0l/2.0l,   TAB_NULLL,   TAB_NULLL,   TAB_NULLL   },
        {   1.0l/2.0l,        0.0l,   1.0l/2.0l,   TAB_NULLL,   TAB_NULLL   },
        {        1.0l,        0.0l,        0.0l,        1.0l,   TAB_NULLL   },
        {   TAB_NULLL,   1.0l/6.0l,   1.0l/3.0l,   1.0l/3.0l,   1.0l/6.0l   }
    };


    const double RK38_TAB[5][5] 
    {
        {        0.0,   TAB_NULL,   TAB_NULL,   TAB_NULL,   TAB_NULL   },
        {    1.0/3.0,    1.0/3.0,   TAB_NULL,   TAB_NULL,   TAB_NULL   },
        {    2.0/3.0,   -1.0/3.0,        1.0,   TAB_NULL,   TAB_NULL   },
        {        1.0,        1.0,       -1.0,        1.0,   TAB_NULL   },
        {   TAB_NULL,    1.0/8.0,    3.0/8.0,    3.0/8.0,    3.0/8.0   }
    };


    const float RK38_TABF[5][5] 
    {
        {        0.0f,   TAB_NULLF,   TAB_NULLF,   TAB_NULLF,   TAB_NULLF   },
        {   1.0f/3.0f,   1.0f/3.0f,   TAB_NULLF,   TAB_NULLF,   TAB_NULLF   },
        {   2.0f/3.0f,   1.0f/3.0f,        1.0f,   TAB_NULLF,   TAB_NULLF   },
        {        1.0f,        1.0f,       -1.0f,        1.0f,   TAB_NULLF   },
        {   TAB_NULLF,   1.0f/8.0f,   3.0f/8.0f,   3.0f/8.0f,   3.0f/8.0f   }
    };


    const long double RK38_TABL[5][5] 
    {
        {        0.0l,   TAB_NULLL,   TAB_NULLL,   TAB_NULLL,   TAB_NULLL   },
        {   1.0l/3.0l,   1.0l/3.0l,   TAB_NULLL,   TAB_NULLL,   TAB_NULLL   },
        {   2.0l/3.0l,   1.0l/3.0l,        1.0l,   TAB_NULLL,   TAB_NULLL   },
        {        1.0l,        1.0l,       -1.0l,        1.0l,   TAB_NULLL   },
        {   TAB_NULLL,   1.0l/8.0l,   3.0l/8.0l,   3.0l/8.0l,   3.0l/8.0l   }
    };


    const double RKF45_TAB[8][7] 
    {
        {         0.0,        TAB_NULL,         TAB_NULL,         TAB_NULL,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {     1.0/4.0,         1.0/4.0,         TAB_NULL,         TAB_NULL,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {     3.0/8.0,        3.0/32.0,         9.0/32.0,         TAB_NULL,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {   12.0/13.0,   1932.0/2197.0,   -7200.0/2197.0,    7296.0/2197.0,          TAB_NULL,     TAB_NULL,   TAB_NULL    },
        {         1.0,     439.0/216.0,             -8.0,     3680.0/513.0,     -845.0/4104.0,     TAB_NULL,   TAB_NULL    },
        {     1.0/2.0,       -8.0/27.0,              2.0,   -3544.0/2565.0,     1859.0/4104.0,   -11.0/40.0,   TAB_NULL    },
        {    TAB_NULL,      25.0/216.0,              0.0,    1408.0/2565.0,     2197.0/1404.0,     -1.0/5.0,        0.0    },
        {    TAB_NULL,      16.0/135.0,              0.0,   6656.0/12825.0,   28561.0/56430.0,    -9.0/50.0,   2.0/55.0    }
    };


    const float RKF45_TABF[8][7] 
    {
        {          0.0f,         TAB_NULLF,          TAB_NULLF,          TAB_NULLF,           TAB_NULLF,      TAB_NULLF,    TAB_NULLF    },
        {     1.0f/4.0f,         1.0f/4.0f,          TAB_NULLF,          TAB_NULLF,           TAB_NULLF,      TAB_NULLF,    TAB_NULLF    },
        {     3.0f/8.0f,        3.0f/32.0f,         9.0f/32.0f,          TAB_NULLF,           TAB_NULLF,      TAB_NULLF,    TAB_NULLF    },
        {   12.0f/13.0f,   1932.0f/2197.0f,   -7200.0f/2197.0f,    7296.0f/2197.0f,           TAB_NULLF,      TAB_NULLF,    TAB_NULLF    },
        {          1.0f,     439.0f/216.0f,              -8.0f,     3680.0f/513.0f,     -845.0f/4104.0f,      TAB_NULLF,    TAB_NULLF    },
        {     1.0f/2.0f,       -8.0f/27.0f,               2.0f,   -3544.0f/2565.0f,     1859.0f/4104.0f,   -11.0f/40.0f,    TAB_NULLF    },
        {     TAB_NULLF,      25.0f/216.0f,               0.0f,    1408.0f/2565.0f,     2197.0f/1404.0f,     -1.0f/5.0f,         0.0f    },
        {     TAB_NULLF,      16.0f/135.0f,               0.0f,   6656.0f/12825.0f,   28561.0f/56430.0f,    -9.0f/50.0f,   2.0f/55.0f    }
    };


    const long double RKF45_TABL[8][7] 
    {
        {          0.0l,         TAB_NULLL,          TAB_NULLL,          TAB_NULLL,           TAB_NULLL,      TAB_NULLL,    TAB_NULLL    },
        {     1.0l/4.0l,         1.0l/4.0l,          TAB_NULLL,          TAB_NULLL,           TAB_NULLL,      TAB_NULLL,    TAB_NULLL    },
        {     3.0l/8.0l,        3.0l/32.0l,         9.0l/32.0l,          TAB_NULLL,           TAB_NULLL,      TAB_NULLL,    TAB_NULLL    },
        {   12.0l/13.0l,   1932.0l/2197.0l,   -7200.0l/2197.0l,    7296.0l/2197.0l,           TAB_NULLL,      TAB_NULLL,    TAB_NULLL    },
        {          1.0l,     439.0l/216.0l,              -8.0l,     3680.0l/513.0l,     -845.0l/4104.0l,      TAB_NULLL,    TAB_NULLL    },
        {     1.0l/2.0l,       -8.0l/27.0l,               2.0l,   -3544.0l/2565.0l,     1859.0l/4104.0l,   -11.0l/40.0l,    TAB_NULLL    },
        {     TAB_NULLL,      25.0l/216.0l,               0.0l,    1408.0l/2565.0l,     2197.0l/1404.0l,     -1.0l/5.0l,         0.0l    },
        {     TAB_NULLL,      16.0l/135.0l,               0.0l,   6656.0l/12825.0l,   28561.0l/56430.0l,    -9.0l/50.0l,   2.0l/55.0l    }
    };

};


#define  DEFAULT_ERROR          0.01
#define  _LOOP_TO_M(i, n)       for (size_t i = n; i <= m; i++)



template <typename T>
DataFrame<T> _EULER(ODESystem<T> & ode) 
{
   timeBound_t<T> tBound = ode.getTimeBound();

    size_t m = ode.getNumEquations();
    size_t row = 0;

    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1);
    res.addRow(ode._iValues.vec);


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
DataFrame<T> _RK4(ODESystem<T>& ode) 
{
    timeBound_t<T> tBound = ode.getTimeBound();

    size_t m = ode.getNumEquations();
    size_t row = 0;

    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1);
    res.addRow(ode._iValues.vec);


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
DataFrame<T> _RKF45(ODESystem<T>& ode, T maxError) 
{
    timeBound_t<T> tBound = ode.getTimeBound();

    size_t m = ode.getNumEquations();
    size_t row = 0;

    T h = ode.getTimeStep();
    T t = tBound.first;

    DataFrame<T> res(0, m + 1);
    res.addRow(ode._iValues.vec);


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

    std::cout << "here" << std::endl;

    using namespace ButcherTableau;
    do
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

        std::vector<T> truncErr, result(res.getRow(row));
        result[0] = t + h;

        _LOOP_TO_M(i, 1)
        {
            result[i] += h * (
                  _RKF45_UNIT(6, 1)
                + _RKF45_UNIT(6, 2)
                + _RKF45_UNIT(6, 3)
                + _RKF45_UNIT(6, 4)
                + _RKF45_UNIT(6, 5)
                + _RKF45_UNIT(6, 6)
            );

            truncErr.push_back( 
                  _RKF45_UNIT(7, 1)
                + _RKF45_UNIT(7, 2)
                + _RKF45_UNIT(7, 3)
                + _RKF45_UNIT(7, 4)
                + _RKF45_UNIT(7, 5)
                + _RKF45_UNIT(7, 6)
            );
        }

        res.addRow(result);
        std::cout << "added result" << std::endl;

        T truncErrMagnitude = _RKF45_SQRT(
              truncErr[0] * truncErr[0] 
            + truncErr[1] * truncErr[1] 
            + truncErr[2] * truncErr[2] 
            + truncErr[3] * truncErr[3] 
            + truncErr[4] * truncErr[4] 
            + truncErr[5] * truncErr[5] 
        );

        h = (T)0.9 * h * _RKF45_POW(maxError / truncErrMagnitude, 0.2);
        std::cout << "h: " << h << std::endl;
        row++;
        if (truncErrMagnitude <= maxError) {
            t += h;
        }
        // break;
    } while (t < tBound.second);

    return res;
}



template <typename T>
DataFrame<T> _RKF45(ODESystem<T>& ode)
{
    return _RKF45(ode, (T)DEFAULT_ERROR);
} 


} // namespace DES


#endif