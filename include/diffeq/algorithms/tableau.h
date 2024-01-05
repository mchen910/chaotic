#ifndef DIFFEQ_ALGORITHMS_TABLEAU_H
#define DIFFEQ_ALGORITHMS_TABLEAU_H


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
        {    TAB_NULL,      25.0/216.0,              0.0,    1408.0/2565.0,     2197.0/4104.0,     -1.0/5.0,        0.0    },
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
        {     TAB_NULLF,      25.0f/216.0f,               0.0f,    1408.0f/2565.0f,     2197.0f/4104.0f,     -1.0f/5.0f,         0.0f    },
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
        {     TAB_NULLL,      25.0l/216.0l,               0.0l,    1408.0l/2565.0l,     2197.0l/4104.0l,     -1.0l/5.0l,         0.0l    },
        {     TAB_NULLL,      16.0l/135.0l,               0.0l,   6656.0l/12825.0l,   28561.0l/56430.0l,    -9.0l/50.0l,   2.0l/55.0l    }
    };

};


}


#endif