/**
 * @file diffeq.h
 * @author Matthew Chen
 * @brief A small header-only library to numerically solve differential 
 * equations. 
 * @version 0.1
 * @date 2024-01-03
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef DIFFEQ_H
#define DIFFEQ_H


#if !defined(DIFFEQ_FLOAT_PRECISION) &&     \
    !defined(DIFFEQ_DOUBLE_PRECISION) &&    \
    !defined(DIFFEQ_LONG_DOUBLE_PRECISION)  
#error Must specify precision
#endif


#define ALGORITHM_EULER     0x001
#define ALGORITHM_RK4       0x002
#define ALGORITHM_RKF45     0x003


#include "diffeq/dataframe.h"
#include "diffeq/ode.h"
#include "diffeq/algorithms/rk.h"


#endif