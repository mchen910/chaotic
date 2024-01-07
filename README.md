# Chaotic Motion #
This repository contains a small collection of simulations of chaotic motion, along with the underlying dynamics and governing equations. Chaotic systems are those in which small changes in the initial conditions give rise to completely different behavior as time progresses.

It also contains a small differential equation solver written in C++, with numerical methods of varying precision implemented. These include:

| Name      | Method                | Notes |
| --------- | --------------------- | ------|
| `euler`   | Euler's Method        | Shouldn't be used in most cases due to low precision. |
| `RK4`     | Runge-Kutta Order 4   | Canonical numerical method with a fixed timestep.
| `RKF45`   | Runge-Kutta-Fehlberg  | Adaptive timestep, so additional interpolation is required for fixed-step computation.
| `TSIT45`  | Tsitouras             | Should be used in most cases to solve non-stiff systems.
| `RB23`    | Rosenbrock            | Used for stiff systems.

The precision may also be specified, by defining one of the following macros *before* importing `diffeq.h`.

| Name      | Precision                | 
| --------- | --------------------- |
| `DIFFEQ_FLOAT_PRECISION`   | Float (7 digits)        |
| `DIFFEQ_DOUBLE_PRECISION`     | Double (15 digits)   |
| `DIFFEQ_LONG_DOUBLE_PRECISION`   | Long Double (17 digits, platform dependent)  |


## Solver Specifics ##
The general process of solving a differential equation is as follows, with the example of a simple harmonic oscillator with frequency $` \omega \equiv 1 `$ for simplicity, and initial conditions $`x(0) \equiv 1`$ and $`\dot{x_0}(0) \equiv 10 `$.


### 1. Setting up the equations ###
The differential equation describing the oscillator's motion is 

$$ \ddot{x} + x = 0 $$

, a second order ordinary differential equation which the solver cannot handle directly. Therefore, it is necessary to break up the original equation into two first order equations to interface with the solver. These are:

$$ \dot{v} = -u $$

$$ \dot{u} = v $$

where the substitution $`v = \dot{x}`$ and $`u = x`$ where made. 


### 2. Create `function_t` functors ###
Next, create the functions required to represent this system, using the type `function_t`. The order of the arguments is especially important - the first argument must be time, and the rest must be in the same order as the initial conditions provided.

### 3. Solve the system
The resulting `function_t` objects can now be passed into an `ODESystem` object and solved, given the time bounds and timestep. There is also an option to step through one timestep only and solve the system interatively through time, which is useful for simulations in real time. Below is the complete example.


```cpp
#define DIFFEQ_FLOAT_PRECISION
#include "diffeq.h"
#define T float     // match with the precision specified above, for convenience

using namespace DES;

int main(int argc, char** argv) 
{
    // The initial conditions are in the order (t, u, v), so the arguments to the functions are also in this order
    iv_t<T> initialConditions = {0.0, 1.0, 10.0};   

    // u' = v
    function_t<T> uPrime([](std::vector<T> args) {
        return args[2];
    });

    // v' = -u
    function_t<T> vPrime([](std::vector<T> args) {
        return -args[1];
    });

    // Time bounds
    timeBound_t<T> bounds = { 0.0, 1.0 };

    // Timestep
    T dT = 0.1;

    // Create the ODESystem
    ODESystem<T> system(initialConditions, { uPrime, vPrime }, bounds, dT);

    // Solve with a given algorithm
    DataFrame<T> sol = solve(system, ALGORITHM_RK4);

    return 0;
}
```


The exact solution of the differential equation above is 

$$ x(t)=\cos(x)+10\sin(x) $$

A comparison between the solver and the exact solution is below.

| t | u | u (exact) | v | v (exact) |
|--|--|--|--|--|
| 0   | 1        | 1           | 10       | 10          |
| 0.1 | 1.993338 | 1.993338332 | 9.850208 | 9.850208236 |
| 0.2 | 2.966758 | 2.966759886 | 9.601996 | 9.601996448 |
| 0.3 | 3.910536 | 3.910538556 | 9.257845 | 9.257844685 |
| 0.4 | 4.815241 | 4.815244417 | 8.821193 | 8.821191598 |
| 0.5 | 5.671834 | 5.671837948 | 8.296402 | 8.296400080 |
| 0.6 | 6.471756 | 6.471760349 | 7.688716 | 7.688713676 |
| 0.7 | 7.207015 | 7.207019060 | 7.004208 | 7.004204186 |
| 0.8 | 7.870264 | 7.870267618 | 6.249715 | 6.249711003 |
| 0.9 | 8.454875 | 8..45487906 | 5.432778 | 5.432772773 |
| 1   | 8.955009 | 8.955012154 | 4.561558 | 4.561552074 |

