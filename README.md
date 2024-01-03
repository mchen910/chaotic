# Chaotic Motion #
This repository contains a small collection of simulations of chaotic motion, along with the underlying dynamics and governing equations. Chaotic systems are those in which small changes in the initial conditions give rise to completely different behavior as time progresses.

It also contains a small differential equation solver written in C++, with numerical methods of varying precision implemented. These include:

| Name | Method | 
| ------|------- |
|`euler` | Euler's method |
|`RK4`   | Runge Kutta Order 4 |
| `RK45` | Runge-Kutta-Fehlberg |
| `TSIT45` | Tsitouras
| `RB23` | Rosenbrock 


## Solver Specifics ##
The general process of solving a differential equation is as follows, with the example of a double pendulum. Defining $ \theta_1 $ and $ \theta_2$ as the angle of the first mass and the second mass, respectively, from the vertical yields the Lagrangian 

$$ \mathcal{L} = \frac{1}{2} (m_1 + m_2) l_1\dot{\theta_1}^2 + \frac{1}{2} m_2l_2^2\dot{\theta_2}^2 + m_2l_1l_2\dot{\theta_1}\dot{\theta_2}\cos(\theta_1-\theta_2) + (m_1 + m_2)gl_1\cos(\theta_1) + m_2gl_2\cos(\theta_2)$$

Applying the Euler-Lagrange equations, the two equations of motion are

$$ \ddot{\theta_1} + \ddot{\theta_2} \left( \frac{m_1}{m_1 + m_2} \frac{l_2}{l_1} \cos(\theta_1-\theta_2)\right) + \dot{\theta_2} \left( \frac{m_1}{m_1 + m_2} \frac{l_2}{l_1} \sin(\theta_1-\theta_2) \right) + \frac{g}{l_1} \sin(\theta_1) = 0 $$

$$ \ddot{\theta_2} + \ddot{\theta_1} \left( \frac{l_1}{l_2} \cos(\theta_1-\theta_2)\right) - \dot{\theta_1}^2\left(\frac{l_1}{l_2} \sin(\theta_1-\theta_2)\right) + \frac{g}{l_2}\sin(\theta_2) = 0 $$

These are second-order nonlinear differential equations, which the solver cannot handle directly. Instead, the equations must first be solved for $\ddot{\theta_1}$ and $\ddot{\theta_2}$ individually, then converted into a system of first-order equations. Doing so yields the following equations:

$$\omega_1 = \dot{\theta_1}$$
$$\omega_2 = \dot{\theta_2}$$
$$\dot{\omega_1} = \frac{−g (2m_1 + m_2)\sin\theta_1 − m_2g\sin(\theta_1 − 2\theta_2) − 2m_2\sin(\theta_1 − \theta_2)(\omega_2^2 l_2 + \omega_1^2 l_1 \cos(\theta_1 − \theta_2))}{l_1(2m_1 + m_2 - m_2\cos(2(\theta_1-\theta_2)))}$$
$$\dot{\omega_2} = \frac{2\sin(\theta_1-\theta_2) (\omega_1^2 l_1(m_1 + m_2) + g(m_1 + m_2)\cos(\theta_1) + \omega_2^2 l_2m_2\cos(\theta_1-\theta_2))}{l_2(2m_1 + m_2 - m_2\cos(2(\theta_1-\theta_2)))}$$

Next, `ODE` objects are created...