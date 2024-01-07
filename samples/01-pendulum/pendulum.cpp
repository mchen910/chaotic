/**
 * @file pendulum.cpp
 * @author Matthew Chen
 * @brief A simple simulation of a double pendulum. See samples/pendulum/README.md 
 * for details.
 * @version 0.1
 * @date 2024-01-05
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#define DIFFEQ_FLOAT_PRECISION
#include "diffeq.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/vec3.hpp>

#include <utils/line.h>


#if defined(DIFFEQ_FLOAT_PRECISION)
    #define SIN(x)      sinf(x)
    #define COS(x)      cosf(x)
#elif defined(DIFFEQ_DOUBLE_PRECISION)
    #define SIN(x)      sin(x)
    #define COS(x)      cos(x)
#elif defined(DIFFEQ_LONG_DOUBLE_PRECISION)
    #define SIN(x)      sinl(x)
    #define COS(x)      cosl(x)
#endif

#define G 9.81f
#define T float


// Allow the user to pause and play the simulation
bool shouldPause = false;


class DoublePendulum
{

private:
    DES::ODESystem<T> system;
    Line line1, line2;

    T l1, l2, m1, m2, theta1, theta2, omega1, omega2;

    bool drawTrail1 = false;
    bool drawTrail2 = false;

    glm::vec3 color;


    void init()
    {
        // functions are of the form f(t, θ₁, θ₂, ω₁, ω₂)
        DES::function_t<T> theta1prime ([&](std::vector<T> args) {
            return args[3];
        });

        DES::function_t<T> theta2prime ([&](std::vector<T> args) {
            return args[4];
        });

        DES::function_t<T> omega1prime ([&](std::vector<T> args) {
            T o1 = args[1], o2 = args[2];
            T w1 = args[3], w2 = args[4];

            T num = -G * (2 * m1 + m2) * SIN(o1) - m2 * G * SIN(o1 - 2 * o2) - 2 * SIN(o1 - o2) * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * COS(o1 - o2));
            T den = l1 * (2 * m1 + m2 - m2 * COS(2 * o1 - 2 * o2));

            return num / den;
        });

        DES::function_t<T> omega2prime ([&](std::vector<T> args) {
            T o1 = args[1], o2 = args[2];
            T w1 = args[3], w2 = args[4];

            T num = 2 * SIN(o1 - o2) * (w1 * w1 * l1 * (m1 + m2) + G * (m1 + m2) * COS(o1) + w2 * w2 * l2 * m2 * COS(o1 - o2));
            T den = l2 * (2 * m1 + m2 - m2 * COS(2 * o1 - 2 * o2));

            return num / den;
        });

        DES::iv_t<T> initialConditions = { theta1, theta2, omega1, omega2 };
        DES::timeBound_t<T> bounds = { 0.0, 0.0 };   // arbitrary, since the system is propagated through time indefinitely
        T timeStep = 0.1;                       

        system = DES::ODESystem<T>(initialConditions, {theta1prime, theta2prime, omega1prime, omega2prime}, bounds, timeStep);


        // Initialize lines
        Shader shader1("samples/shaders/01-vertex.glsl", "samples/shaders/01-fragment.glsl"), 
               shader2("samples/shaders/01-vertex.glsl", "samples/shaders/01-fragment.glsl");

        line1 = Line(
            glm::vec3(0.0, 0.0, 0.0),
            glm::vec3(-l1 * SIN(theta1), l1 * COS(theta1), 0.0),
            shader1
        );

        line2 = Line(
            glm::vec3(-l1 * SIN(theta1), l1 * COS(theta1), 0.0),
            glm::vec3(-l2 * SIN(theta2) - l1 * SIN(theta1), l1 * COS(theta1) + l2 * COS(theta2), 0.0),
            shader2
        );

        line1.draw();
        line2.draw();
    }


public:
    DoublePendulum(T l1, T l2, T m1, T m2, T theta1, T theta2, T omega1, T omega2)
    {
        this->l1 = l1;
        this->l2 = l2;
        this->m1 = m1;
        this->m2 = m2;
        this->theta1 = theta1;
        this->theta2 = theta2;
        this->omega1 = omega1;
        this->omega2 = omega2;

        init();       

    }

    void propagate()
    {
        std::vector<T> sol = DES::solve_i(system, ALGORITHM_RK4);
        std::cout << "solved" << std::endl;
        // unpack (dataframe is of the form (t, θ₁, θ₂, ω₁, ω₂))
        theta1 = sol[1];
        theta2 = sol[2];
        omega1 = sol[3];
        omega2 = sol[4];

        std::cout << theta1 << ", " << theta2 << ", " << omega1 << ", " << omega2 << std::endl;

        // set the new positions
        line1.setVertices(
            glm::vec3(0.0, 0.0, 0.0),
            glm::vec3(-l1 * SIN(theta1), l1 * COS(theta1), 0.0)
        );

        line2.setVertices(
            glm::vec3(-l1 * SIN(theta1), l1 * COS(theta1), 0.0),
            glm::vec3(-l2 * SIN(theta2) - l1 * SIN(theta1), l1 * COS(theta1) + l2 * COS(theta2), 0.0)
        );

        line1.draw();
        line2.draw();
    }
};



void framebufferSizeCallback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    else if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
        shouldPause = !shouldPause;
}



int main(int argc, char const *argv[])
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(800, 600, "01-Pendulum", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSetKeyCallback(window, keyCallback);


    DoublePendulum p1(0.3, 0.3, 1.0, 1.0, 0.3, 0.0, 0.0, 0.0);
    size_t i = 0;

    while (!glfwWindowShouldClose(window))
    {
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        // std::cout << "here" << std::endl;
        // p1.propagate();

        Shader shader("samples/shaders/01-vertex.glsl", "samples/shaders/01-fragment.glsl");

        Line l(glm::vec3(0, 0, 0), glm::vec3(0, 1, 0), shader);
        l.draw();

        if (i == 0) glfwSwapBuffers(window);
        glfwPollEvents();
        i++;

    }

    glfwTerminate();

    return 0;
}


