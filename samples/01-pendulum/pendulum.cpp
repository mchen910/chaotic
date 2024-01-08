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
#include "utils/line.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


#define G 0.30f
#define T float
#define N 1000


// Allow the user to pause and play the simulation
bool shouldPause = false;


class DoublePendulum
{

private:
    DES::ODESystem<T> system;
    Shader shader1, shader2;
    GLuint VAO1, VAO2;
    GLuint VBO1, VBO2, EBO1, EBO2;

    T l1, l2, m1, m2, theta1, theta2, omega1, omega2;
    T theta10, theta20, omega10, omega20;

    bool drawTrail1 = false;
    bool drawTrail2 = false;

    glm::vec3 color;
    float t = 0.005; // line thickness
    

    void init()
    {
        // functions are of the form f(t, θ₁, θ₂, ω₁, ω₂)
        DES::function_t<float> theta1prime ([&](std::vector<T> args) {
            return args[3];
        });

        DES::function_t<float> theta2prime ([&](std::vector<T> args) {
            return args[4];
        });

        DES::function_t<float> omega1prime ([&](std::vector<T> args) {
            float o1 = args[1], o2 = args[2];
            float w1 = args[3], w2 = args[4];

            float num = -G * (2 * m1 + m2) * sinf(o1) - m2 * G * sinf(o1 - 2 * o2) - 2 * sinf(o1 - o2) * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * cosf(o1 - o2));
            float den = l1 * (2 * m1 + m2 - m2 * cosf(2 * o1 - 2 * o2));

            return num / den;
        });

        DES::function_t<float> omega2prime ([&](std::vector<float> args) {
            float o1 = args[1], o2 = args[2];
            float w1 = args[3], w2 = args[4];

            float num = 2 * sinf(o1 - o2) * (w1 * w1 * l1 * (m1 + m2) + G * (m1 + m2) * cosf(o1) + w2 * w2 * l2 * m2 * cosf(o1 - o2));
            float den = l2 * (2 * m1 + m2 - m2 * cosf(2 * o1 - 2 * o2));

            return num / den;
        });

        DES::iv_t<T> initialConditions = { 0.0, theta1, theta2, omega1, omega2 };
        DES::timeBound_t<T> bounds = { 0.0, 0.0 };   // arbitrary, since the system is propagated through time indefinitely
        T timeStep = 0.1;                       

        system = DES::ODESystem<T>(initialConditions, {theta1prime, theta2prime, omega1prime, omega2prime}, bounds, timeStep);


        // Initialize lines (triangles to implement line thickness)
        shader1 = Shader("samples/shaders/01-vertex.glsl", "samples/shaders/01-fragment.glsl");
        shader2 = Shader("samples/shaders/01-vertex.glsl", "samples/shaders/01-fragment.glsl");

        float x1 = l1 * sinf(theta1), y1 = -l1 * cosf(theta1);
        float x2 = l1 * sinf(theta1) + l2 * sinf(theta2), y2 = -l1 * cosf(theta1) - l2 * cosf(theta2);

        float vertices1[] = {
                 t / 2 * cosf(theta1),      t / 2 * sinf(theta1), 0.0,
               - t / 2 * cosf(theta1),    - t / 2 * sinf(theta1), 0.0,
            x1 + t / 2 * cosf(theta1), y1 + t / 2 * sinf(theta1), 0.0,
            x1 - t / 2 * cosf(theta1), y1 - t / 2 * sinf(theta1), 0.0
        };

        unsigned int indices1[] = {
            0, 1, 3, 
            0, 2, 3
        };

        float vertices2[] = {
            x1 + t / 2 * cosf(theta2), y1 + t / 2 * sinf(theta2), 0.0,
            x1 - t / 2 * cosf(theta2), y1 - t / 2 * sinf(theta2), 0.0,
            x2 + t / 2 * cosf(theta2), y2 + t / 2 * sinf(theta2), 0.0,
            x2 - t / 2 * cosf(theta2), y2 - t / 2 * sinf(theta2), 0.0,
        };

        unsigned int indices2[] = {
            0, 1, 3,
            0, 2, 3
        };


        // line 1
        glGenBuffers(1, &VBO1);
        glGenBuffers(1, &EBO1);
        glGenVertexArrays(1, &VAO1);
        glBindVertexArray(VAO1);
        glBindBuffer(GL_ARRAY_BUFFER, VBO1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices1), vertices1, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO1);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices1), indices1, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        // line 2
        glGenBuffers(1, &VBO2);
        glGenBuffers(1, &EBO2);
        glGenVertexArrays(1, &VAO2);
        glBindVertexArray(VAO2);
        glBindBuffer(GL_ARRAY_BUFFER, VBO2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices2), vertices2, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices2), indices2, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        // unbind
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

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
        this->theta10 = theta1;
        this->theta20 = theta2;
        this->omega10 = omega1;
        this->omega20 = omega2;

        init();
    }


    ~DoublePendulum()
    {
        glDeleteVertexArrays(1, &VAO1);
        glDeleteVertexArrays(1, &VAO2);
        glDeleteBuffers(1, &VBO1);
        glDeleteBuffers(1, &EBO1);
        glDeleteBuffers(1, &VBO2);
        glDeleteBuffers(1, &EBO2);
    }


    void propagate()
    {
        std::vector<T> sol = DES::solve_i(system, ALGORITHM_RK4);

        // unpack (dataframe is of the form (t, θ₁, θ₂, ω₁, ω₂))
        theta1 = sol[1];
        theta2 = sol[2];
        omega1 = sol[3];
        omega2 = sol[4];

        // set the new positions via calculating transformation matrices and modifying uniform variables
        // line 1 is just a rotation about (0, 0)
        glm::mat4 matrix1 = glm::mat4(1.0f);
        glm::mat4 trans1 = glm::rotate(matrix1, theta1 - theta10, glm::vec3(0.0, 0.0, 1.0));
        glm::vec3 oldPivot = glm::vec3(l1 * sinf(theta10), -l1 * cosf(theta10), 0.0);
        glm::vec3 newPivot = glm::vec3(trans1 * glm::vec4(oldPivot, 1.0));

        // line 2 is a translation to the origin, a rotation, a translation to the old pivot, then a translation to the new pivot
        glm::mat4 matrix2 = glm::mat4(1.0f);
        glm::mat4 translation2 = glm::translate(matrix2, glm::vec3(l1 * sinf(theta10), -l1 * cos(theta10), 0.0));
        glm::mat4 rot2 = glm::rotate(matrix2, theta2 - theta20, glm::vec3(0.0, 0.0, 1.0));
        glm::mat4 translation2Inv = glm::inverse(translation2);
        glm::mat4 finalTranslation2 = glm::translate(matrix2, newPivot - oldPivot);
        glm::mat4 trans2 = finalTranslation2 * translation2 * rot2 * translation2Inv;

        shader1.use();
        GLuint transform1Loc = glGetUniformLocation(shader1.ID, "transform");
        glUniformMatrix4fv(transform1Loc, 1, GL_FALSE, glm::value_ptr(trans1));
        glBindVertexArray(VAO1);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        shader2.use();
        GLuint transform2Loc = glGetUniformLocation(shader2.ID, "transform");
        glUniformMatrix4fv(transform2Loc, 1, GL_FALSE, glm::value_ptr(trans2));
        glBindVertexArray(VAO2);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        glBindVertexArray(0);
    }


    void setColor(glm::vec4 color) 
    {
        shader1.use();
        GLuint colorLoc1 = glGetUniformLocation(shader1.ID, "color");
        glUniform4f(colorLoc1, color.x, color.y, color.z, color.w);

        shader2.use();
        GLuint colorLoc2 = glGetUniformLocation(shader2.ID, "color");
        glUniform4f(colorLoc2, color.x, color.y, color.z, color.w);
    }
};



void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_RELEASE) return; //only handle press events
    if(key == GLFW_KEY_SPACE) shouldPause = !shouldPause;
    else if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(window, true);
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

    GLFWwindow* window = glfwCreateWindow(600, 600, "01-Pendulum", NULL, NULL);
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

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);

    DoublePendulum* pendulums[N];
    glm::vec4 startColor(0.0, 0.0, 1.0, 1.0);
    glm::vec4 endColor(1.0, 0.0, 1.0, 1.0);

    for (int i = 0; i < N; i++) 
    {
        pendulums[i] = new DoublePendulum(0.5, 0.3, 1.0, 2.0, M_PI / 2 - 0.001 * i, 0.0, 0.0, 0.0);
        pendulums[i]->setColor(startColor + (endColor - startColor) * ((float)i / (N - 1)));
    }

    while (!glfwWindowShouldClose(window))
    {
        if (!shouldPause) {
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            
            for (int i = 0; i < N; i++) 
                pendulums[i]->propagate();
                
            glfwSwapBuffers(window);
        }
        glfwPollEvents();

    }

    glfwTerminate();

    return 0;
}


