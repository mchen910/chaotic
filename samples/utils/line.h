#ifndef UTILS_LINE_H
#define UTILS_LINE_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>

#include <vector>

#include "shader.h"


class Line
{

private:
    Shader shaderProgram;
    GLuint VAO, VBO;

    glm::vec3 start, end;
    glm::mat4 MVP;
    glm::vec3 color;

    std::vector<float> vertices;

public:
    Line() = default;
    Line(glm::vec3 start, glm::vec3 end, Shader& shader);
    ~Line();

    void draw();
    void setMVP(glm::mat4 mvp);
    void setColor(glm::vec3 color);
    void setVertices(glm::vec3 start, glm::vec3 end);
};





Line::Line(glm::vec3 start, glm::vec3 end, Shader& shader) : shaderProgram(shader)
{
    this->start = start;
    this->end = end;
    
    this->color = glm::vec3(1.0, 1.0, 1.0);
    this->MVP = glm::mat4(1.0f);

    this->vertices = {
        start.x,    start.y,    start.z,
        end.x,      end.y,      end.z
    };

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

}

Line::~Line()
{
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram.ID);
}


void Line::draw()
{
    glUseProgram(shaderProgram.ID);
    // glUniformMatrix4fv(glGetUniformLocation(shaderProgram.ID, "MVP"), 1, GL_FALSE, &MVP[0][0]);
    // glUniform3fv(glGetUniformLocation(shaderProgram.ID, "color"), 1, &color[0]);

    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, 2);
}


void Line::setColor(glm::vec3 color) 
{
    this->color = color;
}


void Line::setMVP(glm::mat4 mvp) 
{
    this->MVP = mvp;
}


void Line::setVertices(glm::vec3 start, glm::vec3 end)
{
    this->start = start;
    this->end = end;

    // rebind the VAO/VBO
    this->vertices = {
        start.x,    start.y,    start.z,
        end.x,      end.y,      end.z
    };

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}


#endif