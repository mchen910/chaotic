#ifndef UTILS_SHADER_H
#define UTILS_SHADER_H

#include <GLFW/glfw3.h>

#include <fstream>
#include <string>
#include <sstream>


class Shader
{

private:
    GLuint compileShader(GLuint type, const std::string& source);
    std::string readFile(const std::string& filepath);
    void createShader(const std::string& vertexShader, const std::string& fragmentShader);


public:
    GLuint ID;
    
    Shader() = default;
    Shader(const std::string& vertexShader, const std::string& fragmentShader);
    ~Shader();

};


Shader::Shader(const std::string& vertexShader, const std::string& fragmentShader)
{
    std::string vShader = readFile(vertexShader);
    std::string fShader = readFile(fragmentShader);
    createShader(vShader, fShader);
}



Shader::~Shader()
{
    if (ID != 0) 
        glDeleteProgram(ID);
}


std::string Shader::readFile(const std::string& filePath) 
{
    std::ifstream shaderFile;
    shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try
    {
        shaderFile.open(filePath);
        std::stringstream shaderStream;

        shaderStream << shaderFile.rdbuf();
        shaderFile.close();

        std::string shaderCode = shaderStream.str(); 
        return shaderCode;
    
    } catch(const std::ifstream::failure &e)
    {
        throw std::runtime_error("Failed to read shader");
    }
    
}


GLuint Shader::compileShader(GLuint type, const std::string& source)
{
    GLuint hShader = glCreateShader(type);

    const char* src = source.c_str();
    glShaderSource(hShader, 1, &src, NULL);
    glCompileShader(hShader);

    GLint result;
    glGetShaderiv(hShader, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE)
    {
        char infoLog[512];
        glGetShaderInfoLog(hShader, 512, NULL, infoLog);
        #include <iostream>
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
                glDeleteShader(hShader);

        throw std::runtime_error("Failed to compile shader " + std::string(infoLog));
    }

    return hShader;
}


void Shader::createShader(const std::string& vertexShader, const std::string& fragmentShader)
{
    GLuint vs = compileShader(GL_VERTEX_SHADER, vertexShader);
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, fragmentShader);

    ID = glCreateProgram();

    glAttachShader(ID, vs);
    glAttachShader(ID, fs);

    glLinkProgram(ID);

    GLint result;
    glGetShaderiv(ID, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE)
    {
        char infoLog[512];
        glGetProgramInfoLog(ID, 512, NULL, infoLog);
        #include <iostream>
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        glDeleteProgram(ID);
        throw std::runtime_error("Failed to create shader program");
    }

    glValidateProgram(ID);

    glDeleteShader(vs);
    glDeleteShader(fs);

}


#endif