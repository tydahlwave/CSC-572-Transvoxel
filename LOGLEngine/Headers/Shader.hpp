//
//  Shader.hpp
//  LOGLEngine
//
//  Created by Tyler Dahl on 5/10/17.
//
//

#ifndef Shader_h
#define Shader_h

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

// Include glad to get all the required OpenGL headers
#include <glad/glad.h>

class Shader {
public:
    // The program ID
    GLuint program;
    
    // Constructor reads and builds the shader
    Shader(std::string vertexPath, std::string fragmentPath) {
        // Retrieve the shader source code from filePaths
        std::string vShaderCode = getShaderCode(vertexPath);
        std::string fShaderCode = getShaderCode(fragmentPath);

        // Compile shaders
        GLuint vertex = compileShader(GL_VERTEX_SHADER, vShaderCode);
        GLuint fragment = compileShader(GL_FRAGMENT_SHADER, fShaderCode);

        // Create and link shader program
        linkProgram(vertex, fragment);
    }
    
    // Use the program
    void use() { glUseProgram(this->program); }
    
    // Get the location of a uniform
    GLint getUniform(std::string name) {
        GLint loc = glGetUniformLocation(program, name.c_str());
        if (loc < 0) {
            std::cerr << "WARN: "<< name << " cannot be bound (it either doesn't exist or has been optimized away). safe_glAttrib calls will silently ignore it.\n" << std::endl;
        }
        return loc;
    }
    
    // Get the location of an attribute
    GLint getAttribute(std::string name) {
        GLint loc = glGetAttribLocation(program, name.c_str());
        if (loc < 0) {
            std::cerr << "WARN: "<< name << " cannot be bound (it either doesn't exist or has been optimized away). safe_glAttrib calls will silently ignore it.\n" << std::endl;
        }
        return loc;
    }
private:
    std::string getShaderCode(std::string filePath) {
        std::string shaderCode;
        std::ifstream shaderFile;

        // Ensures ifstream object can throw exceptions
        shaderFile.exceptions (std::ifstream::badbit);
        
        try {
            // Open file
            shaderFile.open(filePath);
            
            std::stringstream shaderStream;
            // Read file's buffer content into stream
            shaderStream << shaderFile.rdbuf();
            // Close file handler
            shaderFile.close();
            // Convert stream into string
            shaderCode = shaderStream.str();
        } catch (std::ifstream::failure e) {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
        }

        return shaderCode;
    }
    
    GLuint compileShader(GLenum shaderType, std::string shaderCode) {
        GLuint shader;
        GLint success;
        GLchar infoLog[512];
        
        // Necessary conversion from string to const GLchar * (see glShaderSource method)
        const GLchar *code = shaderCode.c_str();
        
        // Compile shader
        shader = glCreateShader(shaderType);
        glShaderSource(shader, 1, &code, NULL);
        glCompileShader(shader);
        
        // Print compile errors if any
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(shader, 512, NULL, infoLog);
            std::cout << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
        }
        return shader;
    }
    
    void linkProgram(GLuint vertexShader, GLuint fragmentShader) {
        GLint success;
        GLchar infoLog[512];
        
        // Create and link shader program
        this->program = glCreateProgram();
        glAttachShader(this->program, vertexShader);
        glAttachShader(this->program, fragmentShader);
        glLinkProgram(this->program);
        
        // Print linking errors if any
        glGetProgramiv(this->program, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(this->program, 512, NULL, infoLog);
            std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        }
        
        // Delete the shaders as they're linked into our program now and no longer necessary
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }
};

#endif /* Shader_h */
