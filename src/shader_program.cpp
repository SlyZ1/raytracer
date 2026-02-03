#include "shader_program.hpp"

ShaderProgram::ShaderProgram(){
    shaderProgram = glCreateProgram();
}

string ShaderProgram::getShaderSource(const char* path){
    ifstream file(path);

    if (!file.is_open()) {
        cerr << "Erreur: impossible d'ouvrir le fichier " << path << endl;
        return "";
    }

    stringstream shaderText;
    shaderText << file.rdbuf(); 
    return shaderText.str();
}

void ShaderProgram::load(int type, const char *path){
    //Create shader
    unsigned int shader;
    shader = glCreateShader(type);
    shaders.push_back(shader);
    
    //Compile shader
    string shaderSourceString = getShaderSource(path); 
    const char *shaderSource = shaderSourceString.c_str();
    glShaderSource(shader, 1, &shaderSource, NULL);
    glCompileShader(shader);

    //Check for errors
    int success;
    char infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(shader, 512, NULL, infoLog);
        cerr << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << endl;
    }

    //Add to program
    glAttachShader(shaderProgram, shader);
}

void ShaderProgram::link(){
    //Link
    glLinkProgram(shaderProgram);

    //Check for errors
    int success;
    char infoLog[512];
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        cerr << "ShaderProgram linking failed : " << infoLog << endl;
    }
    
    //Delete shaders
    for(const int shader : shaders){
        glDeleteShader(shader);
    }
}

void ShaderProgram::use(){
    glUseProgram(shaderProgram);
}

void ShaderProgram::destroy(){
    glDeleteProgram(shaderProgram);
}

unsigned int ShaderProgram::id(){
    return shaderProgram;
}