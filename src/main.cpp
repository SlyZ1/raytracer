#include <iostream>
#include <fstream>
#include <sstream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include "app.hpp"
#include "shader_program.hpp"
#include "camera.hpp"

using namespace std;

#define SAMPLES 10

string getShaderSource(const char *filepath){
    ifstream file(filepath);

    if (!file.is_open()) {
        cerr << "Erreur: impossible d'ouvrir le fichier " << filepath << endl;
        return "";
    }

    stringstream shaderText;
    shaderText << file.rdbuf(); 
    return shaderText.str();
}


App app;
int main(){
    app.init(800, 600, "GLSL Project");
    app.toggleCursor(false);

    ShaderProgram rayTraceShader = ShaderProgram();
    rayTraceShader.load(GL_VERTEX_SHADER, "src/shaders/vertex.glsl");
    rayTraceShader.load(GL_FRAGMENT_SHADER, "src/shaders/frag.glsl");
    rayTraceShader.link();

    ShaderProgram accumulationShader = ShaderProgram();
    accumulationShader.load(GL_VERTEX_SHADER, "src/shaders/screenVertex.glsl");
    accumulationShader.load(GL_FRAGMENT_SHADER, "src/shaders/screenFrag.glsl");
    accumulationShader.link();

    Camera camera(0.1, 0.3);
    int camPosLoc = glGetUniformLocation(rayTraceShader.id(), "camera.pos");
    int camDirLoc = glGetUniformLocation(rayTraceShader.id(), "camera.lookDir");
    int fragSamplesLoc = glGetUniformLocation(rayTraceShader.id(), "samples");
    
    vector<float> vertices = {
        1.f,  1.f, 0.0f,
        1.f, -1.f, 0.0f,
        -1.f, -1.f, 0.0f,
        -1.f,  1.f, 0.0f
    };
    vector<unsigned int> indices = {
        0, 1, 3,
        1, 2, 3
    };  
    unsigned int VBO, VAO, EBO;
    tie(VBO, VAO, EBO) = ShaderProgram::addData(vertices, indices);
    ShaderProgram::linkData(3, sizeof(float), 0);
    
    unsigned int FBO;
    glGenFramebuffers(1, &FBO);
    
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, app.width(), app.height(), 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    unsigned int oldTexture;
    glGenTextures(1, &oldTexture);
    glBindTexture(GL_TEXTURE_2D, oldTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, app.width(), app.height(), 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    glDisable(GL_FRAMEBUFFER_SRGB);
    int frameCount = SAMPLES;
    int samples = SAMPLES;
    while(!app.shouldClose())
    {
        //Current frame
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);
        app.startFrame();
        rayTraceShader.use();
        
        camera.move(
            app.keyPressed(GLFW_KEY_W), 
            app.keyPressed(GLFW_KEY_S), 
            app.keyPressed(GLFW_KEY_D), 
            app.keyPressed(GLFW_KEY_A), 
            app.keyPressed(GLFW_KEY_SPACE),
            app.keyPressed(GLFW_KEY_LEFT_CONTROL),
            app.keyPressed(GLFW_KEY_LEFT_SHIFT)
        );
        camera.rotate(app.mouseX(), app.mouseY());
        
        if (camera.getIsMoving()){
            samples = 3;
            frameCount = 3;
        }
        
        glUniform1i(fragSamplesLoc, samples);
        glUniform1f(1, static_cast<float>(glfwGetTime()));
        glUniform2f(2, app.width(), app.height());
        glUniform1i(4, frameCount);
        glUniform3f(camPosLoc, camera.position().x, camera.position().y, camera.position().z);
        glUniform3f(camDirLoc, camera.lookDir().x, camera.lookDir().y, camera.lookDir().z);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, oldTexture);
        glUniform1i(3, 0);

        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        
        //Update oldTexture
        std::swap(texture, oldTexture);

        //Screen display (accumulation)
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        app.startFrame();
        accumulationShader.use();

        glBindTexture(GL_TEXTURE_2D, texture);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        if (app.keyPressed(GLFW_KEY_K))
            cout << "Frame Time: " << glfwGetTime() << "s with " << frameCount << " samples." << endl;

        app.eventAndSwapBuffers();
        samples = SAMPLES;
        frameCount += samples;
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteFramebuffers(1, &FBO);
    rayTraceShader.destroy();
    app.terminate();
    return EXIT_SUCCESS;
}