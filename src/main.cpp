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

#define SAMPLES 2
#define BSDF_TYPES 3

ShaderProgram rayTraceShader;
ShaderProgram accumulationShader;
unsigned int VBO, VAO, EBO;
unsigned int FBO;
unsigned int texture;
unsigned int oldTexture;
int frameCount = 0;
int frameAccumulator = SAMPLES;
int samples = SAMPLES;
int bsdfType = 0;
int bsdfWeighting = 0;
unsigned int texWidth, texHeight;
const unsigned int RES_MUL = 3;
Camera camera(0.1, 0.3);
App app;

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

void genTexture(unsigned int width, unsigned int height){
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenTextures(1, &oldTexture);
    glBindTexture(GL_TEXTURE_2D, oldTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    samples = SAMPLES;
    frameAccumulator = SAMPLES;
    texWidth = width;
    texHeight = height;
}

void framebuffer_size_callback(GLFWwindow*, int width, int height)
{
    glViewport(0, 0, width, height);
    genTexture(width, height);

    rayTraceShader.use();
    int winSizeLoc = glGetUniformLocation(rayTraceShader.id(), "winSize");
    glUniform2f(winSizeLoc, width, height);
}

void init(){
    app.init(800, 600, "GLSL Project", framebuffer_size_callback);
    app.toggleCursor(false);
    
    rayTraceShader.create();
    rayTraceShader.load(GL_VERTEX_SHADER, "src/shaders/vertex.glsl");
    rayTraceShader.load(GL_FRAGMENT_SHADER, "src/shaders/frag.glsl");
    rayTraceShader.link();

    accumulationShader.create();
    accumulationShader.load(GL_VERTEX_SHADER, "src/shaders/screenVertex.glsl");
    accumulationShader.load(GL_FRAGMENT_SHADER, "src/shaders/screenFrag.glsl");
    accumulationShader.link();
    
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
    tie(VBO, VAO, EBO) = ShaderProgram::addData(vertices, indices);
    ShaderProgram::linkData(3, sizeof(float), 0);

    rayTraceShader.use();
    int winSizeLoc = glGetUniformLocation(rayTraceShader.id(), "winSize");
    glUniform2f(winSizeLoc, app.width(), app.height());
    
    glGenFramebuffers(1, &FBO);
    
    genTexture(app.width(), app.height());
    
    glDisable(GL_FRAMEBUFFER_SRGB);
}

void handleCamera(){
    if (!app.cursorIsHidden()){
        camera.hasStoppedMoving();
        return;
    }

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
}

void render(){
    //Current frame
    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);
    rayTraceShader.use();
    
    glUniform1f(1, static_cast<float>(glfwGetTime()));
    glUniform2f(2, texWidth, texHeight);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, oldTexture);
    glUniform1i(3, 0);
    glUniform1i(4, frameAccumulator);
    glUniform1i(5, bsdfType);
    glUniform1i(6, bsdfWeighting);
    glUniform1i(7, samples);
    int camPosLoc = glGetUniformLocation(rayTraceShader.id(), "camera.pos");
    glUniform3f(camPosLoc, camera.position().x, camera.position().y, camera.position().z);
    int camDirLoc = glGetUniformLocation(rayTraceShader.id(), "camera.lookDir");
    glUniform3f(camDirLoc, camera.lookDir().x, camera.lookDir().y, camera.lookDir().z);
    
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    
    //Update oldTexture
    std::swap(texture, oldTexture);
    
    //Screen display (accumulation)
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    accumulationShader.use();
    
    glBindTexture(GL_TEXTURE_2D, oldTexture);
    glUniform1i(1, 0);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
}

void inputs(){
    if (app.keyPressedOnce(GLFW_KEY_K, frameCount))
        cout << "Frame Time: " << glfwGetTime() << "s with " << frameAccumulator << " samples." << endl;
    
    if (app.keyPressedOnce(GLFW_KEY_ESCAPE, frameCount)){
        cout << "Toggling cursor." << endl;
        app.toggleCursor(true);
    }
    
    if (app.keyPressedOnce(GLFW_KEY_ENTER, frameCount))
        app.toggleCursor(false);

    if (app.keyPressedOnce(GLFW_KEY_RIGHT, frameCount)){
        frameAccumulator = samples;
        bsdfType = (bsdfType + 1) % BSDF_TYPES;
        cout << "BSDF type : " << bsdfType << endl;
    }
    
    if (app.keyPressedOnce(GLFW_KEY_LEFT, frameCount)){
        frameAccumulator = samples;
        bsdfType = (bsdfType - 1 + BSDF_TYPES) % BSDF_TYPES;
        cout << "BSDF type : " << bsdfType << endl;
    }
    
    if (app.keyPressedOnce(GLFW_KEY_UP, frameCount)){
        frameAccumulator = samples;
        bsdfWeighting = std::clamp(bsdfWeighting + 1, 0, 8);
        cout << "BSDF weighting : " << bsdfWeighting << endl;
    }
    
    if (app.keyPressedOnce(GLFW_KEY_DOWN, frameCount)){
        frameAccumulator = samples;
        bsdfWeighting = std::clamp(bsdfWeighting - 1, 0, 8);
        cout << "BSDF weighting : " << bsdfWeighting << endl;
    }
}

void end(){
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteFramebuffers(1, &FBO);
    rayTraceShader.destroy();
    app.terminate();
}

int main(){
    init();
    while(!app.shouldClose())
    {
        app.startFrame(frameCount);
        handleCamera();
        
        render();
        inputs();
        
        app.eventAndSwapBuffers();
        samples = SAMPLES;
        frameAccumulator += samples;
        frameCount++;

        if (camera.getIsMoving(frameCount)){
            samples = SAMPLES;
            frameAccumulator = SAMPLES;
            if (texWidth != app.width() / RES_MUL) 
                genTexture(app.width() / RES_MUL, app.height() / RES_MUL);
        }
        else if (texWidth != app.width()){
            samples = SAMPLES;
            frameAccumulator = SAMPLES;
            genTexture(app.width(), app.height());
        }
    }
    end();
    return EXIT_SUCCESS;
}