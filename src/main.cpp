#include <iostream>
#include <fstream>
#include <sstream>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_opengl3.h>
#include <imgui/imgui_impl_glfw.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "app.hpp"
#include "shader_program.hpp"
#include "camera.hpp"
#include "mesh.hpp"

using namespace std;

#define SAMPLES 2


ShaderProgram rayTraceShader;
ShaderProgram accumulationShader;
unsigned int VBO, VAO, EBO;
unsigned int FBO;
unsigned int texture;
unsigned int oldTexture;
unsigned int texWidth, texHeight;
float metalColor[3] = {1,1,1};
float metallic = 0;
float roughness = 0;
int resMultiplier = 2;
int maxBounces = 15;
int frameCount = 0;
int frameAccumulator = SAMPLES;
int samples = SAMPLES;
int bounces = maxBounces;
int bsdfType = 3;
bool showUI = false;
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

void resetFrame(){
    frameAccumulator = samples;
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
    resetFrame();
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
    app.init(1080, 720, "Basic Raytracer", framebuffer_size_callback);
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

    Mesh* mesh = new Mesh();
    mesh->loadFromModel("models/Cube.obj");
    vector<Triangle> triangles = mesh->getTriangles();
    GLuint ssbo;
    glGenBuffers(1, &ssbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, triangles.size() * sizeof(Triangle), triangles.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo); // binding 0
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glUniform1i(10, triangles.size());
    
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
    
    glUniform4f(1, metalColor[0], metalColor[1], metalColor[2], metallic);
    glUniform2f(2, texWidth, texHeight);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, oldTexture);
    glUniform1i(3, 0);
    glUniform1i(4, frameAccumulator);
    glUniform1i(5, bsdfType);
    glUniform1f(6, roughness);
    glUniform1i(7, samples);
    glUniform1i(9, bounces);
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
    if (app.UIInteract()) return;

    if (app.keyPressedOnce(GLFW_KEY_K, frameCount))
        cout << "Frame Time: " << glfwGetTime() << "s with " << frameAccumulator << " samples." << endl;
    
    if (app.keyPressedOnce(GLFW_KEY_P, frameCount)){
        showUI = !showUI;
        app.toggleCursor(showUI);
        camera.resetMousePos(app.mouseX(), app.mouseY());
    }
}

void dynamicResolution(){
    if (app.UIInteract()) samples = 1;
    else samples = SAMPLES;

    if (camera.getIsMoving(frameCount) || app.UIDrag()){
        samples = 1;
        resetFrame();
        bounces = 4;
        if (texWidth != app.width() / resMultiplier) 
            genTexture(app.width() / resMultiplier, app.height() / resMultiplier);
    }
    else if (texWidth != app.width()){
        bounces = maxBounces;
        samples = SAMPLES;
        resetFrame();
        genTexture(app.width(), app.height());
    }
}

void end(){
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteFramebuffers(1, &FBO);
    rayTraceShader.destroy();
    app.terminate();
}

void UI(){
    if (!showUI) return;

    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(300, app.getIo()->DisplaySize.y), ImGuiCond_Always);

    ImGui::Begin("Parameters", (bool*)NULL, ImGuiWindowFlags_MenuBar);

    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("Image"))
        {
            if (ImGui::MenuItem("Render", "R")) app.exportImage();
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    if (ImGui::InputInt("Max Bounces", &maxBounces)) resetFrame();
    ImGui::InputInt("Resolution Divider", &resMultiplier);
    resMultiplier = glm::max(resMultiplier, 1);

    const char* items[] = { 
        "Lambert", 
        "Qualitative Oren-Nayar (QON)", 
        "Fujii-Oren-Nayar (FON)", 
        "Energy-Preserving Oren-Nayar (EON)"
    };
    if (ImGui::Combo("Diffuse Model", &bsdfType, items, IM_ARRAYSIZE(items)))
        resetFrame();

    if (ImGui::SliderFloat("Ball's Roughness", &roughness, 0, 1))
        resetFrame();
    
    if (ImGui::ColorEdit3("Metal Color", metalColor)) 
        resetFrame();

    if (ImGui::SliderFloat("Metallic", &metallic, 0, 1)) 
        resetFrame();

    ImGui::End();
}

int main(){
    init();
    while(!app.shouldClose())
    {
        app.startFrame(frameCount);
        handleCamera();
        
        render();
        inputs();
        
        UI();
        
        samples = SAMPLES;
        frameAccumulator += samples;
        frameCount++;
        dynamicResolution();
        app.endFrame();
    }
    end();
    return EXIT_SUCCESS;
}