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

#define SAMPLES 5


ShaderProgram rayTraceShader;
ShaderProgram accumulationShader;
unsigned int VBO, VAO, EBO;
unsigned int FBO;
unsigned int texture;
unsigned int oldTexture;
unsigned int texWidth, texHeight;
glm::vec3 modelPos = glm::vec3(2,-1,-1.5);
float metalColor[3] = {1,1,1};
float metallic = 0;
float roughness = 0;
int resMultiplier = 3;
int maxBounces = 6;
int frameCount = 0;
int frameAccumulator = 0;
int samples = 1;
int bounces = maxBounces;
int bsdfType = 3;
bool showUI = false;
bool debugBVH = false;
Camera camera(0.1, 0.3);
App app;

bool isRendering = false;
int renderSamples = 2048;

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
    frameAccumulator = 0;
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
    mesh->loadFromModel("models/suzanne.obj");
    vector<Triangle> triangles = mesh->getTriangles();
    GLuint trissbo;
    glGenBuffers(1, &trissbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, trissbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, triangles.size() * sizeof(Triangle), triangles.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, trissbo); // binding 0
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glUniform1i(glGetUniformLocation(rayTraceShader.id(), "numTriangles"), triangles.size());
    
    vector<int> triIndices(triangles.size());
    for (int i = 0; i < triangles.size(); i++) {
        triIndices[i] = i;
    }
    BVHNode* bvh = Mesh::computeBVH(triangles, triIndices, 0, triangles.size());
    vector<linBVHNode> linNodes = Mesh::lineariseBVH(bvh, triangles);
    GLuint bvhssbo;
    glGenBuffers(1, &bvhssbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bvhssbo);
    glBufferData(GL_SHADER_STORAGE_BUFFER, linNodes.size() * sizeof(linBVHNode), linNodes.data(), GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, bvhssbo); // binding 1
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glUniform1i(glGetUniformLocation(rayTraceShader.id(), "numBVHNodes"), linNodes.size());
    
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

    int debugBVHPos = glGetUniformLocation(rayTraceShader.id(), "debugBVH");
    glUniform1i(debugBVHPos, debugBVH ? 1 : 0);
    int modelPosPos = glGetUniformLocation(rayTraceShader.id(), "modelPos");
    glUniform3f(modelPosPos, modelPos.x, modelPos.y, modelPos.z);
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
    if (app.UIInteract() || isRendering) return;

    if (app.keyPressedOnce(GLFW_KEY_K, frameCount))
        cout << "Frame Time: " << glfwGetTime() << "s with " << frameAccumulator << " samples." << endl;
    
    if (app.keyPressedOnce(GLFW_KEY_P, frameCount)){
        showUI = !showUI;
        app.toggleCursor(showUI);
        camera.resetMousePos(app.mouseX(), app.mouseY());
    }
}

void dynamicResolution(){
    if (camera.getIsMoving(frameCount) || app.UIDrag()){
        resetFrame();
        bounces = 4;
        if (texWidth != app.width() / resMultiplier) 
            genTexture(app.width() / resMultiplier, app.height() / resMultiplier);
    }
    else if (texWidth != app.width()){
        bounces = maxBounces;
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

void BeginTwoColumnLayout()
{
    ImGui::BeginTable("##layout", 2, ImGuiTableFlags_SizingStretchProp);
    ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthStretch);
    ImGui::TableSetupColumn("Input", ImGuiTableColumnFlags_WidthFixed, 120);
}

void EndTwoColumnLayout()
{
    ImGui::EndTable();
}

void Label(const char* label)
{
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::Text("%s", label);
    ImGui::TableSetColumnIndex(1);
    ImGui::SetNextItemWidth(-FLT_MIN);
}

void LeftWindowUI(){
    if (!showUI) return;

    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(300, app.getIo()->DisplaySize.y), ImGuiCond_Always);
    ImGui::Begin("Parameters", (bool*)NULL, ImGuiWindowFlags_MenuBar);

    ImGui::BeginDisabled(isRendering);

    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("Image"))
        {
            if (ImGui::MenuItem("Render", "R")) app.exportImage();
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Technical Settings", ImGuiTreeNodeFlags_DefaultOpen)){
        BeginTwoColumnLayout();
        
        Label("Max Bounces");
        if (ImGui::InputInt("##Max Bounces", &maxBounces)) resetFrame();
        
        Label("Resolution Divider");
        ImGui::InputInt("##Resolution Divider", &resMultiplier);
        resMultiplier = glm::max(resMultiplier, 1);
        
        EndTwoColumnLayout();
    }

    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Material Settings", ImGuiTreeNodeFlags_DefaultOpen)){
        BeginTwoColumnLayout();

        const char* items[] = { 
            "Lambert", 
            "(QON) Qualitative Oren-Nayar", 
            "(FON) Fujii-Oren-Nayar", 
            "(EON) Energy-Preserving Oren-Nayar"
        };
        Label("Diffuse Model");
        if (ImGui::Combo("##Diffuse Model", &bsdfType, items, IM_ARRAYSIZE(items)))
            resetFrame();

        Label("Ball's Roughness");
        if (ImGui::SliderFloat("##Ball's Roughness", &roughness, 0, 1))
            resetFrame();
        
        Label("Metal Color");
        if (ImGui::ColorEdit3("##Metal Color", metalColor))
            resetFrame();
    
        Label("Metallic");
        if (ImGui::SliderFloat("##Metallic", &metallic, 0, 1))
            resetFrame();

        EndTwoColumnLayout();
    }

    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Model Settings", ImGuiTreeNodeFlags_DefaultOpen)){
        BeginTwoColumnLayout();

        Label("Debug BVH");
        if (ImGui::Checkbox("##Debug BVH", &debugBVH))
            resetFrame();

        Label("Model Position");
        if (ImGui::DragFloat3("##Model Position", &modelPos[0]))
            resetFrame();

        EndTwoColumnLayout();
    }

    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Render Settings", ImGuiTreeNodeFlags_DefaultOpen)){
        
        if (ImGui::Button("Render", ImVec2(-FLT_MIN, 0))){
            isRendering = true;
            samples = SAMPLES;
            resetFrame();
        }
        
        BeginTwoColumnLayout();
        
        Label("Render Samples");
        ImGui::InputInt("##Render Samples", &renderSamples);
        
        EndTwoColumnLayout();
        
        if (isRendering) {
            ImGui::ProgressBar((float)frameAccumulator / renderSamples, ImVec2(-FLT_MIN, 0));

            string text = to_string(frameAccumulator) + " / " + to_string(renderSamples);
            float windowWidth = ImGui::GetContentRegionAvail().x;
            float textWidth   = ImGui::CalcTextSize(text.c_str()).x;
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (windowWidth - textWidth) * 0.5f);
            ImGui::Text(text.c_str());

            if (frameAccumulator >= renderSamples){
                app.exportImage();
                samples = 1;
                isRendering = false;
            }
        }
    }
    
    ImGui::EndDisabled();
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
        
        LeftWindowUI();
        
        frameAccumulator += samples;
        frameCount++;
        dynamicResolution();
        app.endFrame();
    }
    end();
    return EXIT_SUCCESS;
}