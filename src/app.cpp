#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "app.hpp"

using namespace std;

void framebuffer_size_callback(GLFWwindow*, int width, int height)
{
    glViewport(0, 0, width, height);
}

void App::init(int width, int height, const char *name){
    if (!glfwInit())
    {
        cerr << "Failed to initialize GLFW" << endl;
        exit(1);
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    window = glfwCreateWindow(width, height, "ZMMR", NULL, NULL);
    glfwSetWindowTitle(window, name);
    if (window == NULL)
    {
        cerr << "Failed to open GLFW window" << endl;
        exit(1);
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        cerr << "Failed to initialize GLAD" << endl;
        exit(1);
    }

    glViewport(0, 0, width, height);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
}

void App::setClearColor(float r, float g, float b, float a){
    glClearColor(r, g, b, a);
}

void App::startFrame(){
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    glClear(GL_COLOR_BUFFER_BIT);
}

void App::eventAndSwapBuffers(){
    glfwSwapBuffers(window);
    glfwPollEvents();
}

void App::toggleCursor(bool show){
    glfwSetInputMode(window, GLFW_CURSOR, !show ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
}

bool App::shouldClose(){
    return glfwWindowShouldClose(window);
}

bool App::keyPressed(int key){
    return glfwGetKey(window, key) == GLFW_PRESS;
}

float App::mouseX(){
    double mouseX;
    glfwGetCursorPos(window, &mouseX, nullptr);
    return static_cast<float>(mouseX);
}

float App::mouseY(){
    double mouseY;
    glfwGetCursorPos(window, nullptr, &mouseY);
    return static_cast<float>(mouseY);
}

void App::terminate(){
    glfwTerminate();
}

unsigned int App::width(){
    int width;
    glfwGetWindowSize(window, &width, nullptr);
    return width;
}

unsigned int App::height(){
    int height;
    glfwGetWindowSize(window, nullptr, &height);
    return height;
}