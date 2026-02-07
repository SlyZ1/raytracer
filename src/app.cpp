#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "app.hpp"

using namespace std;

void App::init(int width, int height, const char *name, GLFWframebuffersizefun framebuffer_size_callback){
    if (!glfwInit())
    {
        cerr << "Failed to initialize GLFW" << endl;
        exit(1);
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    m_window = glfwCreateWindow(width, height, "ZMMR", NULL, NULL);
    glfwSetWindowTitle(m_window, name);
    if (m_window == NULL)
    {
        cerr << "Failed to open GLFW window" << endl;
        exit(1);
    }
    glfwMakeContextCurrent(m_window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        cerr << "Failed to initialize GLAD" << endl;
        exit(1);
    }

    glViewport(0, 0, width, height);
    glfwSetFramebufferSizeCallback(m_window, framebuffer_size_callback);
}

void App::setClearColor(float r, float g, float b, float a){
    glClearColor(r, g, b, a);
}

void App::startFrame(int frame){
    if(keyPressedOnce(GLFW_KEY_ESCAPE, frame) && !m_cursorHidden)
        glfwSetWindowShouldClose(m_window, true);
    glClear(GL_COLOR_BUFFER_BIT);
}

void App::eventAndSwapBuffers(){
    glfwSwapBuffers(m_window);
    glfwPollEvents();
}

void App::toggleCursor(bool show){
    m_cursorHidden = !show;
    glfwSetInputMode(m_window, GLFW_CURSOR, !show ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
}

bool App::cursorIsHidden(){
    return m_cursorHidden;
}

bool App::shouldClose(){
    return glfwWindowShouldClose(m_window);
}

bool App::keyPressed(int key){
    return glfwGetKey(m_window, key) == GLFW_PRESS;
}

bool App::keyPressedOnce(int key, int frame){
    static int wasPressed[GLFW_KEY_LAST + 1] = {INT_MAX};

    bool isPressed = glfwGetKey(m_window, key) == GLFW_PRESS;

    if (isPressed && wasPressed[key] >= frame) {
        wasPressed[key] = frame;
        return true;
    }

    if (!isPressed) {
        wasPressed[key] = INT_MAX;
    }

    return false;
}

float App::mouseX(){
    double mouseX;
    glfwGetCursorPos(m_window, &mouseX, nullptr);
    return static_cast<float>(mouseX);
}

float App::mouseY(){
    double mouseY;
    glfwGetCursorPos(m_window, nullptr, &mouseY);
    return static_cast<float>(mouseY);
}

void App::terminate(){
    glfwTerminate();
}

unsigned int App::width(){
    int width;
    glfwGetWindowSize(m_window, &width, nullptr);
    return width;
}

unsigned int App::height(){
    int height;
    glfwGetWindowSize(m_window, nullptr, &height);
    return height;
}