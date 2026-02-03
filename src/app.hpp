#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;

class App {
    private:
        GLFWwindow *window;

    public:
        void init(int width, int height, const char *name);
        void setClearColor(float r, float g, float b, float a);
        void startFrame();
        void eventAndSwapBuffers();
        void toggleCursor(bool hide);
        bool shouldClose();
        bool keyPressed(int key);
        float mouseX();
        float mouseY();
        void terminate();
        unsigned int width();
        unsigned int height();
};