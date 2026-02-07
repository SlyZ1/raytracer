#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;

class App {
    private:
        GLFWwindow *m_window;
        bool m_cursorHidden = false;

    public:
        void init(int width, int height, const char *name, GLFWframebuffersizefun framebuffer_size_callback);
        void setClearColor(float r, float g, float b, float a);
        void startFrame(int frame);
        void eventAndSwapBuffers();
        void toggleCursor(bool hide);
        bool cursorIsHidden();
        bool shouldClose();
        bool keyPressed(int key);
        bool keyPressedOnce(int key, int frame);
        float mouseX();
        float mouseY();
        void terminate();
        unsigned int width();
        unsigned int height();
};