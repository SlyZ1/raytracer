#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;

class App {
    private:
        GLFWwindow *m_window;
        bool m_cursorHidden = false;
        ImGuiIO* m_io = {};

    public:
        void init(int width, int height, const char *name, GLFWframebuffersizefun framebuffer_size_callback);
        void setClearColor(float r, float g, float b, float a);
        void startFrame(int frame);
        void endFrame();
        void toggleCursor(bool hide);
        void exportImage(string additionalPath = "", string name = "image.png");
        bool cursorIsHidden();
        bool shouldClose();
        bool keyPressed(int key);
        bool keyPressedOnce(int key, int frame);
        bool UIInteract();
        bool UIDrag();
        ImGuiIO* getIo();
        float mouseX();
        float mouseY();
        void terminate();
        unsigned int width();
        unsigned int height();
};