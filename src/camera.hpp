#include <glm/glm.hpp>

using namespace glm;

class Camera {
    private:
        float moveSensitivity = 0;
        float lookSensitivity = 0;
        vec3 pos = vec3(0);
        vec2 angles = vec2(0);
        float lastMouseX = 0;
        float lastMouseY = 0;
        bool isMoving = false;
        bool isLooking = false;

    public:
        Camera(float moveSensitivity, float lookSensitivity) 
            : moveSensitivity(moveSensitivity), lookSensitivity(lookSensitivity) {}
        void move(bool forward, bool backward, bool right, bool left, bool up, bool down, bool sprinting);
        void rotate(float mouseX, float mouseY);
        vec3 lookDir();
        vec3 position() {return pos;};
        bool getIsMoving();
};