#include <glm/glm.hpp>

using namespace glm;

class Camera {
    private:
        float m_moveSensitivity = 0;
        float m_lookSensitivity = 0;
        vec3 m_pos = vec3(0);
        vec2 m_angles = vec2(0);
        float m_lastMouseX = 0;
        float m_lastMouseY = 0;
        bool m_isMoving = false;
        bool m_isLooking = false;
        int m_lastMovingFrame = 0;

    public:
        Camera(float moveSensitivity, float lookSensitivity) 
            : m_moveSensitivity(moveSensitivity), m_lookSensitivity(lookSensitivity) {}
        void move(bool forward, bool backward, bool right, bool left, bool up, bool down, bool sprinting);
        void rotate(float mouseX, float mouseY);
        void resetMousePos(float mouseX, float mouseY);
        vec3 lookDir();
        vec3 position() {return m_pos;};
        bool getIsMoving(int frame);
        void hasStoppedMoving() {m_isMoving = false; m_isLooking = false;};
};