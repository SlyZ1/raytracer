#include "camera.hpp"
#include <iostream>

using namespace std;

vec3 Camera::lookDir(){
    vec2 radAngles = vec2(radians(m_angles.x), radians(m_angles.y));
    vec3 lookDir = vec3(- sin(radAngles.x) * cos(radAngles.y), sin(radAngles.y), - cos(radAngles.x) * cos(radAngles.y));
    return normalize(lookDir);
}

void Camera::move(bool forward, bool backward, bool right, bool left, bool up, bool down, bool sprinting){
    vec3 step = vec3(
        (int)right - (int)left,
        (int)up - (int)down,
        (int)forward - (int)backward
    );
    vec3 moveForward = step.z * lookDir();
    vec3 moveUp = vec3(0, step.y, 0);
    vec3 moveRight = step.x * normalize(cross(lookDir(), vec3(0,1,0)));
    m_isMoving = length(step) > 0;
    if (m_isMoving)
        m_pos += (sprinting ? 2 : 1) * m_moveSensitivity * normalize(moveForward + moveUp + moveRight);
}

void Camera::resetMousePos(float mouseX, float mouseY){
    m_lastMouseX = mouseX;
    m_lastMouseY = mouseY;
}

void Camera::rotate(float mouseX, float mouseY){
    float deltaMouseX = mouseX - m_lastMouseX;
    float deltaMouseY = mouseY - m_lastMouseY;
    m_lastMouseX = mouseX;
    m_lastMouseY = mouseY;
    m_angles += vec2(-deltaMouseX, -deltaMouseY) * m_lookSensitivity;
    m_angles.y = std::clamp(m_angles.y, -89.0f, 89.0f);
    m_isLooking = length(vec2(-deltaMouseX, -deltaMouseY)) > 0;
}

bool Camera::getIsMoving(int frame){
    bool result = m_isMoving || m_isLooking;
    if (result) m_lastMovingFrame = frame;
    return frame - m_lastMovingFrame < 10;
}