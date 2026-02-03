#include "camera.hpp"
#include <iostream>

using namespace std;

vec3 Camera::lookDir(){
    vec2 radAngles = vec2(radians(angles.x), radians(angles.y));
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
    isMoving = length(step) > 0;
    if (isMoving)
        pos += (sprinting ? 2 : 1) * moveSensitivity * normalize(moveForward + moveUp + moveRight);
    //cout << pos.x << " " << pos.z << endl;
}

void Camera::rotate(float mouseX, float mouseY){
    float deltaMouseX = mouseX - lastMouseX;
    float deltaMouseY = mouseY - lastMouseY;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
    angles += vec2(-deltaMouseX, -deltaMouseY) * lookSensitivity;
    angles.y = std::clamp(angles.y, -89.0f, 89.0f);
    isLooking = length(vec2(-deltaMouseX, -deltaMouseY)) > 0;
}

bool Camera::getIsMoving(){
    return isMoving || isLooking;
}