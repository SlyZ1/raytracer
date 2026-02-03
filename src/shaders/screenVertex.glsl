#version 430 core
layout (location = 0) in vec3 aPos;
out vec4 vClipPos;

void main()
{
    gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);
    vClipPos = gl_Position;
}