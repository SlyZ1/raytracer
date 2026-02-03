#version 430 core
layout (location = 0) in vec3 aPos;
layout (location = 1) uniform float time;
out vec4 vClipPos;

void main()
{
    gl_Position = vec4(aPos, 1.0);
    vClipPos = gl_Position;
}