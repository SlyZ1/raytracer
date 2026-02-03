#version 430 core

in vec4 vClipPos;
out vec4 FragColor;
uniform sampler2D screenTex;

void main()
{
    vec2 uv = (vClipPos.xy + vec2(1)) * 0.5;
    FragColor = texture(screenTex, uv);
}