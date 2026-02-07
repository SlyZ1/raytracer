#version 430 core

in vec4 vClipPos;
out vec4 FragColor;
layout (location = 1) uniform sampler2D screenTex;

void main()
{
    vec2 uv = (vClipPos.xy + vec2(1)) * 0.5;

    vec3 color = vec3(0);
    int radius = 100;

    for (int i = -radius; i <= radius; ++i) {
        if (i == 0) continue;
        vec2 offset = i * vec2(1,0) * (1/600);
        color += texture(screenTex, uv + offset).rgb / (2 * radius);
    }
    
    FragColor = vec4(color, 1);
}