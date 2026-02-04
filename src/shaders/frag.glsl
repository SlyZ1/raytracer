#version 430 core

struct Camera {
    vec3 pos;
    vec3 lookDir;
};

struct Mat {
    int type;
    vec3 color;
    float intensity;
    float fuzz;
};

struct Sphere {
    vec3 pos;
    float rad;
    Mat mat;
};

struct Light {
    vec3 pos;
    float rad;
    vec3 color;
    float intensity;
};

struct Plane {
    vec3 origin;
    vec3 normal;
    Mat mat;
};

struct Ray {
    vec3 oldNormal;
    vec3 origin;
    vec3 dir;
    vec4 color;
};

struct Hit {
    float t;
    vec3 normal;
    Mat mat;
};

#define NUM_SPHERE 3
#define NUM_PLANE 1
#define NUM_LIGHT 1

struct World {
    Sphere spheres[NUM_SPHERE];
    Plane planes[NUM_PLANE];
    Light lights[NUM_LIGHT];
};

out vec4 FragColor;
layout (location = 1) uniform float time;
layout (location = 2) uniform vec2 winSize;
layout (location = 3) uniform sampler2D screenTex;
layout (location = 4) uniform int frameCount;
uniform int samples;
uniform Camera camera;
in vec4 vClipPos;

#define NOHIT Hit(-1, vec3(0,0,0), Mat(0,vec3(0),0,0))
#define MAX_BOUNCES 1
//#define SAMPLES 1
#define EPS 1e-4
#define PI 3.14159265

#define MAT_DIFF 0
#define MAT_METAL 1
#define MAT_EMIT 2

// UTILS ----------------

vec2 ratio(vec2 vec){
    return vec2(vec.x * winSize.x / winSize.y, vec.y);
}

uint pcg_hash(uint v) {
    v = v * 747796405u + 2891336453u;
    uint word = ((v >> ((v >> 28u) + 4u)) ^ v) * 277803737u;
    return (word >> 22u) ^ word;
}

uint initSeed(uvec2 pos, uint frame) {
    uint v = pos.x + pos.y * 4096u + frame * 1315423911u;
    return pcg_hash(v);
}

float rand(inout uint seed) {
    seed = pcg_hash(seed);
    return float(seed) * (1.0 / 4294967296.0);
}

vec3 randomInSphere(inout uint seed) {
    float u = rand(seed);
    float v = rand(seed);
    float w = rand(seed);

    float theta = 2.0 * PI * u;
    float phi   = acos(2.0 * v - 1.0);
    float r     = pow(w, 1.0/3.0);

    float x = r * sin(phi) * cos(theta);
    float y = r * sin(phi) * sin(theta);
    float z = r * cos(phi);

    return vec3(x, y, z);
}

vec3 randomOnUnitSphere(inout uint seed) {
    float u = rand(seed);
    float v = rand(seed);

    float theta = 2.0 * PI * u;
    float phi   = acos(2.0 * v - 1.0);

    float x = sin(phi) * cos(theta);
    float y = sin(phi) * sin(theta);
    float z = cos(phi);

    return vec3(x, y, z);
}

vec3 randomOnUnitHemiphere(inout uint seed, vec3 normal){
    vec3 dir = randomOnUnitSphere(seed);
    if (dot(dir, normal) < 0)
        dir *= -1;
    return dir;
}

vec3 randomCosineHemisphere(inout uint seed, vec3 normal, float randomizationFactor){
    vec3 val = normalize(normal) * (1 + EPS) + randomOnUnitSphere(seed) * randomizationFactor;
    return normalize(val);
}

// INTERSECTIONS

Hit sphereIntersect(Sphere sphere, Ray ray){
    vec3 oc = ray.origin - sphere.pos;
    float b = dot(oc, ray.dir);
    float c = dot(oc, oc) - sphere.rad * sphere.rad;
    float h = b*b - c;
    if (h < 0.) return NOHIT;
    float t = -b - sqrt(h);
    if (t <= 0) return NOHIT;
    vec3 pos = ray.origin + t * ray.dir;
    return Hit(t, normalize(pos - sphere.pos), sphere.mat);
}

Hit lightIntersect(Light light, Ray ray){
    vec3 oc = ray.origin - light.pos;
    float b = dot(oc, ray.dir);
    float c = dot(oc, oc) - light.rad * light.rad;
    float h = b*b - c;
    if (h < 0.) return NOHIT;
    float t = -b - sqrt(h);
    if (t <= 0) return NOHIT;
    vec3 pos = ray.origin + t * ray.dir;
    return Hit(t, normalize(pos - light.pos), Mat(MAT_EMIT, light.color, light.intensity, 0));
}

Hit planeIntersect(Plane plane, Ray ray){
    vec3 rp = plane.origin - ray.origin;
    float t = dot(rp, plane.normal) / dot(ray.dir, plane.normal);
    vec3 relativePoint = ray.origin + t * ray.dir;
    vec3 difference = relativePoint - plane.origin;
    if (t <= 0 || dot(difference, difference) > 10000){
        return NOHIT;
    }

    if (int((abs(relativePoint.x + 1) * 0.5) + int(abs(relativePoint.z + 1) * 0.5 + 1)) % 2 == 0) 
        plane.mat = Mat(plane.mat.type, vec3(0.83, 0.89, 0.95), 1, 0.);

    return Hit(t, plane.normal, plane.mat);
}

Hit rayIntersection(World world, Ray ray){
    Hit hit = NOHIT;
    hit.t = 100000;
    for(int i = 0; i < NUM_SPHERE; i += 1){
        Hit sphereHit = sphereIntersect(world.spheres[i], ray);
        if (sphereHit.t > 0 && sphereHit.t < hit.t) hit = sphereHit;
    }
    for(int i = 0; i < NUM_SPHERE; i += 1){
        Hit planeHit = planeIntersect(world.planes[i], ray);
        if (planeHit.t > 0 && planeHit.t < hit.t) hit = planeHit;
    }
    for(int i = 0; i < NUM_LIGHT; i += 1){
        Hit lightHit = lightIntersect(world.lights[i], ray);
        if (lightHit.t > 0 && lightHit.t < hit.t) hit = lightHit;
    }

    if(hit.t == 100000) hit.t = -1;

    return hit;
}

// LIGHTS

vec3 f_r(Hit hit){
    return hit.mat.color / PI;
}

float p_direct(Light light, float distance, float cosLight, World world, Ray ray){
    Hit hit = rayIntersection(world, ray);
    if (hit.t < 0 || hit.t - distance > 0.01) return 0;
    return distance * distance / (cosLight * 2 * PI * light.rad * light.rad);
}

float p_bsdf(vec3 normal, vec3 w){
    return max(EPS, dot(normal, w)) / PI;
}

float power_heuristic(float p1, float p2, float beta){
    return pow(p1, beta) / (pow(p1, beta) + pow(p2, beta));
}

void diffuse(World world, inout Hit hit, inout Ray ray, inout uint seed){
    Light light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
    ray.origin = ray.origin + hit.t * ray.dir + EPS * hit.normal;
    float r = rand(seed);
    float wdirect = 1;
    float wbsdf = 0.0;
    if (r < wdirect){
        vec3 nLight = randomOnUnitHemiphere(seed, ray.origin - light.pos);
        vec3 lightPoint = light.pos + nLight * light.rad;
        vec3 newDir = normalize(lightPoint - ray.origin);
        ray.dir = newDir;

        float cosR = max(dot(hit.normal, newDir), 0.0);
        float cosL = abs(dot(-newDir, nLight));
        vec3 Le = light.color * light.intensity;
        
        float pdirect = p_direct(light, length(lightPoint - ray.origin), cosL, world, ray);
        if (pdirect > 0){
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = power_heuristic(wdirect * pdirect, wbsdf * pbsdf, 2);

            vec3 contrib = Le * f_r(hit) * cosR / (pdirect * wdirect);
            ray.color += vec4(contrib * weight, 1);
            hit.t = -2;
            return;
        }
    }
    else{
        vec3 newDir = randomCosineHemisphere(seed, hit.normal, 1);
        ray.dir = newDir;
        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * nextHit.mat.intensity;
            float cosR = max(dot(hit.normal, newDir), 0.0);
            float cosL = abs(dot(-newDir, nextHit.normal));

            float pdirect = p_direct(light, nextHit.t, cosL, world, ray);
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = power_heuristic(wbsdf * pbsdf, wdirect * pdirect, 2);

            vec3 contrib = Le * f_r(hit) * cosR / (pbsdf * wbsdf);
            ray.color += vec4(contrib * weight, 1);
            hit.t = -2;
            return;
        }
    }

    vec3 newDir = randomCosineHemisphere(seed, hit.normal, 1);
    ray.dir = newDir;
    ray.color *= vec4(f_r(hit), 1);
    ray.color = vec4(1,1,0,1);
}

void metal(Hit hit, in out Ray ray, in out uint seed){
    vec3 newDir = ray.dir - 2 * dot(ray.dir, hit.normal) * hit.normal;
    newDir = randomCosineHemisphere(seed, newDir, hit.mat.fuzz); 
    ray.origin = ray.origin + hit.t * ray.dir + hit.normal * EPS;
    ray.dir = newDir;
    ray.color *= vec4(hit.mat.color, 1);
}

void emit(inout Hit hit, inout Ray ray){
    ray.color += vec4(hit.mat.color * hit.mat.intensity,1);
    hit.t = -2;
}

void computeLighting(World world, in out Hit hit, in out Ray ray, in out uint seed){
    if (hit.t < 0) return;

    switch (hit.mat.type){
        case MAT_DIFF:
            diffuse(world, hit, ray, seed);
            break;
        case MAT_METAL:
            metal(hit, ray, seed);
            break;
        case MAT_EMIT:
            emit(hit, ray);
            break;
    }
}

// RAY TRACING --------------------

Ray fovRay(vec2 pos, Ray ray){
    float fov = 90;
    float angle = radians(fov / 2);

    vec3 right = normalize(cross(camera.lookDir, vec3(0,1,0)));
    vec3 up = normalize(cross(right, camera.lookDir));
    vec3 deviation = right * pos.x + up * pos.y;
    deviation = deviation * tan(angle);
    ray.dir = camera.lookDir + deviation;
    ray.dir = normalize(ray.dir);

    return ray;
}

vec4 sky(vec3 lookingAt){
    float a = (dot(normalize(lookingAt), vec3(0,1,0)) + 1) * 0.5;
    vec3 top = vec3(0.17, 0.24, 0.31) * 1.1;
    vec3 bot = vec3(1.0, 0.49, 0.37);
    return mix(vec4(bot, 1), vec4(top, 1), a);
}

vec4 rayColor(World world, in out uint seed, Ray ray){
    Ray tracedRay = ray;
    for (int i = 0; i < MAX_BOUNCES; i++){

        Hit hit = rayIntersection(world, tracedRay);
        computeLighting(world, hit, tracedRay, seed);

        if (hit.t < 0){
            //tracedRay.color *= sky(tracedRay.dir);
            return tracedRay.color;
        }// Cas ou on tombe sur une light ou sur rien
    }

    return vec4(0);
}

void main()
{
    uint seed = initSeed(uvec2(gl_FragCoord.xy), frameCount);
    vec2 pos = ratio(vClipPos.xy);
    vec2 uv = (vClipPos.xy + vec2(1)) * 0.5;

    Ray ray = Ray(
        vec3(0),
        camera.pos, 
        camera.lookDir, 
        vec4(0.0)
    );
    ray = fovRay(pos, ray);

    Mat sphereMat = Mat(MAT_DIFF, vec3(1, 0, 0), 1, 0);
    Mat sphereMat2 = Mat(MAT_METAL, vec3(1.), 1, 0);
    Mat planeMat = Mat(MAT_DIFF, vec3(0.7), 1, 0);

    Sphere spheres[NUM_SPHERE];
    spheres[0] = Sphere(vec3(0,EPS,-5), 2, sphereMat);
    spheres[1] = Sphere(vec3(0,EPS,-10), 2, planeMat);
    spheres[2] = Sphere(vec3(6, 2+EPS, -5), 4, sphereMat2);
    
    Plane planes[NUM_PLANE];
    planes[0] = Plane(vec3(0,-2,-3), vec3(0,1,0), planeMat);
    Light lights[NUM_LIGHT];
    lights[0] = Light(vec3(0,EPS,2), 1, vec3(1), 10);
    World world = World(spheres, planes, lights);

    vec4 currentColor;
    for(int i = 0; i < samples; i++){
        ray.origin = ray.origin + 1.5*randomInSphere(seed) / winSize.y; // slight Anti-aliasing
        currentColor += rayColor(world, seed, ray);
    }

    FragColor = currentColor + float(frameCount - samples) * texture(screenTex, uv);
    FragColor /= frameCount;
}