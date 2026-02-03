#version 430 core

struct Camera {
    vec3 pos;
    vec3 lookDir;
};

struct Mat {
    int type;
    vec3 color;
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
int bounces = 0;

#define NOHIT Hit(-1, vec3(0,0,0), Mat(0,vec3(0),0))
#define MAX_BOUNCES 10
//#define SAMPLES 1
#define EPS 1e-2

#define MAT_DIFF 0
#define MAT_METAL 1
#define MAT_EMIT 2

// UTILS ----------------

vec2 ratio(vec2 vec){
    return vec2(vec.x * winSize.x / winSize.y, vec.y);
}

vec3 initSeed() {
    return vec3(
        fract(sin(dot(vClipPos.xy , vec2(12.9898,78.233))) * 43758.5453 + time),
        fract(sin(dot(vClipPos.xy , vec2(93.9898,67.345))) * 43758.5453 + time),
        fract(sin(dot(vClipPos.xy , vec2(56.1234,12.345))) * 43758.5453 + time)
    );
}

float hash13(vec3 p3) {
    p3  = fract(p3 * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

float rand(inout vec3 seed) {
    float r = hash13(seed);
    seed += vec3(1.0, 1.0, 1.0);
    return r;
}

vec3 randomInSphere(inout vec3 seed) {
    float u = rand(seed);
    float v = rand(seed);
    float w = rand(seed);

    float theta = 2.0 * 3.14159265 * u;
    float phi   = acos(2.0 * v - 1.0);
    float r     = pow(w, 1.0/3.0);

    float x = r * sin(phi) * cos(theta);
    float y = r * sin(phi) * sin(theta);
    float z = r * cos(phi);

    return vec3(x, y, z);
}

vec3 randomInHemisphere(inout vec3 seed, vec3 normal) {
    vec3 rand = randomInSphere(seed);
    if (dot(rand, normal) < 0.0f)
        rand *= -1.0f;
    return rand;
}

vec3 randomOnLight(inout vec3 seed, vec3 direction, Light light){
    vec3 randomDir = randomInHemisphere(seed, direction);
    return light.pos + randomDir * light.rad;
}

vec3 randomCosineHemisphere(inout vec3 seed, vec3 normal, float randomizationFactor){
    vec3 val = normalize(normal) * (1 + EPS) + randomInSphere(seed) * randomizationFactor;
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
    return Hit(t, normalize(pos - light.pos), Mat(MAT_EMIT, light.color * light.intensity, 0));
}

Hit planeIntersect(Plane plane, Ray ray){
    vec3 rp = plane.origin - ray.origin;
    float t = dot(rp, plane.normal) / dot(ray.dir, plane.normal);
    vec3 relativePoint = ray.origin + t * ray.dir;
    if (t <= 0 || length(relativePoint - plane.origin) > 100){
        return NOHIT;
    }

    if (int((abs(relativePoint.x + 1) * 0.5) + int(abs(relativePoint.z + 1) * 0.5 + 1)) % 2 == 0) 
        plane.mat = Mat(plane.mat.type, vec3(0.83, 0.89, 0.95), 0.);

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

/*float p_direct(Ray ray, Light light){
    Hit obstacleHit = lightIntersect(light, ray);

    vec3 hitPoint = ray.origin + obstacleHit.t * ray.dir;
    if (obstacleHit.t > 0) {
        float A = 2 * light.rad * light.rad; // 2 pi R²
        float cos = dot(normalize(hitPoint - light.pos), -ray.dir);
        float r2 = dot(hitPoint - ray.origin, hitPoint - ray.origin);
        return r2 / (max(EPS,cos) * A);
    }
    else return 0;
}

float p_brdf(vec3 normal, vec3 rayDir){
    return max(EPS, dot(normal, rayDir));
}*/

void diffuse(World world, inout Hit hit, inout Ray ray, inout vec3 seed){
    Light light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
    ray.origin = ray.origin + hit.t * ray.dir + EPS * hit.normal;

    vec3 lightDir = normalize(light.pos - ray.origin);
    vec3 randomPoint = randomOnLight(seed, -lightDir, light);
    ray.dir = normalize(randomPoint - ray.origin);

    vec3 newDir = randomCosineHemisphere(seed, hit.normal, 1);
    ray.dir = newDir;
    ray.color *= 0.9 * vec4(hit.mat.color, 1);
}

void metal(Hit hit, in out Ray ray, in out vec3 seed){
    vec3 newDir = ray.dir - 2 * dot(ray.dir, hit.normal) * hit.normal;
    newDir = randomCosineHemisphere(seed, newDir, hit.mat.fuzz); 
    ray.origin = ray.origin + hit.t * ray.dir + hit.normal * EPS;
    ray.dir = newDir;
    ray.color *= 0.9 * vec4(hit.mat.color, 1);
}

void emit(inout Hit hit, inout Ray ray){
    ray.color *= vec4(hit.mat.color,1);
    hit.t = -2;
}

void computeLighting(World world, in out Hit hit, in out Ray ray, in out vec3 seed){
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
    //return vec4(0.7);
    return mix(vec4(bot, 1), vec4(top, 1), a);
}

vec4 rayColor(World world, in out vec3 seed, Ray ray){
    bounces = 0;
    Ray tracedRay = ray;
    while (bounces < MAX_BOUNCES){

        Hit hit = rayIntersection(world, tracedRay);
        computeLighting(world, hit, tracedRay, seed);

        if (ray.color == vec4(0,0,0,1)) break;

        if (length(ray.color.xyz) <= 0.01) break;

        if (hit.t < 0){
            if (hit.t > -2) tracedRay.color *= sky(tracedRay.dir);
            break; 
        }// Cas ou on tombe sur une light ou sur rien

        bounces += 1;
    }

    return tracedRay.color;
}

void main()
{
    vec3 seed = initSeed();
    vec2 pos = ratio(vClipPos.xy);
    vec2 uv = (vClipPos.xy + vec2(1)) * 0.5;

    Ray ray = Ray(
        camera.pos, 
        camera.lookDir, 
        vec4(1)
    );
    ray = fovRay(pos, ray);

    Mat sphereMat = Mat(MAT_DIFF, vec3(1, 0, 0), 0);
    Mat sphereMat2 = Mat(MAT_METAL, vec3(1.), 0);
    Mat planeMat = Mat(MAT_DIFF, vec3(0.7), 0);

    Sphere spheres[NUM_SPHERE];
    spheres[0] = Sphere(vec3(0,EPS,-5), 2, sphereMat);
    spheres[1] = Sphere(vec3(0,EPS,-9 - EPS), 2, planeMat);
    spheres[2] = Sphere(vec3(6, 2+EPS, -5), 4, sphereMat2);
    //spheres[3] = Sphere(vec3(0, EPS, -1 + EPS), 2, sphereMat);
    Plane planes[NUM_PLANE];
    planes[0] = Plane(vec3(0,-2,-3), vec3(0,1,0), planeMat);
    Light lights[NUM_LIGHT];
    lights[0] = Light(vec3(0,EPS,2), 1, vec3(1), 10);
    World world = World(spheres, planes, lights);

    vec4 currentColor;
    for(int i = 0; i < samples; i++){
        ray.origin = ray.origin + randomInSphere(seed) / winSize.y; // slight Anti-aliasing
        currentColor += rayColor(world, seed, ray);
    }

    FragColor = currentColor + float(frameCount - samples) * texture(screenTex, uv);
    FragColor /= frameCount;
}