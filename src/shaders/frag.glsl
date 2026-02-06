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
    vec3 origin;
    vec3 dir;
    vec3 throughput;
    vec3 radiance;
};

struct Hit {
    float t;
    vec3 normal;
    Mat mat;
    bool frontface;
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
layout (location = 5) uniform int bsdfType;
layout (location = 6) uniform int bsdfWeighting;
layout (location = 7) uniform int samples;
uniform Camera camera;
in vec4 vClipPos;

#define NOHIT Hit(-1, vec3(0,0,0), Mat(0,vec3(0),0,0), true)
#define MAX_BOUNCES 7
//#define SAMPLES 1
#define EPS 1e-4
#define PI 3.14159265

#define MAT_DIFF 0
#define MAT_METAL 1
#define MAT_EMIT 2

#define BSDF_LAMBERT 0
#define BSDF_LAMBERT_WRAP 1
#define BSDF_OREN_NAYAR 2

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
    bool frontface = true;
    if (h < 0.) return NOHIT;
    float t = -b - sqrt(h);
    if (t <= 0){
        /*t = -b + sqrt(h);
        frontface = false;*/
        if (t <= 0) return NOHIT;
    }
    vec3 pos = ray.origin + t * ray.dir;
    vec3 normal = normalize(pos - sphere.pos);
    if (!frontface) normal *= -1;
    return Hit(t, normal, sphere.mat, frontface);
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
    return Hit(t, normalize(pos - light.pos), Mat(MAT_EMIT, light.color, light.intensity, 0), true);
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

    return Hit(t, plane.normal, plane.mat, true);
}

Hit rayIntersection(World world, Ray ray){
    Hit hit = NOHIT;
    hit.t = 100000;
    for(int i = 0; i < NUM_SPHERE; i += 1){
        Hit sphereHit = sphereIntersect(world.spheres[i], ray);
        if (sphereHit.t > 0 && sphereHit.t < hit.t) hit = sphereHit;
    }
    for(int i = 0; i < NUM_PLANE; i += 1){
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

vec3 sky(vec3 lookingAt){
    float a = (dot(normalize(lookingAt), vec3(0,1,0)) + 1) * 0.5;
    vec3 top = vec3(0.17, 0.24, 0.31) * 1.1;
    vec3 bot = vec3(1.0, 0.49, 0.37);
    return mix(bot, top, a);
}

vec3 lambert(Hit hit){
    return hit.mat.color / PI;
}

float wrap(float ndotl, float wrapFactor){
    float wrap = max(0, (ndotl + wrapFactor) / (1 + wrapFactor));
    return wrap;
}

vec3 oren_nayar(Hit hit, vec3 normal, vec3 lightDir, vec3 viewDir, float roughness){
    float sigma2 = roughness * roughness;
    float A = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
    float B = 0.45 * sigma2 / (sigma2 + 0.09);

    float LdotN = max(dot(normal, lightDir), 0.0);
    float VdotN = max(dot(normal, viewDir), 0.0);
    float alpha = max(acos(LdotN), acos(VdotN));
    float beta = min(acos(LdotN), acos(VdotN));
    float gamma = dot(viewDir, lightDir);

    float angleTerm = A + B * max(0.0, gamma) * sin(alpha) * tan(beta);

    return hit.mat.color * angleTerm * LdotN / PI;
}

float p_direct(Light light, float distance, float cosLight){
    return distance * distance / (cosLight * 2 * PI * light.rad * light.rad);
}

float shadow_hit(Light light, World world, Ray ray){
    Hit hit = rayIntersection(world, ray);
    vec3 hitPos = ray.origin + ray.dir * hit.t;
    if (hit.t > 0 && length(hitPos - light.pos) > light.rad + sqrt(EPS)) return 0;
    return 1;
}

float p_bsdf(vec3 normal, vec3 w){
    return max(EPS, dot(normal, w)) / PI;
}

float p_metal(vec3 reflect, vec3 w, float fuzz){
    return max(EPS, dot(reflect, w)) / (PI * pow(fuzz, 2));
}

void diffuse(World world, inout Hit hit, inout Ray ray, inout uint seed){
    ray.origin = ray.origin + hit.t * ray.dir + EPS * hit.normal;

    Light light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
    vec3 nLight = randomOnUnitHemiphere(seed, ray.origin - light.pos);
    vec3 lightPoint = light.pos + nLight * light.rad;
    vec3 lightDir = normalize(lightPoint - ray.origin);

    float wdirect = bsdfWeighting / 8.0;
    float wbsdf = 1 - bsdfWeighting / 8.0;
    float r = rand(seed);

    if (r < wbsdf){
        // BSDF sampling
        vec3 newDir = randomCosineHemisphere(seed, hit.normal, 1);
        vec3 f_r = vec3(0);
        if (bsdfType == BSDF_OREN_NAYAR) f_r = oren_nayar(hit, hit.normal, lightDir, ray.dir, 1);
        if (bsdfType == BSDF_LAMBERT || bsdfType == BSDF_LAMBERT_WRAP) f_r = lambert(hit);
        ray.throughput *= f_r * PI / wbsdf;
        ray.dir = newDir;

        if (bsdfType == BSDF_LAMBERT_WRAP){
            float wrap = wrap(dot(hit.normal, lightDir), 0.5);
            float cos = max(dot(hit.normal, newDir), wrap);
            ray.throughput *= cos / dot(newDir, hit.normal);
        }
        // Check if we hit a light with the BSDF sampling
        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * nextHit.mat.intensity;
            float cosL = max(dot(-newDir, nextHit.normal), 1);

            float pdirect = p_direct(light, nextHit.t, cosL)
                            * shadow_hit(light, world, ray);
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = wbsdf * pbsdf / (wbsdf * pbsdf + wdirect * pdirect);

            ray.radiance += ray.throughput * Le * weight;
            hit.t = -2;
        }
    } 
    else {
        // Direct lighting
        vec3 newDir = lightDir;
        float cosR = max(dot(hit.normal, newDir), 0.0);
        float cosL = max(dot(-newDir, nLight), 0);
        vec3 f_r = vec3(0);
        if (bsdfType == BSDF_OREN_NAYAR) f_r = oren_nayar(hit, hit.normal, lightDir, ray.dir, 1);
        if (bsdfType == BSDF_LAMBERT || bsdfType == BSDF_LAMBERT_WRAP) f_r = lambert(hit);
        ray.throughput *= f_r * cosR;
        ray.dir = newDir;

        vec3 Le = light.color * light.intensity;
        
        float distance = length(lightPoint - ray.origin);
        float pdirect = p_direct(light, distance, cosL);

        if (shadow_hit(light, world, ray) > 0){
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = 1.0 / (wdirect * pdirect + wbsdf * pbsdf);

            if (bsdfType == BSDF_LAMBERT_WRAP){
                ray.throughput *= wrap(dot(hit.normal, lightDir), 0.5) / cosR;
            }
            ray.radiance += ray.throughput * Le * weight;
        }
        else{
            ray.throughput /= pdirect * wdirect;
        }
        hit.t = -2;
    }
}

void metal(World world, Hit hit, in out Ray ray, in out uint seed){
    hit.t = -2;
    return;
    Light light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
    ray.origin = ray.origin + hit.t * ray.dir + hit.normal * EPS;
    vec3 reflection = ray.dir - 2 * dot(ray.dir, hit.normal) * hit.normal;
    float wdirect = 0.5 * hit.mat.fuzz;
    float wmetal = 0.5 + 0.5 * (1 - hit.mat.fuzz);
    wdirect = 0;
    wmetal = 1;
    float r = rand(seed);

    if (r < wmetal){
        vec3 newDir = randomCosineHemisphere(seed, reflection, hit.mat.fuzz);
        ray.dir = newDir;
        ray.throughput *= hit.mat.color / wmetal;

        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * nextHit.mat.intensity;
            float cosL = max(dot(-newDir, nextHit.normal), 0);

            float pdirect = p_direct(light, nextHit.t, cosL)
                            * shadow_hit(light, world, ray);
            float pmetal = p_metal(reflection, newDir, hit.mat.fuzz);
            float weight = wmetal * pmetal / (wmetal * pmetal + wdirect * pdirect);

            ray.radiance += ray.throughput * Le * weight;
            hit.t = -2;
            return;
        }
    }
    else{
        // Direct lighting
        vec3 nLight = randomOnUnitHemiphere(seed, ray.origin - light.pos);
        vec3 lightPoint = light.pos + nLight * light.rad;
        vec3 newDir = normalize(lightPoint - ray.origin);
        ray.dir = newDir;

        float cosR = max(dot(hit.normal, newDir), 0.0);
        float cosL = max(dot(-newDir, nLight), 0);
        vec3 Le = light.color * light.intensity;
        
        float distance = length(lightPoint - ray.origin);
        float pdirect = p_direct(light, distance, cosL);
        ray.throughput *= lambert(hit) * cosR;

        if (shadow_hit(light, world, ray) > 0){
            float pmetal = p_metal(reflection, newDir, hit.mat.fuzz);
            float weight = 1.0 / (wdirect * pdirect + wmetal * pmetal);

            ray.radiance += ray.throughput * Le * weight;
        }
        else{
            ray.throughput /= pdirect * wdirect;
        }
        hit.t = -2;
    }
}

void emit(inout Hit hit, inout Ray ray){
    ray.radiance = hit.mat.color * hit.mat.intensity;
    hit.t = -1;
}

void computeLighting(World world, in out Hit hit, in out Ray ray, in out uint seed){
    if (hit.t < 0) return;

    switch (hit.mat.type){
        case MAT_DIFF:
            diffuse(world, hit, ray, seed);
            break;
        case MAT_METAL:
            metal(world, hit, ray, seed);
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



vec4 rayColor(World world, in out uint seed, Ray ray){
    Ray tracedRay = ray;
    for (int i = 0; i < MAX_BOUNCES; i++){

        Hit hit = rayIntersection(world, tracedRay);
        computeLighting(world, hit, tracedRay, seed);

        if (hit.t < 0){
            //tracedRay.radiance += tracedRay.throughput * sky(tracedRay.dir);
            return vec4(tracedRay.radiance, 1);
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
        camera.pos, 
        camera.lookDir, 
        vec3(1),
        vec3(0)
    );
    ray = fovRay(pos, ray);

    Mat sphereMat = Mat(MAT_DIFF, vec3(1, 0, 0), 1, 0);
    Mat sphereMat2 = Mat(MAT_METAL, vec3(1.), 1, 1);
    Mat planeMat = Mat(MAT_DIFF, vec3(0.7), 1, 0);

    Sphere spheres[NUM_SPHERE];
    spheres[0] = Sphere(vec3(0,EPS,-5), 2, sphereMat);
    spheres[1] = Sphere(vec3(0,EPS,-10), 2, planeMat);
    spheres[2] = Sphere(vec3(6, 2-EPS, -5), 4, sphereMat2);
    
    Plane planes[NUM_PLANE];
    planes[0] = Plane(vec3(0,-2,-3), vec3(0,1,0), planeMat);
    Light lights[NUM_LIGHT];
    lights[0] = Light(vec3(0,3,2), 1.5, vec3(1), 15);
    World world = World(spheres, planes, lights);

    vec4 radiance;
    for(int i = 0; i < samples; i++){
        ray.origin = ray.origin + 1.5*randomInSphere(seed) / winSize.y; // slight Anti-aliasing
        radiance += max(rayColor(world, seed, ray), 0);
    }

    FragColor = radiance + float(frameCount - samples) * texture(screenTex, uv);
    FragColor /= frameCount;
}