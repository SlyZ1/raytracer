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
#define MAX_BOUNCES 4
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

void createTangentBasis(in vec3 N, out vec3 T, out vec3 B) {
    vec3 up = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);

    T = normalize(cross(up, N));
    B = cross(N, T);
}

vec3 randomGGXHemisphere(inout uint seed, vec3 normal, float alpha){
    float xi1 = rand(seed);
    float xi2 = rand(seed);
    float phi = 2.0 * PI * xi2;
    float theta = atan(alpha * sqrt(xi1) / sqrt(1.0 - xi1));

    vec3 T, B;
    createTangentBasis(normal, T, B);
    vec3 H_local = vec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    vec3 H = normalize(H_local.x * T + H_local.y * B + H_local.z * normal);

    return H;
}

vec3 reflect(vec3 I, vec3 N) {
    return I - 2.0 * dot(I, N) * N;
}

void stop(inout Hit hit){
    hit.t = -2;
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

// brdfs

vec3 lambert(Hit hit){
    return hit.mat.color / PI;
}

float wrap(float ndotl, float wrapFactor){
    float wrap = max(0, (ndotl + wrapFactor) / (1 + wrapFactor));
    return wrap;
}

float cosPhiDiff(vec3 normal, vec3 lightDir, vec3 viewDir){
    vec3 Lp = normalize(lightDir - normal * dot(lightDir, normal));
    vec3 Vp = normalize(viewDir - normal * dot(viewDir, normal));

    return dot(Lp, Vp);
}

vec3 oren_nayar(Hit hit, vec3 normal, vec3 lightDir, vec3 viewDir, float roughness){
    return hit.mat.color / PI;
}

vec3 schlickFresnel(float VdotH, vec3 F0)
{
    return F0 + (vec3(1) - F0) * pow(1 - VdotH, 5);
}

float DGTR(float a2, float NdotH){
    float denom = NdotH * NdotH * (a2 - 1) + 1;
    return a2 / max(PI * denom * denom, EPS);
}

float G1GTR(float a2, float NdotW){
    float denom = NdotW + sqrt(a2 + (1.0 - a2) * NdotW * NdotW);
    return 2.0 * NdotW / max(denom, EPS);
}

vec3 cookTorrance(Hit hit, vec3 viewDir, vec3 lightDir, float roughness){
    vec3 h = normalize(lightDir + viewDir);
    float NdotH = max(dot(hit.normal, h), 0);
    float NdotL = max(dot(hit.normal, lightDir), 0);
    float NdotV = max(dot(viewDir, hit.normal), 0);
    float VdotH = max(dot(h, viewDir), 0);
    float alpha = roughness * roughness;
    float a2 = alpha * alpha;
    float D = DGTR(a2, NdotH);
    float G = G1GTR(a2, NdotV) * G1GTR(a2, NdotL);
    vec3 F = schlickFresnel(VdotH, hit.mat.color);

    return F * D * G / max(4 * NdotV * NdotL, EPS);
}

// pdfs

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

float p_GGX(vec3 normal, vec3 l, vec3 v, float alpha){
    vec3 h = normalize(l + v);
    float NdotH = max(dot(normal, h), 0.0);
    float VdotH = max(dot(v, h), 0.0);
    float D = DGTR(alpha * alpha, NdotH);
    return D * NdotH / max(4.0 * VdotH, EPS);
}

// lighting functions

void russianRoulette(inout Ray ray, inout Hit hit, inout uint seed){
    float prob = max(ray.throughput.r, ray.throughput.g);
    prob = max(ray.throughput.b, prob);
    if (rand(seed) > prob) {
        ray.throughput = vec3(0);
        hit.t = -1;
    }  
    else { ray.throughput /= prob; }
}

void diffuse(World world, inout Hit hit, inout Ray ray, inout uint seed){
    ray.origin = ray.origin + hit.t * ray.dir + EPS * hit.normal;

    Light light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
    vec3 nLight = randomOnUnitHemiphere(seed, ray.origin - light.pos);
    vec3 lightPoint = light.pos + nLight * light.rad;
    vec3 lightDir = normalize(lightPoint - ray.origin);

    // MIS
    float wdirect = bsdfWeighting / 8.0;
    float wbsdf = 1 - bsdfWeighting / 8.0;
    wdirect = 3/8.0;
    wbsdf = 5/8.0;
    float r = rand(seed);
    float orenNayarRoughness = sqrt(0.3);

    if (r < wbsdf){
        // BSDF sampling
        vec3 newDir = randomCosineHemisphere(seed, hit.normal, 1);
        vec3 f_r = vec3(0);
        if (bsdfType == BSDF_OREN_NAYAR) 
            f_r = oren_nayar(hit, hit.normal, newDir, ray.dir, orenNayarRoughness);
        if (bsdfType == BSDF_LAMBERT || bsdfType == BSDF_LAMBERT_WRAP) f_r = lambert(hit);
        ray.throughput *= f_r * PI / wbsdf;
        ray.dir = newDir;

        if (bsdfType == BSDF_LAMBERT_WRAP){
            float wrap = wrap(dot(hit.normal, newDir), 0.5);
            //float cos = max(dot(hit.normal, newDir), wrap);
            ray.throughput *= wrap / dot(newDir, hit.normal);
        }
        // Check if we hit a light with the BSDF sampling
        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * nextHit.mat.intensity;
            float LdotNl = max(dot(-newDir, nextHit.normal), 1);

            float pdirect = p_direct(light, nextHit.t, LdotNl)
                            * shadow_hit(light, world, ray);
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = wbsdf * pbsdf / (wbsdf * pbsdf + wdirect * pdirect);

            ray.radiance += ray.throughput * Le * weight;
            stop(hit);
        }
    } 
    else {
        // Direct lighting
        vec3 newDir = lightDir;
        float NdotL = max(dot(hit.normal, newDir), 0.0);
        float LdotNl = max(dot(-newDir, nLight), 0);
        vec3 f_r = vec3(0);
        if (bsdfType == BSDF_OREN_NAYAR) f_r = oren_nayar(hit, hit.normal, newDir, ray.dir, orenNayarRoughness);
        if (bsdfType == BSDF_LAMBERT || bsdfType == BSDF_LAMBERT_WRAP) f_r = lambert(hit);
        ray.throughput *= f_r * NdotL;
        ray.dir = newDir;

        vec3 Le = light.color * light.intensity;
        
        float distance = length(lightPoint - ray.origin);
        float pdirect = p_direct(light, distance, LdotNl);

        if (shadow_hit(light, world, ray) > 0){
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = 1.0 / (wdirect * pdirect + wbsdf * pbsdf);

            if (bsdfType == BSDF_LAMBERT_WRAP){
                ray.throughput *= wrap(dot(hit.normal, lightDir), 0.5) / NdotL;
            }
            ray.radiance += ray.throughput * Le * weight;
        }
        else{
            ray.throughput /= pdirect * wdirect;
        }
        stop(hit);
    }
    //russianRoulette(ray, hit, seed);
}

void metal(World world, inout Hit hit, in out Ray ray, in out uint seed){
    Light light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
    ray.origin = ray.origin + hit.t * ray.dir + hit.normal * EPS;
    
    if (hit.mat.fuzz < EPS){
        vec3 newDir = reflect(ray.dir, hit.normal);
        float VdotN = max(dot(-ray.dir, hit.normal), 0.0);
        vec3 f_r = schlickFresnel(VdotN, hit.mat.color);
        ray.throughput *= f_r;
        ray.dir = newDir;
        return;
    }

    // MIS
    float t = clamp(hit.mat.fuzz / 0.2, 0, 1);
    float wdirect = 3/8.0 * t;
    float wGGX = 5/8.0 + 3/8.0 * (1 - t);
    float r = rand(seed);
    if (r < wGGX){
        float alpha = hit.mat.fuzz * hit.mat.fuzz;
        vec3 h = randomGGXHemisphere(seed, hit.normal, alpha);
        vec3 newDir = reflect(ray.dir, h);
        
        float NdotL = max(dot(hit.normal, newDir), 0.0);
        float pGGX = p_GGX(hit.normal, newDir, -ray.dir, alpha);
        vec3 f_r = cookTorrance(hit, -ray.dir, newDir, hit.mat.fuzz);
        ray.throughput *= f_r * NdotL / (pGGX * wGGX);
        ray.dir = newDir;

        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * nextHit.mat.intensity;
            float LdotNl = max(dot(-newDir, nextHit.normal), 0);
            float pdirect = p_direct(light, nextHit.t, LdotNl)
                            * shadow_hit(light, world, ray);
            float weight = wGGX * pGGX / (wGGX * pGGX + wdirect * pdirect);

            ray.radiance += ray.throughput * Le * weight;
            stop(hit);
        }
    }
    else{
        // Direct lighting
        vec3 nLight = randomOnUnitHemiphere(seed, ray.origin - light.pos);
        vec3 lightPoint = light.pos + nLight * light.rad;
        vec3 newDir = normalize(lightPoint - ray.origin);
        vec3 viewDir = -ray.dir;
        ray.dir = newDir;

        float NdotL = max(dot(hit.normal, newDir), 0.0);
        float LdotNl = max(dot(-newDir, nLight), 0);
        vec3 Le = light.color * light.intensity;
        
        float distance = length(lightPoint - ray.origin);
        float pdirect = p_direct(light, distance, LdotNl);
        vec3 f_r = cookTorrance(hit, viewDir, newDir, hit.mat.fuzz);
        ray.throughput *= f_r * NdotL / (pdirect * wdirect);

        if (shadow_hit(light, world, ray) > 0){
            float pGGX = p_GGX(hit.normal, newDir, viewDir, hit.mat.fuzz * hit.mat.fuzz);
            float weight = (pdirect * wdirect) / (wdirect * pdirect + wGGX * pGGX);

            ray.radiance += ray.throughput * Le * weight;
        }
        stop(hit);
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
            return;
        case MAT_METAL:
            metal(world, hit, ray, seed);
            return;
        case MAT_EMIT:
            emit(hit, ray);
            return;
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
            tracedRay.radiance += tracedRay.throughput * sky(tracedRay.dir);
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
    Mat sphereMat2 = Mat(MAT_METAL, vec3(1.), 1, bsdfWeighting / 8.0);
    Mat planeMat = Mat(MAT_DIFF, vec3(0.7), 1, 0);

    Sphere spheres[NUM_SPHERE];
    spheres[0] = Sphere(vec3(0,EPS,-5), 2, sphereMat);
    spheres[1] = Sphere(vec3(0,EPS,-10), 2, planeMat);
    spheres[2] = Sphere(vec3(6, 2-EPS, -5), 4, sphereMat2);
    
    Plane planes[NUM_PLANE];
    planes[0] = Plane(vec3(0,-2,-3), vec3(0,1,0), planeMat);
    Light lights[NUM_LIGHT];
    lights[0] = Light(vec3(0,5,10), 2.3, vec3(1), 10);
    World world = World(spheres, planes, lights);

    vec4 radiance;
    for(int i = 0; i < samples; i++){
        ray.origin = ray.origin + 1.5*randomInSphere(seed) / winSize.y; // slight Anti-aliasing
        radiance += max(rayColor(world, seed, ray), 0);
    }

    FragColor = radiance + float(frameCount - samples) * texture(screenTex, uv);
    FragColor /= frameCount;
}