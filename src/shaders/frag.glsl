#version 430 core

struct Camera {
    vec3 pos;
    vec3 lookDir;
};

struct Mat {
    int type;
    vec3 color;
    float[2] data;
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
    bool inside;
};

#define NUM_SPHERE 6
#define NUM_PLANE 1
#define NUM_LIGHT 1

struct World {
    Sphere spheres[NUM_SPHERE];
    Plane planes[NUM_PLANE];
#if NUM_LIGHT > 0
    Light lights[NUM_LIGHT];
#endif
};

struct RaycastData {
    Hit hit;
    Ray ray;
    uint seed;
};

out vec4 FragColor;
layout (location = 1) uniform float time;
layout (location = 2) uniform vec2 texSize;
layout (location = 3) uniform sampler2D screenTex;
layout (location = 4) uniform int frameCount;
layout (location = 5) uniform int bsdfType;
layout (location = 6) uniform int bsdfWeighting;
layout (location = 7) uniform int samples;
layout (location = 8) uniform vec2 winSize;
layout (location = 9) uniform int maxBounces;
uniform Camera camera;
in vec4 vClipPos;

//#define SAMPLES 1
#define EPS 1e-4
#define PROBA_EPS 1e-6
#define MAX_NONEMIT_BOUNCE 10
#define PI 3.14159265

#define MAT_DIFF 0
#define MAT_PBR 1
#define MAT_EMIT 2
#define MAT_GLASS 3

#define updateData(data) data = RaycastData(hit, ray, seed)
#define unwrapData(data) ray = data.ray; hit = data.hit; seed = data.seed

#define mData(d0, d1) float[2](d0, d1)
#define mData0(d) float[2](d, 0)
#define mData1(d) float[2](0, d)
#define mNoData() float[2](0, 0)

#define pbrFuzz(m) m.data[0]              
#define pbrMetallic(m) m.data[1]
#define emitIntensity(m) m.data[0]          
#define glassIndex(m) m.data[0]

#define BSDF_LAMBERT 0
#define BSDF_LAMBERT_WRAP 1
#define BSDF_OREN_NAYAR 2

#define QON 0
#define FON 1
#define OREN_NAYAR_APPROX 1

// -------------------- UTILS


vec2 ratio(vec2 vec){
    return vec2(vec.x * texSize.x / texSize.y, vec.y);
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
    vec3 H_local = vec3(
        sin(theta) * cos(phi), 
        sin(theta) * sin(phi), 
        cos(theta)
    );

    vec3 T, B;
    createTangentBasis(normal, T, B);
    vec3 H = normalize(H_local.x * T + H_local.y * B + H_local.z * normal);

    return H;
}


vec3 reflect(vec3 I, vec3 N) {
    return I - 2.0 * dot(I, N) * N;
}

vec3 refract(vec3 I, vec3 N, float n) {
    vec3 r_perp = n * (I - N * dot(N, I));
    vec3 r_para = - sqrt(abs(1 - dot(r_perp, r_perp))) * N;
    return normalize(r_perp + r_para);
}

float luminanceMean(vec3 c){
    return dot(c, vec3(0.2126, 0.7152, 0.0722));
}

void stop(inout Hit hit, bool touchedLight){
    hit.t = touchedLight ? -2 : -1;
}


// -------------------- INTERSECTIONS

Hit sphereIntersect(Sphere sphere, Ray ray){
    vec3 oc = ray.origin - sphere.pos;
    float b = dot(oc, ray.dir);
    float c = dot(oc, oc) - sphere.rad * sphere.rad;
    float h = b*b - c;
    Hit hit;
    hit.t = -1;
    if (h < 0.) return hit;
    bool inside = false;
    float sqrtH = sqrt(h);
    float t = -b - sqrtH;
    if (t <= 0){
        t = -b + sqrtH;
        inside = true;
        if (t <= 0) return hit;
    }
    vec3 pos = ray.origin + t * ray.dir;
    vec3 normal = normalize(pos - sphere.pos);
    return Hit(t, normal, sphere.mat, inside);
}

Hit lightIntersect(Light light, Ray ray){
    vec3 oc = ray.origin - light.pos;
    float b = dot(oc, ray.dir);
    float c = dot(oc, oc) - light.rad * light.rad;
    float h = b*b - c;
    Hit hit;
    hit.t = -1;
    if (h < 0.) return hit;
    float t = -b - sqrt(h);
    if (t <= 0) return hit;
    vec3 pos = ray.origin + t * ray.dir;
    Mat mat = Mat(MAT_EMIT, light.color, mData0(light.intensity));
    return Hit(t, normalize(pos - light.pos), mat, false);
}

Hit planeIntersect(Plane plane, Ray ray){
    vec3 rp = plane.origin - ray.origin;
    float t = dot(rp, plane.normal) / dot(ray.dir, plane.normal);
    vec3 relativePoint = ray.origin + t * ray.dir;
    vec3 difference = relativePoint - plane.origin;
    Hit hit;
    hit.t = -1;
    if (t <= 0 || dot(difference, difference) > 10000) 
        return hit;

    if (int((abs(relativePoint.x + 1) * 0.5) + int(abs(relativePoint.z + 1) * 0.5 + 1)) % 2 == 0) 
        plane.mat = Mat(plane.mat.type, vec3(0.83), mNoData());

    return Hit(t, plane.normal, plane.mat, false);
}

Hit rayIntersection(World world, Ray ray){
    Hit hit;
    hit.t = 100000;
    for(int i = 0; i < NUM_SPHERE; i += 1){
        if (dot(world.spheres[i].pos - ray.origin, ray.dir) < 0 && false) continue;
        Hit sphereHit = sphereIntersect(world.spheres[i], ray);
        if (sphereHit.t > 0 && sphereHit.t < hit.t) hit = sphereHit;
    }
    for(int i = 0; i < NUM_PLANE; i += 1){
        if (dot(world.planes[i].normal, ray.dir) >= 0) continue;
        Hit planeHit = planeIntersect(world.planes[i], ray);
        if (planeHit.t > 0 && planeHit.t < hit.t) hit = planeHit;
    }
#if NUM_LIGHT > 0
    for(int i = 0; i < NUM_LIGHT; i += 1){
        Hit lightHit = lightIntersect(world.lights[i], ray);
        if (lightHit.t > 0 && lightHit.t < hit.t) hit = lightHit;
    }
#endif

    if(hit.t == 100000) hit.t = -1;

    return hit;
}

// -------------------- LIGHTS

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
    vec3 lambert = hit.mat.color / PI;

    float NdotL = max(dot(normal, lightDir), 0);
    float NdotV = max(dot(normal, viewDir), 0);

#if OREN_NAYAR_APPROX == QON
    float sigma2 = pow(roughness * PI * 0.5, 2);
    float A = 1 - 0.5 * sigma2 / (sigma2 + 0.33);
    float B = 0.45 * sigma2 / (sigma2 + 0.09);

    float cosPhi = max(0, cosPhiDiff(normal, lightDir, viewDir));
    float s = cosPhi * sqrt(1 - NdotL * NdotL) * sqrt(1 - NdotV * NdotV);
    float oneOverTq = s <= 0 ? 0 : 1 / max(NdotL, NdotV);
    float angleTerm = A + B * s * oneOverTq;

#elif OREN_NAYAR_APPROX == FON
    float A = 1 / (1 + (0.5 - 2 / (3 * PI)) * roughness);
    float B = roughness * A;

    float cosPhi = max(0, cosPhiDiff(normal, lightDir, viewDir));
    float s = cosPhi * sqrt(1 - NdotL * NdotL) * sqrt(1 - NdotV * NdotV);
    float tq = s <= 0 ? 1 : max(NdotL, NdotV);
    float angleTerm = A + B * s / tq;
#endif

    return lambert * angleTerm;
}

vec3 schlickFresnel(float VdotH, vec3 F0)
{
    return F0 + (vec3(1) - F0) * pow(1 - VdotH, 5);
}

float DGTR(float a2, float NdotH){
    float denom = NdotH * NdotH * (a2 - 1) + 1;
    return a2 / max(PI * denom * denom, PROBA_EPS);
}

float G1GTR(float a2, float NdotW){
    if (NdotW <= 0.0) return 0.0;

    float NdotW2 = NdotW * NdotW;
    float tan2 = (1.0 - NdotW2) / max(NdotW2, PROBA_EPS);
    return 2.0 / (1.0 + sqrt(1.0 + a2 * tan2));
}

vec3 cookTorrance(Hit hit, vec3 viewDir, vec3 lightDir, float alpha, bool simplify){
    vec3 h = normalize(lightDir + viewDir);
    float NdotL = max(dot(hit.normal, lightDir), 0);
    float NdotV = max(dot(viewDir, hit.normal), 0);
    float VdotH = max(dot(h, viewDir), 0);
    float a2 = alpha * alpha;
    float D = 1; // D is simplified with the PDF
    if (!simplify){
        float NdotH = max(dot(hit.normal, h), 0);
        D = DGTR(a2, NdotH);
    }
    float G = G1GTR(a2, NdotV) * G1GTR(a2, NdotL);
    vec3 F = schlickFresnel(VdotH, hit.mat.color);

    return F * D * G / max(4 * NdotV * NdotL, PROBA_EPS);
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

float p_GGX(vec3 normal, vec3 l, vec3 v, float alpha, bool simplify){
    vec3 h = normalize(l + v);
    float NdotH = max(dot(normal, h), 0.0);
    float VdotH = max(dot(v, h), 0.0);
    float D = 1;
    if (!simplify) D = DGTR(alpha * alpha, NdotH);
    return D * NdotH / max(4.0 * VdotH, PROBA_EPS);
}

// lighting functions

void russianRoulette(inout RaycastData data){
    Ray ray; Hit hit; uint seed;
    unwrapData(data);
    
    float prob = max(max(ray.throughput.r, ray.throughput.g), ray.throughput.b);
    if (rand(seed) > prob) {
        stop(hit, false);
    }  
    else { ray.throughput /= max(min(prob, 1), EPS); }

    updateData(data);
}

float sampleLight(inout RaycastData data, Light light){
    Ray ray; Hit hit; uint seed;
    unwrapData(data);

    vec3 nLight = randomOnUnitHemiphere(seed, ray.origin - light.pos);
    vec3 lightPoint = light.pos + nLight * light.rad;
    vec3 newDir = normalize(lightPoint - ray.origin);
    ray.dir = newDir;

    float NdotL = max(dot(hit.normal, newDir), 0.0);
    float LdotNl = max(dot(-newDir, nLight), 0);
    
    float distance = length(lightPoint - ray.origin);
    float pdirect = p_direct(light, distance, LdotNl);
    ray.throughput *= NdotL;

    updateData(data);
    return pdirect;
}

void diffuse(World world, inout RaycastData data){
    Ray ray; Hit hit; uint seed;
    unwrapData(data);
    ray.origin += hit.t * ray.dir + EPS * hit.normal;

    Light light;
#if NUM_LIGHT > 0
    light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
#endif

    // MIS
    float wdirect = bsdfWeighting / 8.0;
    float wbsdf = 1 - bsdfWeighting / 8.0;
    wdirect = 3/8.0;
    wbsdf = 5/8.0;
    wdirect = NUM_LIGHT > 0 ? wdirect : 0;
    wbsdf = NUM_LIGHT > 0 ? wbsdf : 1;
    float r = rand(seed);
    float orenNayarRoughness = bsdfWeighting / 8.0;
    vec3 viewDir = -ray.dir;

    if (r < wbsdf){
        // BSDF sampling
        vec3 newDir = randomCosineHemisphere(seed, hit.normal, 1);
        vec3 f_r = vec3(0);
        if (bsdfType == BSDF_OREN_NAYAR) 
            f_r = oren_nayar(hit, hit.normal, newDir, viewDir, orenNayarRoughness);
        if (bsdfType == BSDF_LAMBERT || bsdfType == BSDF_LAMBERT_WRAP) f_r = lambert(hit);
        ray.throughput *= f_r * PI / wbsdf;
        ray.dir = newDir;

        // Check if we hit a light with the BSDF sampling
        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * emitIntensity(nextHit.mat);
            float LdotNl = max(dot(-newDir, nextHit.normal), 1);

            float pdirect = p_direct(light, nextHit.t, LdotNl)
                            * shadow_hit(light, world, ray);
            float pbsdf = p_bsdf(hit.normal, newDir);
            float weight = wbsdf * pbsdf / (wbsdf * pbsdf + wdirect * pdirect);

            ray.radiance += ray.throughput * Le * weight;
            stop(hit, true);
        }
        else{
            ray.throughput = clamp(ray.throughput, 0, MAX_NONEMIT_BOUNCE);
        }
    } 
    else {
        // Direct lighting
        updateData(data);
        float pdirect = sampleLight(data, light);
        unwrapData(data);
        vec3 f_r = vec3(0);
        if (bsdfType == BSDF_OREN_NAYAR) f_r = oren_nayar(hit, hit.normal, ray.dir, viewDir, orenNayarRoughness);
        if (bsdfType == BSDF_LAMBERT || bsdfType == BSDF_LAMBERT_WRAP) f_r = lambert(hit);
        ray.throughput *= f_r;

        if (shadow_hit(light, world, ray) > 0){
            float pbsdf = p_bsdf(hit.normal, ray.dir);
            float weight = 1.0 / (wdirect * pdirect + wbsdf * pbsdf);
            vec3 Le = light.intensity * light.color;
            ray.radiance += ray.throughput * Le * weight;
            stop(hit, true);
        }
        else{
            ray.throughput /= pdirect * wdirect;
            stop(hit, false);
        }
    }

    updateData(data);
    russianRoulette(data);
}

void metal(World world, inout RaycastData data){
    Ray ray; Hit hit; uint seed;
    unwrapData(data);
    Light light;
#if NUM_LIGHT > 0
    light = world.lights[int(rand(seed) * (NUM_LIGHT - 1))];
#endif
    ray.origin += hit.t * ray.dir + hit.normal * EPS;
    
    if (pbrFuzz(hit.mat) < EPS){
        vec3 newDir = reflect(ray.dir, hit.normal);
        float VdotN = max(dot(-ray.dir, hit.normal), 0.0);
        vec3 f_r = schlickFresnel(VdotN, hit.mat.color);
        ray.throughput *= f_r;
        ray.dir = newDir;
        data = RaycastData(hit, ray, seed);
        return;
    }

    // MIS
    float fuzz = pbrFuzz(hit.mat);
    float alpha = fuzz * fuzz;
    float wdirect = 3.0/8.0;
    float wGGX = 5.0/8.0;
    wdirect = NUM_LIGHT > 0 ? wdirect : 0;
    wGGX = NUM_LIGHT > 0 ? wGGX : 1;
    float r = rand(seed);
    if (r < wGGX){
        vec3 viewDir = -ray.dir;
        vec3 h = randomGGXHemisphere(seed, hit.normal, alpha);
        vec3 newDir = reflect(-viewDir, h);
        
        float NdotL = max(dot(hit.normal, newDir), 0.0);
        float pGGX = p_GGX(hit.normal, newDir, viewDir, alpha, true);
        vec3 f_r = cookTorrance(hit, viewDir, newDir, alpha, true);
        ray.throughput *= f_r * NdotL / max(pGGX * wGGX, EPS);
        ray.dir = newDir;

        Hit nextHit = rayIntersection(world, ray);
        if (nextHit.t > 0 && nextHit.mat.type == MAT_EMIT){
            vec3 Le = nextHit.mat.color * emitIntensity(nextHit.mat);
            float LdotNl = max(dot(-newDir, nextHit.normal), 0);
            float pdirect = p_direct(light, nextHit.t, LdotNl)
                            * shadow_hit(light, world, ray);
            float weight = wGGX * pGGX / max(wGGX * pGGX + wdirect * pdirect, EPS);

            ray.radiance += ray.throughput * Le * weight;
            stop(hit, true);
        }
        else{
            ray.throughput = clamp(ray.throughput, 0, MAX_NONEMIT_BOUNCE);
        }
    }
    else{
        // Direct lighting
        vec3 viewDir = -ray.dir;
        updateData(data);
        float pdirect = sampleLight(data, light);
        unwrapData(data);
        vec3 f_r = cookTorrance(hit, viewDir, ray.dir, alpha, false);
        ray.throughput *= f_r;

        if (shadow_hit(light, world, ray) > 0){
            float pGGX = p_GGX(hit.normal, ray.dir, viewDir, alpha, false);
            float weight = 1.0 / max(wdirect * pdirect + wGGX * pGGX, EPS);
            vec3 Le = light.color * light.intensity;
            ray.radiance += ray.throughput * Le * weight;
            stop(hit, true);
        }
        else{
            ray.throughput /= (wdirect * pdirect);
            stop(hit, false);
        }
    }

    updateData(data);
    russianRoulette(data);
}

void glass(inout RaycastData data){
    Ray ray; Hit hit; uint seed;
    unwrapData(data);

    float n = 1 / glassIndex(hit.mat);
    vec3 normal = hit.normal;
    if (hit.inside){
        n = 1 / n;
        normal *= -1;
    }

    float NdotV = dot(normal, -ray.dir);
    bool totalReflection = n * n * (1 - NdotV * NdotV) > 1;
    float f0 = (1 - n) / (1 + n);
    f0 = f0 * f0;
    float reflectance = schlickFresnel(NdotV, vec3(f0)).x;
    if (totalReflection || rand(seed) < reflectance){
        ray.origin += hit.t * ray.dir + EPS * normal;
        ray.dir = reflect(ray.dir, hit.normal);
    }
    else{
        ray.origin += hit.t * ray.dir - EPS * normal;
        ray.dir = refract(ray.dir, normal, n);
    }

    if (hit.inside){
        vec3 absorption = exp(-(vec3(1) - hit.mat.color) * hit.t); // Beer-Lambert
        ray.throughput *= absorption;
    }
    updateData(data);
}

void pbr(World world, inout RaycastData data){
    Ray ray; Hit hit; uint seed;
    unwrapData(data);

    float metallic = pbrMetallic(hit.mat);
    vec3 F0 = mix(vec3(0.04), hit.mat.color, metallic);
    float F = luminanceMean(schlickFresnel(dot(-ray.dir, hit.normal), F0));
    float pSpec = mix(F, 1, metallic);
    if (rand(seed) < pSpec){
        metal(world, data);
    }
    else{
        diffuse(world, data);
    }
}

void emit(inout RaycastData data){
    data.ray.radiance += data.ray.throughput * data.hit.mat.color * emitIntensity(data.hit.mat);
    data.hit.t = -1;
}

void computeLighting(World world, in out Hit hit, in out Ray ray, in out uint seed){
    if (hit.t < 0) return;

    RaycastData data = RaycastData(hit, ray, seed);
    switch (hit.mat.type){
        case MAT_DIFF:
            diffuse(world, data);
            break;
        case MAT_PBR:
            pbr(world, data);
            break;
        case MAT_GLASS:
            glass(data);
            break;
        case MAT_EMIT:
            emit(data);
            break;
    }
    unwrapData(data);
}

// RAY TRACING --------------------

Ray fovRay(vec2 pos, Ray ray){
    float fov = radians(mix(50, 90, 0 / 8.0));
    vec3 forward = normalize(camera.lookDir);

    vec3 worldUp = abs(forward.y) < 0.999
                 ? vec3(0,1,0)
                 : vec3(0,0,1);

    vec3 right = normalize(cross(forward, worldUp));
    vec3 up    = cross(right, forward);

    float tanHalfFov = tan(fov * 0.5);

    vec3 dir = forward + (right * pos.x + up * pos.y) * tanHalfFov;
    ray.dir = normalize(dir);
    return ray;
}

vec3 sun(vec3 lookingAt)
{
    vec3 sunDir = normalize(vec3(0.2, 0.6, 0.7));
    float sunAngularRadius = radians(4);
    float sunIntensity = 50.0;
    float cosAngle = dot(lookingAt, sunDir);
    float sun = smoothstep(
        cos(sunAngularRadius),
        cos(sunAngularRadius * 0.8),
        cosAngle
    );
    return vec3(sunIntensity) * sun;
}

vec3 horizon(vec3 lookingAt){
    float a = (dot(normalize(lookingAt), vec3(0,1,0)) + 1) * 0.5;
    vec3 top = vec3(0.32, 0.55, 0.78);
    vec3 bot = vec3(0.62, 0.72, 0.85);
    return mix(bot, top, a) * 0.7;
}

vec3 sky(vec3 lookingAt){
    return /*sun(lookingAt) + */horizon(lookingAt);
}

vec4 rayColor(World world, in out uint seed, Ray ray){
    Ray tracedRay = ray;
    for (int i = 0; i < maxBounces; i++){

        Hit hit = rayIntersection(world, tracedRay);
        computeLighting(world, hit, tracedRay, seed);

        if (hit.t < 0){
            if (hit.t > -2) tracedRay.radiance += tracedRay.throughput * sky(tracedRay.dir);
            return vec4(tracedRay.radiance, 1);
        }// Cas ou on tombe sur une light ou sur rien
    }

    return vec4(0);
}

void main()
{
    float resolutionFactor = winSize.x / texSize.x;
    uint seed = initSeed(uvec2(gl_FragCoord.xy * resolutionFactor), frameCount);
    vec2 pos = ratio(vClipPos.xy) * resolutionFactor + ratio(vec2(1)) * (resolutionFactor - 1);
    vec2 uv = (vClipPos.xy + vec2(1)) * 0.5;

    Ray ray = Ray(
        camera.pos, 
        camera.lookDir, 
        vec3(1),
        vec3(0)
    );
    ray = fovRay(pos, ray);

    Mat planeMat = Mat(MAT_DIFF,  vec3(0.36, 0.20, 0.09), mNoData());
    Mat sphereMat = Mat(MAT_DIFF, vec3(1, 0, 0), mNoData());
    Mat glassMat = Mat(MAT_GLASS, vec3(0.75, 0.90, 0.92), mData0(1.33));

    float metallic = 0.4;
    Mat metalMat1 = Mat(MAT_PBR, vec3(1.0, 0.4, 0.4), mData(0 / 8.0, metallic));
    Mat metalMat2 = Mat(MAT_PBR, vec3(1.0, 0.4, 0.4), mData(2 / 8.0, metallic));
    Mat metalMat3 = Mat(MAT_PBR, vec3(1.0, 0.4, 0.4), mData(4 / 8.0, metallic));
    Mat metalMat4 = Mat(MAT_PBR, vec3(1.0, 0.4, 0.4), mData(6 / 8.0, metallic));

    Sphere spheres[NUM_SPHERE];
    spheres[0] = Sphere(vec3(0,EPS,-5), 2, sphereMat);
    spheres[1] = Sphere(vec3(0,EPS-0.5,9), 1, glassMat);
    spheres[2] = Sphere(vec3(4, EPS, -5), 2, metalMat1);
    spheres[3] = Sphere(vec3(10, EPS, -5), 2, metalMat2);
    spheres[4] = Sphere(vec3(16, EPS, -5), 2, metalMat3);
    spheres[5] = Sphere(vec3(22, EPS, -5), 2, metalMat4);
    
    Plane planes[NUM_PLANE];
    planes[0] = Plane(vec3(0,-2,-3), vec3(0,1,0), planeMat);

#if NUM_LIGHT > 0
    Light lights[NUM_LIGHT];
    lights[0] = Light(vec3(0,5,10), 2.3, vec3(1), 10);
    World world = World(spheres, planes, lights);
#else
    World world = World(spheres, planes);
#endif

    vec4 radiance;
    for(int i = 0; i < samples; i++){
        ray.origin = ray.origin + 5 * randomInSphere(seed) / winSize.y; // slight Anti-aliasing
        radiance += max(rayColor(world, seed, ray), 0);
    }

    FragColor = radiance + max(frameCount - samples, 0) * texture(screenTex, uv);
    FragColor /= max(frameCount, samples);
}