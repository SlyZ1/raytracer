#ifndef MESH
#define MESH

#include <glm/glm.hpp>
#include <vector>

using namespace std;

struct Triangle {
    glm::vec3 v1; 
    float pad0;
    glm::vec3 v2;
    float pad1;
    glm::vec3 v3;
    float pad2;

    glm::vec3 centroid() const {
        return {(v1.x + v3.x + v2.x) / 3.0f,
                (v1.y + v3.y + v2.y) / 3.0f,
                (v1.z + v3.z + v2.z) / 3.0f};
    }
};

struct AABB {
    glm::vec3 min;
    float pad0;
    glm::vec3 max;
    float pad1;

    void expand(const glm::vec3& p) {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
        max.z = std::max(max.z, p.z);
    }

    void expand(const AABB& box) {
        expand(box.min);
        expand(box.max);
    }
};

struct BVHNode {
    AABB bounds;
    BVHNode* left = nullptr;
    BVHNode* right = nullptr;
    int triangle = -1;
};

struct linBVHNode {
    AABB bounds;
    int left = -1;
    int right = -1;
    int triangle = -1;
    int pad = -1;
};

class Mesh {
public:
    Mesh();
    void loadFromModel(const char* path);
    vector<Triangle> getTriangles() const;
    static BVHNode* computeBVH(vector<Triangle>& triangles, vector<int>& indices, int begin, int end);
    static vector<linBVHNode> lineariseBVH(BVHNode* node, const vector<Triangle>& triangles);

private:
    vector<Triangle> m_triangles = {};
    static AABB triangleBounds(const Triangle& tri);
    static AABB computeBounds(const vector<Triangle>& triangles, int start, int end);
};

#endif