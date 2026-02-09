#include "mesh.hpp"
#include <tinyobjloader/tiny_obj_loader.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

Mesh::Mesh() {}

void Mesh::loadFromModel(const char* path){
    ifstream ifs(path);
    if (!ifs) {
        cerr << "Cannot open OBJ file\n" << endl;
        return;
    }

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    tinyobj::LoadObj(
        &attrib,
        &shapes,
        &materials,
        &warn,
        &err,
        &ifs,
        nullptr
    );

    if (!warn.empty()) cout << warn << endl;
    if (!err.empty()){
        cerr << err << endl;
        return;
    }
    
    auto getPos = [&](tinyobj::index_t i){
        return glm::vec3(
            attrib.vertices[3 * i.vertex_index + 0],
            attrib.vertices[3 * i.vertex_index + 1],
            attrib.vertices[3 * i.vertex_index + 2]
        );
    };
    
    for(const auto &shape : shapes){
        const auto &mesh = shape.mesh;
        int numVertex = 3;
        int size = mesh.num_face_vertices.size();
        size_t index_offset = 0;
        for (size_t i = 0; i < size; i++)
        {
            int fv = mesh.num_face_vertices[i];
            if (fv == 3){
                auto i0 = mesh.indices[numVertex * i + 0];
                auto i1 = mesh.indices[numVertex * i + 1];
                auto i2 = mesh.indices[numVertex * i + 2];
    
                Triangle newTriangle = Triangle();
                newTriangle.v1 = getPos(i0);
                newTriangle.v2 = getPos(i1);
                newTriangle.v3 = getPos(i2);
                
                m_triangles.push_back(newTriangle);   
            }

            index_offset += fv;
        }
    }
}

vector<Triangle> Mesh::getTriangles() const {
    return m_triangles;
}

AABB Mesh::triangleBounds(const Triangle& tri) const {
    AABB box;
    box.min =  {std::min({tri.v3.x, tri.v1.x, tri.v2.x}),
                std::min({tri.v3.y, tri.v1.y, tri.v2.y}),
                std::min({tri.v3.z, tri.v1.z, tri.v2.z})};
    box.max =  {std::max({tri.v3.x, tri.v1.x, tri.v2.x}),
                std::max({tri.v3.y, tri.v1.y, tri.v2.y}),
                std::max({tri.v3.z, tri.v1.z, tri.v2.z})};
    return box;
}

AABB Mesh::computeBounds(const std::vector<Triangle>& triangles, int start, int end) const {
    AABB box = triangleBounds(triangles[start]);
    for(int i = start + 1; i < end; ++i)
        box.expand(triangleBounds(triangles[i]));
    return box;
}

BVHNode* Mesh::computeBVH(const std::vector<Triangle>& triangles, int start, int end) const {
    BVHNode* node = new BVHNode();
    return node;
}