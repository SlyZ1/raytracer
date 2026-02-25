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
        int index_offset = 0;
        for (int i = 0; i < size; i++)
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

AABB Mesh::triangleBounds(const Triangle& tri) {
    AABB box;
    box.min =  {std::min({tri.v3.x, tri.v1.x, tri.v2.x}),
                std::min({tri.v3.y, tri.v1.y, tri.v2.y}),
                std::min({tri.v3.z, tri.v1.z, tri.v2.z})};
    box.max =  {std::max({tri.v3.x, tri.v1.x, tri.v2.x}),
                std::max({tri.v3.y, tri.v1.y, tri.v2.y}),
                std::max({tri.v3.z, tri.v1.z, tri.v2.z})};
    return box;
}

AABB Mesh::computeBounds(const vector<Triangle>& triangles, int start, int end) {
    AABB box = triangleBounds(triangles[start]);
    for(int i = start + 1; i < end; ++i)
        box.expand(triangleBounds(triangles[i]));
    return box;
}

BVHNode* Mesh::computeBVH(vector<Triangle>& triangles, 
                          vector<int>& indices,
                          int begin, int end) {
    BVHNode* node = new BVHNode();

    node->bounds = triangleBounds(triangles[indices[begin]]);
    for (int i = begin + 1; i < end; ++i)
        node->bounds.expand(triangleBounds(triangles[indices[i]]));

    const int MAX_TRIANGLES_PER_LEAF = 1;
    if (end - begin <= MAX_TRIANGLES_PER_LEAF) {
        node->triangle = indices[begin]; 
        return node;
    }

    AABB centroidBounds;
    centroidBounds.min = centroidBounds.max = triangles[indices[begin]].centroid();
    for (int i = begin + 1; i < end; ++i)
        centroidBounds.expand(triangles[indices[i]].centroid());

    glm::vec3 extent = centroidBounds.max - centroidBounds.min;
    int axis = 0;
    if (extent.y > extent.x && extent.y > extent.z)
        axis = 1;
    else if (extent.z > extent.x)
        axis = 2;

    std::sort(indices.begin() + begin, indices.begin() + end,
        [&triangles, axis](int a, int b) {
            return triangles[a].centroid()[axis] < triangles[b].centroid()[axis];
        });

    int mid = begin + (end - begin) / 2;

    node->left = computeBVH(triangles, indices, begin, mid);
    node->right = computeBVH(triangles, indices, mid, end);

    return node;
}

int indexOfTriangle(const vector<Triangle>& triangles, const Triangle& triangle){
    int size = triangles.size();
    for(int i = 0; i < size; i++){
        Triangle tri = triangles[i];
        if (tri.v1 == triangle.v1 && tri.v2 == triangle.v2 && tri.v3 == triangle.v3){
            return i;
        }
    }
    return -1;
}

int lineariseRec(BVHNode* node, vector<linBVHNode>& nodes, const vector<Triangle>& triangles){
    if (node == nullptr){
        return -1;
    }

    int left = lineariseRec(node->left, nodes, triangles);
    int right = lineariseRec(node->right, nodes, triangles);
    int triangle = node->triangle;

    linBVHNode linNode = {
        node->bounds,
        left,
        right,
        triangle
    };
    int i = nodes.size();
    nodes.push_back(linNode);
    return i;
}

vector<linBVHNode> Mesh::lineariseBVH(BVHNode* node, const vector<Triangle>& triangles){
    vector<linBVHNode> nodes = {};
    lineariseRec(node, nodes, triangles);
    return nodes;
}