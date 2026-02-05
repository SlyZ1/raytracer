#include <iostream>
#include <fstream>
#include <sstream>
#include <glad/glad.h>
#include <vector>

using namespace std;

class ShaderProgram {
    private:
        unsigned int shaderProgram = 0;
        vector<unsigned int> shaders = {};

        string getShaderSource(const char *path);

    public:
        ShaderProgram();
        unsigned int id();
        void create();
        void load(int type, const char *path);
        void link();
        void use();
        void destroy();

        template<typename T>
        static tuple<unsigned int, unsigned int, unsigned int> 
        addData(const vector<T>& data, const vector<unsigned int>& indices){
            unsigned int VBO, VAO, EBO;
            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);
            glGenBuffers(1, &EBO);

            glBindVertexArray(VAO);

            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(T), data.data(), GL_STATIC_DRAW);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

            return {VBO, VAO, EBO};
        }

        static void linkData(int numCoords, int typesize, int layout, int glType = GL_FLOAT, int normalize = GL_FALSE){
            glVertexAttribPointer(layout, numCoords, glType, normalize, numCoords * typesize, (void*)0);
            glEnableVertexAttribArray(layout);

            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);
        }
};