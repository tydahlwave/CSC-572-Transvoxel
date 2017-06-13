//
//  Voxel.h
//  LOGLEngine
//
//  Created by Tyler Dahl on 5/28/17.
//
//

#ifndef Voxel_h
#define Voxel_h

struct Voxel {
    char vertexCount = 0;
    float isovalue = 0;
//    glm::vec3 vertices[15] = {glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0)};
//    glm::vec3 normals[15] = {glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0),glm::vec3(0)};
};

struct VoxelCube {
    glm::vec3 pos;
    float isovalue;
};

#endif /* Voxel_h */
