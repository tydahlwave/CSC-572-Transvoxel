//
//  Scene.h
//  LOGLEngine
//
//  Created by Tyler Dahl on 5/28/17.
//
//

#ifndef Scene_h
#define Scene_h

#include <vector>
#include <iostream>

#include "Voxel.h"

class Scene {
public:
    int size;
    
    Scene(int size) :size(size) {
        // Create 3D volume initialized to 0s
        Voxel example;//; = {0, {glm::vec3(0)}, {glm::vec3(0)}};
//        example.isovalue = 0;
//        for (int i = 0; i < 15; i++) {
//            example.vertices[i] = glm::vec3(0);
//            example.normals[i] = glm::vec3(0);
//        }
        volumeData.resize(size*size*size, example);
    }
    
//private:
    std::vector<Voxel> volumeData;
};

#endif /* Scene_h */
