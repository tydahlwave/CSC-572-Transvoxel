// Local Headers
#include "loglengine.hpp"
#include "Shader.hpp"
#include "Scene.h"
#include "MatrixStack.h"
#include "Noise/OpenSimplexNoise.hpp"
#include "MarchingCubes.h"
#include "Transvoxel.h"
#include "Isovalues.h"

// System Headers
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// Standard Headers
#include <cstdio>
#include <cstdlib>

GLuint VBO, VAO;

GLuint volumeVBO, volumeVAO;
std::shared_ptr<Shader> volumeShader;
std::vector<float> trianglesVector;
std::vector<float> normalsVector;
float t = 0;
int shapeIndex = 1;
int voxelScale = 2;
int editingIsovalue = 10;
float cameraSpeedFactor = 1.0f;
bool isRotating = false;
bool isWireFrame = false;
bool isCursorDisabled = true;

typedef struct {
    glm::vec3 pos;
    glm::vec3 lookAt;
    glm::vec3 up;
    float velW;
    float velU;
} Camera;

Camera camera;

//void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode) {
//    // When a user presses the escape key, we set the WindowShouldClose property to true,
//    // closing the application
//    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
//        glfwSetWindowShouldClose(window, GL_TRUE);
//    }
//}

void resetVolumeData(Scene &scene);
void setupVolumeTriangles(Scene &scene);
glm::vec3 interpolateNormal(Scene &scene, int x, int y, int z, int isorange, int corner1, int corner2);
Scene *tempScene;

typedef struct Ray {
    glm::vec3 origin, dir;
    glm::vec3 invdir;
    int sign[3];
    
    Ray(glm::vec3 origin, glm::vec3 dir) : origin(origin), dir(dir) {
        invdir = 1.0f / dir;
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
    }
} Ray;

typedef struct Cube {
    //    glm::vec3 center;
    //    glm::vec3 extents;
    glm::vec3 min;
    glm::vec3 max;
    glm::vec3 bounds[2];
    
    //    Cube(glm::vec3 center, glm::vec3 extents) :center(center), extents(extents) {}
    Cube(glm::vec3 min, glm::vec3 max) :min(min), max(max) {
        bounds[0] = min;
        bounds[1] = max;
    }
    
    //    glm::vec3 getMinExtents() { return center - extents; }
    //    glm::vec3 getMaxExtents() { return center + extents; }
    //    glm::vec3 getSize() { return extents * 2.0f; }
    
    // Credit: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
    bool intersects2(Ray r) {
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        
        tmin = (bounds[r.sign[0]].x - r.origin.x) * r.invdir.x;
        tmax = (bounds[1-r.sign[0]].x - r.origin.x) * r.invdir.x;
        tymin = (bounds[r.sign[1]].y - r.origin.y) * r.invdir.y;
        tymax = (bounds[1-r.sign[1]].y - r.origin.y) * r.invdir.y;
        
        if ((tmin > tymax) || (tymin > tmax))
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;
        
        tzmin = (bounds[r.sign[2]].z - r.origin.z) * r.invdir.z;
        tzmax = (bounds[1-r.sign[2]].z - r.origin.z) * r.invdir.z;
        
        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;
        
        return true;
    }
    
    bool intersects(Ray ray) {
        // Transform the ray
//        glm::vec4 rayPosTransformed = invTransform * Vector4(ray->position[0], ray->position[1], ray->position[2], 1);
//        glm::vec4 rayDirTransformed = invTransform * Vector4(ray->direction[0], ray->direction[1], ray->direction[2], 0);
//        glm::vec3 rayPosition(rayPosTransformed[0], rayPosTransformed[1], rayPosTransformed[2]);
//        glm::vec3 rayDirection(rayDirTransformed[0], rayDirTransformed[1], rayDirTransformed[2]);
        auto M = std::make_shared<MatrixStack>();
        M->rotate(t, glm::vec3(0, 1, 0));
        M->translate(glm::vec3(-VOLUME_SIZE/2.0f, -VOLUME_SIZE/2.0f, -VOLUME_SIZE/2.0f));
        glm::vec4 rayPosTransformed = glm::inverse(M->topMatrix()) * glm::vec4(ray.origin, 1);
        glm::vec4 rayDirTransformed = glm::inverse(M->topMatrix()) * glm::vec4(ray.dir, 0);
        glm::vec3 rayPosition(rayPosTransformed[0], rayPosTransformed[1], rayPosTransformed[2]);
        glm::vec3 rayDirection(rayDirTransformed[0], rayDirTransformed[1], rayDirTransformed[2]);
        
        // Switch the min and max coords if the ray is coming from the back side of the box
        glm::vec3 newMin = min;
        glm::vec3 newMax = max;
        for (int i = 0; i < 3; i++) {
            if (fabs(newMax[i] - rayPosition[i]) < fabs(newMin[i] - rayPosition[i])) {
                float temp = newMin[i];
                newMin[i] = newMax[i];
                newMax[i] = temp;
            }
        }
        
        // Short-circuits if the ray is parallel to the box
        int skipDim = -1;
        if (rayDirection[0] == 0) {
            if (rayPosition[0] < min[0] || rayPosition[0] > max[0]/* || rayPosition[1] >= max[1] || rayPosition[2] >= max[2]*/) return false;
            skipDim = 0;
        }
        if (rayDirection[1] == 0) {
            if (rayPosition[1] < min[1] || rayPosition[1] > max[1]/* || rayPosition[0] >= max[0] || rayPosition[2] >= max[2]*/) return false;
            skipDim = 1;
        }
        if (rayDirection[2] == 0) {
            if (rayPosition[2] < min[2] || rayPosition[2] > max[2]/* || rayPosition[0] >= max[0] || rayPosition[1] >= max[1]*/) return false;
            skipDim = 2;
        }
        
        // Calculate intersections with the box
        float tgmin = FLT_MIN;
        float tgmax = FLT_MAX;
        float normalDimMin = -1;
        float normalDimMax = -1;
        // Check every dimension
        for (int i = 0; i < 3; i++) {
            if (i == skipDim) continue;
            float t1 = (newMin[i] - rayPosition[i]) / rayDirection[i];
            float t2 = (newMax[i] - rayPosition[i]) / rayDirection[i];
            if (t2 < t1) {
                float temp = t1;
                t1 = t2;
                t2 = temp;
            }
            if (t1 > tgmin) {
                tgmin = t1;
                normalDimMin = i;
            }
            if (t2 < tgmax) {
                tgmax = t2;
                normalDimMax = i;
            }
            if (tgmin > tgmax) return false; // no hit
            if (tgmax < 0.0001) return false; // hit is behind start of ray
        }
        return true;
    }
} Cube;

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        } else if (key == GLFW_KEY_SPACE) {
            shapeIndex += 1;
            resetVolumeData(*tempScene);
        } else if (key == GLFW_KEY_R) {
            isRotating = !isRotating;
        } else if (key == GLFW_KEY_T) {
            isWireFrame = !isWireFrame;
            glPolygonMode(GL_FRONT_AND_BACK, isWireFrame ? GL_LINE : GL_FILL);
        } else if (key == GLFW_KEY_M) {
            isCursorDisabled = !isCursorDisabled;
            glfwSetInputMode(window, GLFW_CURSOR, isCursorDisabled ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
        } else if (key == GLFW_KEY_I) {
            voxelScale = max(voxelScale / 2, 1);
            setupVolumeTriangles(*tempScene);
        } else if (key == GLFW_KEY_O) {
            voxelScale *= 2;
            setupVolumeTriangles(*tempScene);
        } else if (key == GLFW_KEY_Q) {
            editingIsovalue = -editingIsovalue;
        } else if (key == GLFW_KEY_W) {
            camera.velW = -CAMERA_SPEED * cameraSpeedFactor;
        } else if (key == GLFW_KEY_S) {
            camera.velW = CAMERA_SPEED * cameraSpeedFactor;
        } else if (key == GLFW_KEY_A) {
            camera.velU = -CAMERA_SPEED * cameraSpeedFactor;
        } else if (key == GLFW_KEY_D) {
            camera.velU = CAMERA_SPEED * cameraSpeedFactor;
        } else if (key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT) {
            cameraSpeedFactor = 3.0f;
            if (abs(camera.velU) > CAMERA_STOPPED_THRESHOLD)
                camera.velU = camera.velU / abs(camera.velU) * CAMERA_SPEED * cameraSpeedFactor;
            if (abs(camera.velW) > CAMERA_STOPPED_THRESHOLD)
                camera.velW = camera.velW / abs(camera.velW) * CAMERA_SPEED * cameraSpeedFactor;
        }
    } else if (action == GLFW_RELEASE) {
        if (key == GLFW_KEY_W) {
            if (camera.velW < -CAMERA_STOPPED_THRESHOLD) {
                camera.velW = 0;
            } else {
                camera.velW = CAMERA_SPEED * cameraSpeedFactor;
            }
        } else if (key == GLFW_KEY_S) {
            if (camera.velW > CAMERA_STOPPED_THRESHOLD) {
                camera.velW = 0;
            } else {
                camera.velW = -CAMERA_SPEED * cameraSpeedFactor;
            }
        } else if (key == GLFW_KEY_A) {
            if (camera.velU < -CAMERA_STOPPED_THRESHOLD) {
                camera.velU = 0;
            } else {
                camera.velU = CAMERA_SPEED * cameraSpeedFactor;
            }
        } else if (key == GLFW_KEY_D) {
            if (camera.velU > CAMERA_STOPPED_THRESHOLD) {
                camera.velU = 0;
            } else {
                camera.velU = -CAMERA_SPEED * cameraSpeedFactor;
            }
        } else if (key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT) {
            cameraSpeedFactor = 1.0f;
            if (abs(camera.velU) > CAMERA_STOPPED_THRESHOLD)
                camera.velU = camera.velU / abs(camera.velU) * CAMERA_SPEED * cameraSpeedFactor;
            if (abs(camera.velW) > CAMERA_STOPPED_THRESHOLD)
                camera.velW = camera.velW / abs(camera.velW) * CAMERA_SPEED * cameraSpeedFactor;
        }
    }
}

float alpha = 0;
float beta = -M_PI/2;
static void mouse_move_callback(GLFWwindow *window, double posX, double posY) {
    
    // Get current window size.
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    
    // Compute new alpha and beta for the camera lookat point
    double alpha2 = alpha + ((posY + height/2.0) / height * M_PI*16/18) - M_PI*8/18;
    double beta2 = beta + ((posX + width/2.0) / width * M_PI*2) - M_PI;
    
    // Constrain the view (up and down constrained to (-80,80) degrees)
    if (alpha2 > M_PI*8/18) alpha2 = M_PI*8/18;
    if (alpha2 < -M_PI*8/18) alpha2 = -M_PI*8/18;
    
    // Compute camera lookat point
    camera.lookAt[0] = camera.pos[0] + cos(alpha2) * cos(beta2);
    camera.lookAt[1] = camera.pos[1] + -sin(alpha2);
    camera.lookAt[2] = camera.pos[2] + cos(alpha2) * cos(M_PI/2 - beta2);
}

//void recomputeVertices(int x, int y, int z) {
//    int index = (x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z);
//    for (int i = 0; i < 15; i++) {
//        tempScene->volumeData[index].vertices[i] = glm::vec3(0);
//        tempScene->volumeData[index].normals[i] = glm::vec3(0);
//    }
//    
//    glm::vec3 cornerPositions[8] = {
//        glm::vec3(x-1, y-1, z-1),
//        glm::vec3(x,   y-1, z-1),
//        glm::vec3(x-1, y,   z-1),
//        glm::vec3(x,   y,   z-1),
//        glm::vec3(x-1, y-1, z),
//        glm::vec3(x,   y-1, z),
//        glm::vec3(x-1, y,   z),
//        glm::vec3(x,   y,   z)
//    };
//    
//    int cornerIndices[8] = {
//        (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z-1),
//        (x  )*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z-1),
//        (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z-1),
//        (x  )*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z-1),
//        (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z  ),
//        (x  )*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z  ),
//        (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z  ),
//        (x  )*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z  )
//    };
//    
//    unsigned long cubeIndex = 0;
//    float isorange = 0;
//    if ((tempScene->volumeData[cornerIndices[0]].isovalue) <= isorange) cubeIndex |= 1;
//    if ((tempScene->volumeData[cornerIndices[1]].isovalue) <= isorange) cubeIndex |= 2;
//    if ((tempScene->volumeData[cornerIndices[2]].isovalue) <= isorange) cubeIndex |= 4;
//    if ((tempScene->volumeData[cornerIndices[3]].isovalue) <= isorange) cubeIndex |= 8;
//    if ((tempScene->volumeData[cornerIndices[4]].isovalue) <= isorange) cubeIndex |= 16;
//    if ((tempScene->volumeData[cornerIndices[5]].isovalue) <= isorange) cubeIndex |= 32;
//    if ((tempScene->volumeData[cornerIndices[6]].isovalue) <= isorange) cubeIndex |= 64;
//    if ((tempScene->volumeData[cornerIndices[7]].isovalue) <= isorange) cubeIndex |= 128;
//    
//    /* Cube is entirely in/out of the surface */
//    if (edgeTable[cubeIndex] == 0) return;
//    
//    /* Transvoxel */
//    auto cellClass = regularCellClass[cubeIndex];
//    auto cellData = regularCellData[cellClass];
//    unsigned short vertexData[12] = {0};
//    for (int i = 0; i < 12; i++) {
//        vertexData[i] = regularVertexData[cubeIndex][i];
//    }
//    
//    /* Find the vertices where the surface intersects the cube for each edge */
//    std::vector<glm::vec3> edgeVertices;
//    std::vector<glm::vec3> edgeNormals;
//    for (auto vertex : vertexData) {
//        if (!vertex) break;
//        unsigned char corner1 = (vertex >> 4) & 0x000F;
//        unsigned char corner2 = vertex & 0x000F;
//        float isovalue1 = tempScene->volumeData[cornerIndices[corner1]].isovalue;
//        float isovalue2 = tempScene->volumeData[cornerIndices[corner2]].isovalue;
//        float t = isovalue2 / (isovalue2 - isovalue1);
//        glm::vec3 vertexPos = cornerPositions[corner1] * t + cornerPositions[corner2] * (1-t);
//        edgeVertices.push_back(vertexPos);
//        edgeNormals.push_back(interpolateNormal(*tempScene, x, y, z, isorange, corner1, corner2));
//    }
//    
//    /* Create the triangles */
//    int voxelIndex = x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z;
//    int indexCount = 0;
//    for (auto vertexIndex : cellData.vertexIndex) {
//        //                    for (int i = 0; i < 3; i++) {
//        //                        trianglesVector.push_back(edgeVertices[vertexIndex][i]);
//        //                        normalsVector.push_back(edgeNormals[vertexIndex][i]);
//        //                    }
//        tempScene->volumeData[voxelIndex].vertices[indexCount] = edgeVertices[vertexIndex];
//        tempScene->volumeData[voxelIndex].normals[indexCount] = edgeNormals[vertexIndex];
//        indexCount += 1;
//    }
//    tempScene->volumeData[voxelIndex].vertexCount = indexCount;
//}

void updateIsovalues(int x, int y, int z, float isovalue, int radius) {
    for (int xd = -radius; xd <= radius; xd++) {
        for (int yd = -radius; yd <= radius; yd++) {
            for (int zd = -radius; zd <= radius; zd++) {
                if (xd*xd + yd*yd + zd*zd > radius*radius) continue;
                int index = ((x+xd)*VOLUME_SIZE*VOLUME_SIZE + (y+yd)*VOLUME_SIZE + (z+zd));
                if (index >= 0 && index < VOLUME_SIZE*VOLUME_SIZE*VOLUME_SIZE) {
                    tempScene->volumeData[index].isovalue = isovalue;
                }
            }
        }
    }
//    for (int xd = -(radius+1); xd <= radius+1; xd++) {
//        for (int yd = -(radius+1); yd <= radius+1; yd++) {
//            for (int zd = -(radius+1); zd <= radius+1; zd++) {
//                recomputeVertices(x + xd, y + yd, z + zd);
//            }
//        }
//    }
}

bool hitOctreeNode(const Ray &ray, glm::vec3 min, glm::vec3 max) {
    if (glm::length(max - min) <= 2.0f) {
        std::cout << "Hit Node (" << ((int)min[0]) << "," << ((int)min[1]) << "," << ((int)min[2]) << ")" << std::endl;
        int x = (int)min[0];
        int y = (int)min[1];
        int z = (int)min[2];
        updateIsovalues(x, y, z, editingIsovalue, 2);
        return true;
    }
    
    Cube cube(min, max);
    
    if (cube.intersects(ray)) {
        // Use integer division to ensure it is at a voxel corner
        glm::vec3 center((max[0] + min[0]) / 2,
                         (max[1] + min[1]) / 2,
                         (max[2] + min[2]) / 2);
        
        hitOctreeNode(ray, min, center);
        hitOctreeNode(ray, glm::vec3(min[0],    min[1],    center[2]), glm::vec3(center[0], center[1], max[2]));
        hitOctreeNode(ray, glm::vec3(center[0], min[1],    min[2]),    glm::vec3(max[0],    center[1], center[2]));
        hitOctreeNode(ray, glm::vec3(center[0], min[1],    center[2]), glm::vec3(max[0],    center[1], max[2]));
        hitOctreeNode(ray, glm::vec3(min[0],    center[1], min[2]),    glm::vec3(center[0], max[1],    center[2]));
        hitOctreeNode(ray, glm::vec3(min[0],    center[1], center[2]), glm::vec3(center[0], max[1],    max[2]));
        hitOctreeNode(ray, glm::vec3(center[0], center[1], min[2]),    glm::vec3(max[0],    max[1],    center[2]));
        hitOctreeNode(ray, center, max);
        return true;
    }
    return false;
}

static void mouse_press_callback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        std::cout << "Mouse button pressed" << std::endl;
        std::cout << "Camera Pos:    (" << camera.pos.x << "," << camera.pos.y << "," << camera.pos.z << ")" << std::endl;
        std::cout << "Camera LookAt: (" << camera.lookAt.x << "," << camera.lookAt.y << "," << camera.lookAt.z << ")" << std::endl;
        Ray ray(camera.pos, camera.lookAt - camera.pos);
        
        double time1 = glfwGetTime();
        hitOctreeNode(ray, glm::vec3(0, 0, 0), glm::vec3(VOLUME_SIZE-1, VOLUME_SIZE-1, VOLUME_SIZE-1));
        double time2 = glfwGetTime();
        std::cout << "Time to traverse octree: " << (time2 - time1) << std::endl;
        
        setupVolumeTriangles(*tempScene);
        double time3 = glfwGetTime();
        std::cout << "Time to setup triangles: " << (time3 - time2) << std::endl;
        std::cout << "Total Edit Time: " << (time3 - time1) << std::endl;
    }
}

VoxelCube voxelForVolumePos(Scene &scene, int x, int y, int z) {
    glm::vec3 pos(x, y, z);
    Voxel voxel = scene.volumeData[x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z];
    return { pos, voxel.isovalue };
}

glm::vec3 gradientForPoint(Scene &scene, int x, int y, int z) {
    int minx = (x-voxelScale >= 0) ? x-voxelScale : x;
    int maxx = (x+voxelScale <= VOLUME_SIZE-1) ? x+voxelScale : x;
    float divisorx = (x == 0 || x == VOLUME_SIZE-1) ? 1.0f : 2.0f;
    float gradx = (voxelForVolumePos(scene, maxx, y, z).isovalue - voxelForVolumePos(scene, minx, y, z).isovalue) / divisorx;
    // grad[0] = (maxx - minx) / divisorx;
    
    int miny = (y-voxelScale >= 0) ? y-voxelScale : y;
    int maxy = (y+voxelScale <= VOLUME_SIZE-1) ? y+voxelScale : y;
    float divisory = (y == 0 || y == VOLUME_SIZE-1) ? 1.0f : 2.0f;
    float grady = (voxelForVolumePos(scene, x, maxy, z).isovalue - voxelForVolumePos(scene, x, miny, z).isovalue) / divisory;
    // grad[1] = (maxy - miny) / divisory;
    
    int minz = (z-voxelScale >= 0) ? z-voxelScale : z;
    int maxz = (z+voxelScale <= VOLUME_SIZE-1) ? z+voxelScale : z;
    float divisorz = (z == 0 || z == VOLUME_SIZE-1) ? 1.0f : 2.0f;
    float gradz = (voxelForVolumePos(scene, x, y, maxz).isovalue - voxelForVolumePos(scene, x, y, minz).isovalue) / divisorz;
    // grad[2] = (maxz - minz) / divisorz;
    
    return glm::vec3(gradx, grady, gradz);
}

glm::vec3 interpolateVertex(float isorange, VoxelCube v1, VoxelCube v2) {
    glm::vec3 p1 = v1.pos;
    glm::vec3 p2 = v2.pos;
    float iso1 = v1.isovalue;
    float iso2 = v2.isovalue;
    if (p2[0] < p1[0] || p2[1] < p1[1] || p2[2] < p1[2]) {
        glm::vec3 temp = p1;
        p1 = p2;
        p2 = temp;
        float temp2 = iso1;
        iso1 = iso2;
        iso2 = temp2;
    }
    
    if (abs(iso1 - iso2) > 0.00001) {
        return p1 + (p2 - p1) / (iso2 - iso1) * (isorange - iso1);
    } else {
        return p1;
    }
}

glm::vec3 interpolateNormal(Scene &scene, int x, int y, int z, int isorange, int corner1, int corner2) {
    // Get gradient for each corner
    glm::vec3 grad[2];
    VoxelCube voxels[2];
    for (int i = 0; i < 2; i++) {
        int corner = (i == 0) ? corner1 : corner2;
        int cornerOffset[3] = {0};
        cornerOffset[0] = cornerOffsets[corner][0] * voxelScale;
        cornerOffset[1] = cornerOffsets[corner][1] * voxelScale;
        cornerOffset[2] = cornerOffsets[corner][2] * voxelScale;
        grad[i] = gradientForPoint(scene, x + cornerOffsets[corner][0] * voxelScale, y + cornerOffsets[corner][1] * voxelScale, z + cornerOffsets[corner][2] * voxelScale);
        voxels[i] = voxelForVolumePos(scene, x + cornerOffsets[corner][0] * voxelScale, y + cornerOffsets[corner][1] * voxelScale, z + cornerOffsets[corner][2] * voxelScale);
    }
    
    glm::vec3 grad1 = grad[0];
    glm::vec3 grad2 = grad[1];
    float iso1 = voxels[0].isovalue;
    float iso2 = voxels[1].isovalue;
    if (iso2 < iso1) {
        float temp = iso1;
        iso1 = iso2;
        iso2 = temp;
        glm::vec3 temp2 = grad1;
        grad1 = grad2;
        grad2 = temp2;
    }
    
    if (abs(iso1 - iso2) > 0.00001) {
        return grad1 + (grad2 - grad1) / (iso2 - iso1) * (isorange - iso1);
    } else {
        return grad1;
    }
}

void setupVolumeIsovalues(Scene &scene) {
    OSN::Noise<3> noise(time(0));
    for (int x = 0; x < VOLUME_SIZE; x++) {
        for (int y = 0; y < VOLUME_SIZE; y++) {
            for (int z = 0; z < VOLUME_SIZE; z++) {
                int index = (x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z);
                scene.volumeData[index].isovalue = getIsovalue(noise, x, y, z, shapeIndex);
//                for (int i = 0; i < 15; i++) {
//                    scene.volumeData[index].vertices[i] = glm::vec3(0);
//                    scene.volumeData[index].normals[i] = glm::vec3(0);
//                }
            }
        }
    }
//    for (int x = 1; x < VOLUME_SIZE; x++) {
//        for (int y = 1; y < VOLUME_SIZE; y++) {
//            for (int z = 1; z < VOLUME_SIZE; z++) {
//                recomputeVertices(x, y, z);
//            }
//        }
//    }
}

void computeTrianglesForVoxel(Scene &scene, int x, int y, int z) {
    glm::vec3 cornerPositions[8] = {
        glm::vec3(x-voxelScale, y-voxelScale, z-voxelScale),
        glm::vec3(x,            y-voxelScale, z-voxelScale),
        glm::vec3(x-voxelScale, y,            z-voxelScale),
        glm::vec3(x,            y,            z-voxelScale),
        glm::vec3(x-voxelScale, y-voxelScale, z),
        glm::vec3(x,            y-voxelScale, z),
        glm::vec3(x-voxelScale, y,            z),
        glm::vec3(x,            y,            z)
    };
    
    int cornerIndices[8] = {
        (x-voxelScale)*VOLUME_SIZE*VOLUME_SIZE + (y-voxelScale)*VOLUME_SIZE + (z-voxelScale),
        (x           )*VOLUME_SIZE*VOLUME_SIZE + (y-voxelScale)*VOLUME_SIZE + (z-voxelScale),
        (x-voxelScale)*VOLUME_SIZE*VOLUME_SIZE + (y           )*VOLUME_SIZE + (z-voxelScale),
        (x           )*VOLUME_SIZE*VOLUME_SIZE + (y           )*VOLUME_SIZE + (z-voxelScale),
        (x-voxelScale)*VOLUME_SIZE*VOLUME_SIZE + (y-voxelScale)*VOLUME_SIZE + (z),
        (x           )*VOLUME_SIZE*VOLUME_SIZE + (y-voxelScale)*VOLUME_SIZE + (z),
        (x-voxelScale)*VOLUME_SIZE*VOLUME_SIZE + (y           )*VOLUME_SIZE + (z),
        (x           )*VOLUME_SIZE*VOLUME_SIZE + (y           )*VOLUME_SIZE + (z)
    };
    
    unsigned long cubeIndex = 0;
    float isorange = 0;
    if ((scene.volumeData[cornerIndices[0]].isovalue) <= isorange) cubeIndex |= 1;
    if ((scene.volumeData[cornerIndices[1]].isovalue) <= isorange) cubeIndex |= 2;
    if ((scene.volumeData[cornerIndices[2]].isovalue) <= isorange) cubeIndex |= 4;
    if ((scene.volumeData[cornerIndices[3]].isovalue) <= isorange) cubeIndex |= 8;
    if ((scene.volumeData[cornerIndices[4]].isovalue) <= isorange) cubeIndex |= 16;
    if ((scene.volumeData[cornerIndices[5]].isovalue) <= isorange) cubeIndex |= 32;
    if ((scene.volumeData[cornerIndices[6]].isovalue) <= isorange) cubeIndex |= 64;
    if ((scene.volumeData[cornerIndices[7]].isovalue) <= isorange) cubeIndex |= 128;
    
    /* Cube is entirely in/out of the surface */
    if (edgeTable[cubeIndex] == 0) return;
    
    /* Transvoxel */
    auto cellClass = regularCellClass[cubeIndex];
    auto cellData = regularCellData[cellClass];
    unsigned short vertexData[12] = {0};
    for (int i = 0; i < 12; i++) {
        vertexData[i] = regularVertexData[cubeIndex][i];
    }
    
    /* Find the vertices where the surface intersects the cube for each edge */
    std::vector<glm::vec3> edgeVertices;
    std::vector<glm::vec3> edgeNormals;
    for (auto vertex : vertexData) {
        if (!vertex) break;
        unsigned char corner1 = (vertex >> 4) & 0x000F;
        unsigned char corner2 = vertex & 0x000F;
        float isovalue1 = scene.volumeData[cornerIndices[corner1]].isovalue;
        float isovalue2 = scene.volumeData[cornerIndices[corner2]].isovalue;
        float t = isovalue2 / (isovalue2 - isovalue1);
        glm::vec3 vertexPos = cornerPositions[corner1] * t + cornerPositions[corner2] * (1-t);
        edgeVertices.push_back(vertexPos);
        edgeNormals.push_back(interpolateNormal(scene, x, y, z, isorange, corner1, corner2));
    }
    
    /* Create the triangles */
    for (auto vertexIndex : cellData.vertexIndex) {
        for (int i = 0; i < 3; i++) {
            trianglesVector.push_back(edgeVertices[vertexIndex][i]);
            normalsVector.push_back(edgeNormals[vertexIndex][i]);
        }
    }
}

void setupVolumeTriangles(Scene &scene) {
    double time1 = glfwGetTime();
    trianglesVector.clear();
    normalsVector.clear();
    for (int x = voxelScale; x < VOLUME_SIZE; x += voxelScale) {
        for (int y = voxelScale; y < VOLUME_SIZE; y += voxelScale) {
            for (int z = voxelScale; z < VOLUME_SIZE / 2 - voxelScale; z += voxelScale) {
                computeTrianglesForVoxel(scene, x, y, z);
            }
        }
    }
    voxelScale *= 2;
    for (int x = voxelScale; x < VOLUME_SIZE; x += voxelScale) {
        for (int y = voxelScale; y < VOLUME_SIZE; y += voxelScale) {
            for (int z = VOLUME_SIZE / 2; z < VOLUME_SIZE; z += voxelScale) {
                computeTrianglesForVoxel(scene, x, y, z);
            }
        }
    }
    voxelScale /= 2;
    
//    double time1 = glfwGetTime();
//    trianglesVector.clear();
//    normalsVector.clear();
//    for (int i = 0; i < scene.volumeData.size(); i++) {
//        for (int j = 0; j < scene.volumeData[i].vertexCount; j++) {
//            trianglesVector.push_back(scene.volumeData[i].vertices[j][0]);
//            trianglesVector.push_back(scene.volumeData[i].vertices[j][1]);
//            trianglesVector.push_back(scene.volumeData[i].vertices[j][2]);
//            normalsVector.push_back(scene.volumeData[i].normals[j][0]);
//            normalsVector.push_back(scene.volumeData[i].normals[j][1]);
//            normalsVector.push_back(scene.volumeData[i].normals[j][2]);
//        }
//    }
    double time2 = glfwGetTime();
    std::cout << "Time to re-add all vertices to the vertex and normal buffers: " << (time2 - time1) << std::endl;
    
    // Initialize the vertex array object
    glGenVertexArrays(1, &volumeVAO);
    glBindVertexArray(volumeVAO);
    
    std::cout << "Triangle Vector Size: " << trianglesVector.size() / 3 << std::endl;
//    std::cout << "Triangle Vector Size: " << scene.volumeData.size() * 30 << std::endl;
    //generate vertex buffer to hand off to OGL
    glGenBuffers(1, &volumeVBO);
    glBindBuffer(GL_ARRAY_BUFFER, volumeVBO);
    glBufferData(GL_ARRAY_BUFFER, trianglesVector.size() * sizeof(float), trianglesVector.data(), GL_DYNAMIC_DRAW);
//    glBufferData(GL_ARRAY_BUFFER, scene.volumeData.size() * sizeof(Voxel), scene.volumeData.data(), GL_DYNAMIC_DRAW);
//    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Voxel), (GLvoid*)offsetof(Voxel, vertices));
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    std::cout << "Normals Vector Size: " << normalsVector.size() / 3 << std::endl;
//    std::cout << "Normals Vector Size: " << scene.volumeData.size() * 30 << std::endl;
    glGenBuffers(1, &volumeVBO);
    glBindBuffer(GL_ARRAY_BUFFER, volumeVBO);
    glBufferData(GL_ARRAY_BUFFER, normalsVector.size() * sizeof(float), normalsVector.data(), GL_DYNAMIC_DRAW);
//    glBufferData(GL_ARRAY_BUFFER, scene.volumeData.size() * sizeof(Voxel), scene.volumeData.data(), GL_DYNAMIC_DRAW);
//    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Voxel), (GLvoid*)offsetof(Voxel, normals));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (GLvoid*)0);
    glEnableVertexAttribArray(1);
}

void setupVolumeData(Scene &scene) {
    setupVolumeIsovalues(scene);
    setupVolumeTriangles(scene);
}

void resetVolumeData(Scene &scene) {
    trianglesVector.clear();
    normalsVector.clear();
    setupVolumeData(scene);
}

void drawVolumeData(float aspect) {
    
    volumeShader->use();
    
//    auto P = std::make_shared<MatrixStack>();
    auto P = std::make_shared<MatrixStack>();
    P->perspective(45.0f, aspect, 0.01f, 500.0f);
    glUniformMatrix4fv(glGetUniformLocation(volumeShader->program, "P"), 1, GL_FALSE, value_ptr(P->topMatrix()));
    // Send model matrix to GPU
    auto M = std::make_shared<MatrixStack>();
    // M->rotate(0, Vector3f(0, 1, 0));
    // M->scale(0.0625);
    M->rotate(t, glm::vec3(0, 1, 0));
    // M->scale(voxelScale/2.0f);
    M->translate(glm::vec3(-VOLUME_SIZE/2.0f, -VOLUME_SIZE/2.0f, -VOLUME_SIZE/2.0f));
    glUniformMatrix4fv(glGetUniformLocation(volumeShader->program, "M"), 1, GL_FALSE, value_ptr(M->topMatrix()));
    // Send view matrix to GPU
    auto V = std::make_shared<MatrixStack>();
    V->lookAt(camera.pos, camera.lookAt, camera.up);
    glUniformMatrix4fv(glGetUniformLocation(volumeShader->program, "V"), 1, GL_FALSE, value_ptr(V->topMatrix()));
    
//    std::cout << "Triangle Vector Size: " << trianglesVector.size() << std::endl;
    glBindVertexArray(volumeVAO);
    glDrawArrays(GL_TRIANGLES, 0, trianglesVector.size());
    glBindVertexArray(0);
    
    if (isRotating) {
        t += 0.01;
    }
}

static void updateCamera() {
    glm::vec3 gaze = camera.lookAt - camera.pos;
    glm::vec3 w = glm::normalize(-gaze);
    glm::vec3 u = glm::normalize(glm::cross(camera.up, w));
    if (abs(camera.velW) > CAMERA_STOPPED_THRESHOLD) {
        camera.pos += camera.velW * w;
        camera.lookAt += camera.velW * w;
    }
    if (abs(camera.velU) > CAMERA_STOPPED_THRESHOLD) {
        camera.pos += camera.velU * u;
        camera.lookAt += camera.velU * u;
    }
}

int main() {

    // Load GLFW and create a window
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    auto window = glfwCreateWindow(mWidth, mHeight, "OpenGL", nullptr, nullptr);

    // Check for valid context
    if (!window) {
        fprintf(stderr, "Failed to Create OpenGL Context");
        return EXIT_FAILURE;
    }

    // Create context and load OpenGL functions
    glfwMakeContextCurrent(window);
    gladLoadGL();
    fprintf(stderr, "OpenGL %s\n", glGetString(GL_VERSION));

    // Disable cursor (allows unlimited scrolling)
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    
    // Set vsync.
    glfwSwapInterval(1);
    // Register callbacks for events
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_move_callback);
    glfwSetMouseButtonCallback(window, mouse_press_callback);
    
    // Build and compile our shader program
    Shader simpleShader("../../../LOGLEngine/Shaders/basic.vert", "../../../LOGLEngine/Shaders/basic.frag");
    volumeShader = std::make_shared<Shader>("../../../LOGLEngine/Shaders/volume.vert", "../../../LOGLEngine/Shaders/volume.frag");
    
    Scene scene(VOLUME_SIZE);
    tempScene = &scene;
    setupVolumeData(scene);
    
    camera.pos = glm::vec3(0, 0, -VOLUME_SIZE*1.5f);
    camera.lookAt = glm::vec3(0, 0, 0);
    camera.up = glm::vec3(0, 1, 0);
    
    // Enable z-buffer test.
    glEnable(GL_DEPTH_TEST);
    // Background fill color
    glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
    
    // Rendering loop
    while (!glfwWindowShouldClose(window)) {
        // Check and call events
        glfwPollEvents();

        // Clear color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Set viewport size for rendering
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);

        drawVolumeData(width/(float)height);
        
        // Update camera;
        updateCamera();
        
        // Swap the buffers
        glfwSwapBuffers(window);
    }
    
    // Destroy the window
    glfwTerminate();
    return EXIT_SUCCESS;
}
