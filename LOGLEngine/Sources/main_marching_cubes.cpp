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
float voxelScale = 1.0f;
float t = 0;
int shapeIndex = 0;

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
Scene *tempScene;

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        } else if (key == GLFW_KEY_SPACE) {
            shapeIndex += 1;
            resetVolumeData(*tempScene);
        } else if (key == GLFW_KEY_W) {
            camera.velW = -CAMERA_SPEED;
        } else if (key == GLFW_KEY_S) {
            camera.velW = CAMERA_SPEED;
        } else if (key == GLFW_KEY_A) {
            camera.velU = -CAMERA_SPEED;
        } else if (key == GLFW_KEY_D) {
            camera.velU = CAMERA_SPEED;
        }
    } else if (action == GLFW_RELEASE) {
        if (key == GLFW_KEY_W) {
            if (camera.velW < -CAMERA_STOPPED_THRESHOLD) {
                camera.velW = 0;
            } else {
                camera.velW = CAMERA_SPEED;
            }
        } else if (key == GLFW_KEY_S) {
            if (camera.velW > CAMERA_STOPPED_THRESHOLD) {
                camera.velW = 0;
            } else {
                camera.velW = -CAMERA_SPEED;
            }
        } else if (key == GLFW_KEY_A) {
            if (camera.velU < -CAMERA_STOPPED_THRESHOLD) {
                camera.velU = 0;
            } else {
                camera.velU = CAMERA_SPEED;
            }
        } else if (key == GLFW_KEY_D) {
            if (camera.velU > CAMERA_STOPPED_THRESHOLD) {
                camera.velU = 0;
            } else {
                camera.velU = -CAMERA_SPEED;
            }
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

void triangleSetup() {
    // Set up vertex data (and buffer(s)) and attribute pointers
    GLfloat vertices[] = {
        // Positions         // Colors
        0.5f, -0.5f, 0.0f,   1.0f, 0.0f, 0.0f,  // Bottom Right
        -0.5f, -0.5f, 0.0f,   0.0f, 1.0f, 0.0f,  // Bottom Left
        0.0f,  0.5f, 0.0f,   0.0f, 0.0f, 1.0f   // Top
    };
    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    // Color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
    
    glBindVertexArray(0); // Unbind VAO
}

VoxelCube voxelForVolumePos(Scene &scene, int x, int y, int z) {
    glm::vec3 pos(x, y, z);
    Voxel voxel = scene.volumeData[x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z];
    return { pos, voxel.isovalue };
}

// float vertexForEdge(int edge, int x, int y, int z) {
//    int index = x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z;
//    int corners[2] = edgeVertices[edge];
//    return
// }

glm::vec3 gradientForPoint(Scene &scene, int x, int y, int z) {
    int minx = (x-1 >= 0) ? x-1 : x;
    int maxx = (x+1 <= VOLUME_SIZE-1) ? x+1 : x;
    float divisorx = (x == 0 || x == VOLUME_SIZE-1) ? 1.0f : 2.0f;
    float gradx = (voxelForVolumePos(scene, maxx, y, z).isovalue - voxelForVolumePos(scene, minx, y, z).isovalue) / divisorx;
    // grad[0] = (maxx - minx) / divisorx;
    
    int miny = (y-1 >= 0) ? y-1 : y;
    int maxy = (y+1 <= VOLUME_SIZE-1) ? y+1 : y;
    float divisory = (y == 0 || y == VOLUME_SIZE-1) ? 1.0f : 2.0f;
    float grady = (voxelForVolumePos(scene, x, maxy, z).isovalue - voxelForVolumePos(scene, x, miny, z).isovalue) / divisory;
    // grad[1] = (maxy - miny) / divisory;
    
    int minz = (z-1 >= 0) ? z-1 : z;
    int maxz = (z+1 <= VOLUME_SIZE-1) ? z+1 : z;
    float divisorz = (z == 0 || z == VOLUME_SIZE-1) ? 1.0f : 2.0f;
    float gradz = (voxelForVolumePos(scene, x, y, maxz).isovalue - voxelForVolumePos(scene, x, y, minz).isovalue) / divisorz;
    // grad[2] = (maxz - minz) / divisorz;
    
    return glm::vec3(gradx, grady, gradz);
}

//void interpolateGradient(Scene &scene, int x, int y, int z, int isorange, int corners[2], float finalGrad[3]) {
//    // Get gradient for each corner
//    float grad[2][3];
//    VoxelCube voxels[2];
//    for (int i = 0; i < 2; i++) {
//        int corner = corners[i];
//        int cornerOffset[3] = {0};
//        cornerOffset[0] = cornerOffsets[corner][0];
//        cornerOffset[1] = cornerOffsets[corner][1];
//        cornerOffset[2] = cornerOffsets[corner][2];
//        gradientForPoint(scene, x + cornerOffsets[corner][0], y + cornerOffsets[corner][1], z + cornerOffsets[corner][2], grad[i]);
//        voxels[i] = voxelForVolumePos(scene, x + cornerOffsets[corner][0], y + cornerOffsets[corner][1], z + cornerOffsets[corner][2]);
//    }
//
//    float *grad1 = grad[0];
//    float *grad2 = grad[1];
//    float iso1 = voxels[0].isovalue;
//    float iso2 = voxels[1].isovalue;
//    if (iso2 < iso1) {
//        float temp = iso1;
//        iso1 = iso2;
//        iso2 = temp;
//        float *temp2 = grad1;
//        grad1 = grad2;
//        grad2 = temp2;
//    }
//
//    if (abs(iso1 - iso2) > 0.00001) {
//        finalGrad[0] = grad1[0] + (grad2[0] - grad1[0]) / (iso2 - iso1) * (isorange - iso1);
//        finalGrad[1] = grad1[1] + (grad2[1] - grad1[1]) / (iso2 - iso1) * (isorange - iso1);
//        finalGrad[2] = grad1[2] + (grad2[2] - grad1[2]) / (iso2 - iso1) * (isorange - iso1);
//    } else {
//        finalGrad[0] = grad1[0];
//        finalGrad[1] = grad1[1];
//        finalGrad[2] = grad1[2];
//    }
//
//    // finalGrad[0] = (grad[0][0] + grad[1][0]) / 2.0f;
//    // finalGrad[1] = (grad[0][1] + grad[1][1]) / 2.0f;
//    // finalGrad[2] = (grad[0][2] + grad[1][2]) / 2.0f;
//}
//
//void interpolatePosition(float isorange, VoxelCube v1, VoxelCube v2, float finalPos[3]) {
//    glm::vec3 p1 = v1.pos;
//    glm::vec3 p2 = v2.pos;
//    float iso1 = v1.isovalue;
//    float iso2 = v2.isovalue;
//    if (p2[0] < p1[0] || p2[1] < p1[1] || p2[2] < p1[2]) {
//        glm::vec3 temp = p1;
//        p1 = p2;
//        p2 = temp;
//        float temp2 = iso1;
//        iso1 = iso2;
//        iso2 = temp2;
//    }
//
//    if (abs(iso1 - iso2) > 0.00001) {
//        finalPos[0] = p1[0] + (p2[0] - p1[0]) / (iso2 - iso1) * (isorange - iso1);
//        finalPos[1] = p1[1] + (p2[1] - p1[1]) / (iso2 - iso1) * (isorange - iso1);
//        finalPos[2] = p1[2] + (p2[2] - p1[2]) / (iso2 - iso1) * (isorange - iso1);
//    } else {
//        finalPos[0] = p1[0];
//        finalPos[1] = p1[1];
//        finalPos[2] = p1[2];
//    }
//}

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
        cornerOffset[0] = cornerOffsets[corner][0];
        cornerOffset[1] = cornerOffsets[corner][1];
        cornerOffset[2] = cornerOffsets[corner][2];
        grad[i] = gradientForPoint(scene, x + cornerOffsets[corner][0], y + cornerOffsets[corner][1], z + cornerOffsets[corner][2]);
        voxels[i] = voxelForVolumePos(scene, x + cornerOffsets[corner][0], y + cornerOffsets[corner][1], z + cornerOffsets[corner][2]);
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

void setupVolumeData(Scene &scene) {
    OSN::Noise<3> noise(time(0));
    for (int x = 0; x < VOLUME_SIZE; x++) {
        for (int y = 0; y < VOLUME_SIZE; y++) {
            for (int z = 0; z < VOLUME_SIZE; z++) {
                int index = (x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z);
                //                scene.volumeData[index].pos[0] = x*voxelScale;
                //                scene.volumeData[index].pos[1] = y*voxelScale;
                //                scene.volumeData[index].pos[2] = z*voxelScale;
                scene.volumeData[index].isovalue = getIsovalue(noise, x, y, z, shapeIndex);
            }
        }
    }
    
    for (int x = 1; x < VOLUME_SIZE; x++) {
        for (int y = 1; y < VOLUME_SIZE; y++) {
            for (int z = 1; z < VOLUME_SIZE; z++) {
                // int index = (x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z)*5*3;
                
                glm::vec3 cornerPositions[8] = {
                    glm::vec3(x-1, y-1, z-1),
                    glm::vec3(x,   y-1, z-1),
                    glm::vec3(x,   y-1, z),
                    glm::vec3(x-1, y-1, z),
                    glm::vec3(x-1, y,   z-1),
                    glm::vec3(x,   y,   z-1),
                    glm::vec3(x,   y,   z),
                    glm::vec3(x-1, y,   z)
                };
                
                int cornerIndices[8] = {
                    (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z-1),
                    (x  )*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z-1),
                    (x  )*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z  ),
                    (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y-1)*VOLUME_SIZE + (z  ),
                    (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z-1),
                    (x  )*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z-1),
                    (x  )*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z  ),
                    (x-1)*VOLUME_SIZE*VOLUME_SIZE + (y  )*VOLUME_SIZE + (z  )
                };
                
                int cubeIndex = 0;
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
                if (edgeTable[cubeIndex] == 0)
                    continue;
                
                /* Find the vertices where the surface intersects the cube for each edge */
                glm::vec3 vertlist[12];
                glm::vec3 normlist[12];
                for (int i = 0; i < 12; i++) {
                    if (edgeTable[cubeIndex] & (1 << i)) {
                        VoxelCube voxel1 = { cornerPositions[edgeVertices[i][0]], scene.volumeData[cornerIndices[edgeVertices[i][0]]].isovalue };
                        VoxelCube voxel2 = { cornerPositions[edgeVertices[i][1]], scene.volumeData[cornerIndices[edgeVertices[i][1]]].isovalue };
                        vertlist[i] = interpolateVertex(isorange, voxel1, voxel2);
                        normlist[i] = interpolateNormal(scene, x, y, z, isorange, edgeVertices[i][0], edgeVertices[i][1]);
                    }
                }
                
                /* Create the triangles */
                for (int i = 0; triTable[cubeIndex][i] != -1; i += 3) {
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {
                            trianglesVector.push_back(vertlist[triTable[cubeIndex][i+j]][k]);
                            normalsVector.push_back(normlist[triTable[cubeIndex][i+j]][k]);
                        }
                    }
                }
            }
        }
    }
    
    // Initialize the vertex array object
    glGenVertexArrays(1, &volumeVAO);
    glBindVertexArray(volumeVAO);
    std::cout << "Triangle Vector Size: " << trianglesVector.size() << std::endl;
    //generate vertex buffer to hand off to OGL
    glGenBuffers(1, &volumeVBO);
    glBindBuffer(GL_ARRAY_BUFFER, volumeVBO);
    glBufferData(GL_ARRAY_BUFFER, trianglesVector.size() * sizeof(float), trianglesVector.data(), GL_DYNAMIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    std::cout << "Normals Vector Size: " << normalsVector.size() << std::endl;
    glGenBuffers(1, &volumeVBO);
    glBindBuffer(GL_ARRAY_BUFFER, volumeVBO);
    glBufferData(GL_ARRAY_BUFFER, normalsVector.size() * sizeof(float), normalsVector.data(), GL_DYNAMIC_DRAW);
    
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(1);
    
    // for (int i = 0; i < VOLUME_SIZE*VOLUME_SIZE*VOLUME_SIZE*5*3; i++) {
    //    if (triangles[i] != 0) {
    //       cout << triangles[i] << endl;
    //    }
    // }
    // Iterate through volume, getting cubes via 8 corner points
    // Create 8-bit value to look up triangle vertices
    // Draw triangles
}

void resetVolumeData(Scene &scene) {
    trianglesVector.clear();
    normalsVector.clear();
    setupVolumeData(scene);
}

void drawVolumeData(float aspect) {
    
    // Bind the points GLSL program
    // volumeProg->bind();
    // glUniform2f(pointProg->getUniform("screenCenter"), g_width, g_height);
    // glUniform1f(pointProg->getUniform("t"), t);
    // glUniform1f(pointProg->getUniform("radius_world"), RADIUS);
    // glUniform1f(pointProg->getUniform("radius_pixel"), g_width * RADIUS);
    
    // // Send perspective matrix to GPU
    // auto P = make_shared<MatrixStack>();
    // P->perspective(45.0f, aspect, 0.01f, 500.0f);
    // // glUniformMatrix4fv(volumeProg->getUniform("P"), 1, GL_FALSE, P->topMatrix().data());
    // // Send model matrix to GPU
    // auto M = make_shared<MatrixStack>();
    // // // M->rotate(0, Vector3f(0, 1, 0));
    // // // M->scale(0.0625);
    // // glUniformMatrix4fv(volumeProg->getUniform("M"), 1, GL_FALSE, M->topMatrix().data());
    // // // Send view matrix to GPU
    // // auto V = make_shared<MatrixStack>();
    // // V->lookAt(camera.pos, camera.lookAt, camera.up);
    // // glUniformMatrix4fv(volumeProg->getUniform("V"), 1, GL_FALSE, V->topMatrix().data());
    
    // // // glEnableVertexAttribArray(0);
    // // // glBindBuffer(GL_ARRAY_BUFFER, volumeVBO);
    // // // glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, (void*) 0);
    
    // // // glPointSize(POINT_SIZE);
    // // // glDrawArrays(GL_POINTS, 0, VOLUME_SIZE*VOLUME_SIZE*VOLUME_SIZE*4);
    
    // // glDisableVertexAttribArray(0);
    // // volumeProg->unbind();
    
    //    volumeProg->bind();
    volumeShader->use();
    // for (int x = 0; x < VOLUME_SIZE; x++) {
    //    for (int y = 0; y < VOLUME_SIZE; y++) {
    //       for (int z = 0; z < VOLUME_SIZE; z++) {
    //          int index = (x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z);
    //          float isovalue = volume_data[index].isovalue;
    //          if (abs(isovalue) < 10) {
    //             M->pushMatrix();
    //             M->rotate(t*100, Vector3f(0, 1, 0));
    //             M->translate(Vector3f(volume_data[index].pos[0]-VOLUME_SIZE/2.0f/16.0f, volume_data[index].pos[1]-VOLUME_SIZE/2.0f/16.0f, volume_data[index].pos[2]-VOLUME_SIZE/2.0f/16.0f));
    //             M->scale(1.0f/2/16);
    //             drawShape(cube, phongProg, P, M, volume_data[index].pos[0], volume_data[index].pos[1], volume_data[index].pos[2], 4);
    //             M->popMatrix();
    //          }
    //       }
    //    }
    // }
    
    // glUniform3f(phongProg->getUniform("lightPos"), light.pos[0] + lightPosOffset, light.pos[1], light.pos[2]);
    // glUniform3f(phongProg->getUniform("lightColor"), light.color[0], light.color[1], light.color[2]);
    // glUniform3f(phongProg->getUniform("sunDir"), sun.dir[0], sun.dir[1], sun.dir[2]);
    // glUniform3f(phongProg->getUniform("sunColor"), sun.color[0], sun.color[1], sun.color[2]);
    // SetMaterial(4, phongProg);
    
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
    M->translate(glm::vec3(-VOLUME_SIZE/2.0f*voxelScale, -VOLUME_SIZE/2.0f*voxelScale, -VOLUME_SIZE/2.0f*voxelScale));
    glUniformMatrix4fv(glGetUniformLocation(volumeShader->program, "M"), 1, GL_FALSE, value_ptr(M->topMatrix()));
    // Send view matrix to GPU
    auto V = std::make_shared<MatrixStack>();
    V->lookAt(camera.pos, camera.lookAt, camera.up);
    glUniformMatrix4fv(glGetUniformLocation(volumeShader->program, "V"), 1, GL_FALSE, value_ptr(V->topMatrix()));
    
    //    std::cout << "Triangle Vector Size: " << trianglesVector.size() << std::endl;
    glBindVertexArray(volumeVAO);
    glDrawArrays(GL_TRIANGLES, 0, trianglesVector.size());
    glBindVertexArray(0);
    
    // Draw our first triangle
    // glBindVertexArray(VAO);
    // glDrawArrays(GL_TRIANGLES, 0, 3);
    // glBindVertexArray(0);
    
    //    volumeProg->unbind();
    t += 0.01;
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
    // Set the mouse move call back
    glfwSetCursorPosCallback(window, mouse_move_callback);
    
    // Setup triangle vertex attributes and send vertices to GPU
    triangleSetup();
    // Build and compile our shader program
    Shader simpleShader("../../../LOGLEngine/Shaders/basic.vert", "../../../LOGLEngine/Shaders/basic.frag");
    volumeShader = std::make_shared<Shader>("../../../LOGLEngine/Shaders/volume.vert", "../../../LOGLEngine/Shaders/volume.frag");
    
    Scene scene(VOLUME_SIZE);
    tempScene = &scene;
    setupVolumeData(scene);
    
    //    for (int z = 0; z < VOLUME_SIZE; z++) {
    //        for (int y = 0; y < VOLUME_SIZE; y++) {
    //            for (int x = 0; x < VOLUME_SIZE; x++) {
    //                std::cout << scene.volumeData[x*VOLUME_SIZE*VOLUME_SIZE + y*VOLUME_SIZE + z].isovalue << ", ";
    //            }
    //            std::cout << std::endl;
    //        }
    //        std::cout << std::endl << std::endl;
    //    }
    
    //    volumeProg = std::make_shared<Program>();
    //    volumeProg->setVerbose(true);
    //    volumeProg->setShaderNames(RESOURCE_DIR + "pointVert.glsl", RESOURCE_DIR + "pointFrag.glsl");
    //    volumeProg->init();
    //    volumeProg->addUniform("P");
    //    volumeProg->addUniform("M");
    //    volumeProg->addUniform("V");
    //    volumeProg->addAttribute("vertPos");
    //    volumeProg->addAttribute("vertNor");
    
    
    camera.pos = glm::vec3(0, 0, -VOLUME_SIZE*voxelScale*1.5f);
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
        
        // Draw the triangle
        //        simpleShader.use();
        //        glBindVertexArray(VAO);
        //        glDrawArrays(GL_TRIANGLES, 0, 3);
        //        glBindVertexArray(0);
        
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
