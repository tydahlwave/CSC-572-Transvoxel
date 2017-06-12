// Preprocessor Directives
#ifndef LOGL_ENGINE
#define LOGL_ENGINE
#pragma once

// System Headers
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/scene.h>
#include <btBulletDynamicsCommon.h>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// Reference: https://github.com/nothings/stb/blob/master/stb_image.h#L4
// To use stb_image, add this in *one* C++ source file.
//     #define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

// Define Some Constants
#define LEFT 0
#define RIGHT 1
#define CAMERA_SPEED 0.2
#define CAMERA_STOPPED_THRESHOLD 0.1
#define NUM_MODELS 5
#define VOLUME_SIZE 2
#define POINT_SIZE 50
#define SPHERE_RADIUS 0.5

const int mWidth = 1080;
const int mHeight = 920;

#endif //~ LOGLEngine Header
