//
//  Noise.h
//  Shepherd
//
//  Created by Tyler Dahl on 2/3/17.
//
//

#ifndef Noise_h
#define Noise_h

#include <vector>
#include <time.h>

#include "OpenSimplexNoise.hpp"
#include "glm/glm.hpp"

#include "NoiseProperties.h"

class Noise {
public:

/// Generate fractal noise using the Diamond-Square algorithm.
/// For best results, use sizes in the form of (2^pow)+1.
static std::vector<std::vector<float>> GenerateDiamondSquare(int size) {
    int maxStepSize = 256;
    std::vector<std::vector<float>> map;

    // Initialize map to all 0s
    for (int row = 0; row < size; row++) {
        std::vector<float> rowVector;
        for (int col = 0; col < size; col++) {
            rowVector.push_back(0);
        }
        map.push_back(rowVector);
    }
    
    // Randomize corners
    for (int row = 0; row < size; row += maxStepSize) { // 256
        for (int col = 0; col < size; col += maxStepSize) { // 256
            map[row][col] = (float)(rand() % maxStepSize);
        }
    }
    
    // Perform diamond-square
    //float max = 0;
    for (int stepSize = maxStepSize; stepSize > 1; stepSize /= 2) {
        int halfStepSize = stepSize/2;
        
        // Squares
        for (int row = stepSize; row < size; row += stepSize) { // 256
            for (int col = stepSize; col < size; col += stepSize) { // 256
                float h1 = map[row-stepSize][col-stepSize]; // 0, 0
                float h2 = map[row-stepSize][col]; // 0, 256
                float h3 = map[row][col-stepSize]; // 256, 0
                float h4 = map[row][col]; // 256, 256
                map[row-stepSize/2][col-stepSize/2] = (h1+h2+h3+h4)/4 + (rand() % halfStepSize)-halfStepSize/2;
            }
        }
        
        // Diamonds - top and left of every cell
        for (int row = stepSize; row < size; row += stepSize) { // 256
            for (int col = stepSize; col < size; col += stepSize) { // 256
                float h1 = map[row-stepSize][col-stepSize]; // 0, 0
                float h2 = map[row-stepSize][col]; // 0, 256
                float h3 = map[row][col-stepSize]; // 256, 0
//                float h4 = map[row][col]; // 256, 256
                float middle = map[row-stepSize/2][col-stepSize/2];
                
                float topNeighbor = 0;
                int topCount = 3;
                if (col-3*stepSize/2 >= 0) {
                    topNeighbor = map[row-stepSize/2][col-3*stepSize/2];
                    topCount++;
                }
                float leftNeighbor = 0;
                int leftCount = 3;
                if (row-3*stepSize/2 >= 0) {
                    leftNeighbor = map[row-3*stepSize/2][col-stepSize/2];
                    leftCount++;
                }
                
                // top
                map[row-stepSize/2][col-stepSize] = (h1+middle+h3+topNeighbor)/topCount + (rand() % halfStepSize)-halfStepSize/2;
                // left
                map[row-stepSize][col-stepSize/2] = (h1+middle+h2+leftNeighbor)/leftCount + (rand() % halfStepSize)-halfStepSize/2;
            }
        }
        
        // Diamonds - bottom and right edges of map
        int lastRow = size-1;
        for (int col = stepSize; col < size; col += stepSize) {
            float h3 = map[lastRow][col-stepSize]; // 256, 0
            float h4 = map[lastRow][col]; // 256, 256
            float middle = map[lastRow-stepSize/2][col-stepSize/2];
            map[lastRow][col-stepSize/2] = (h3+middle+h4)/3 + (rand() % halfStepSize)-halfStepSize/2;
        }
        int lastCol = size-1;
        for (int row = stepSize; row < size; row += stepSize) {
            float h2 = map[row-stepSize][lastCol]; // 0, 256
            float h4 = map[row][lastCol]; // 256, 256
            float middle = map[row-stepSize/2][lastCol-stepSize/2];
            map[row-stepSize/2][lastCol] = (h2+middle+h4)/3 + (rand() % halfStepSize)-halfStepSize/2;
        }
    }
    
    // Re-scale
    float min = FLT_MAX;
    float max = FLT_MIN;
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (map[row][col] < min) min = map[row][col];
            if (map[row][col] > max) max = map[row][col];
        }
    }
    float range = max - min;
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            map[row][col] = (map[row][col] - min) / range * 255;
        }
    }
    
    return map;
}

//glm::mat2 Noise::generatePerlin(int size) {
//    
//}
    
    static std::vector<std::vector<float>> GenerateSimplex(NoiseProperties &properties, int size) {
        OSN::Noise<2> noise(properties.seed);
        std::vector<std::vector<float>> map;
        
        // Initialize map to all 0s
        for (int row = 0; row < size; row++) {
            std::vector<float> rowVector;
            for (int col = 0; col < size; col++) {
                rowVector.push_back(0);
            }
            map.push_back(rowVector);
        }
        
        float range = 256.0f;
        float halfScale = (size / range / 2);
        
        // Set map values
        float frequency = properties.frequency;
        float initialHeight = properties.octaveHeight / properties.frequency;
        float octaves = properties.octaves;
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                float nx = ((float)x/range - halfScale) * frequency;
                float ny = ((float)y/range - halfScale) * frequency;
                
                float scale = pow(2, 1);
                map[y][x] += initialHeight/scale * (0.4+noise.eval(scale * nx, scale * ny));
                
                // Increase octaves for more detailed terrain
                for (int oct = 1; oct < octaves; oct++) {
                    float scale = pow(2, oct);
                    map[y][x] += map[y][x]/scale * (0.4+noise.eval(scale * nx, scale * ny));
                }
                //            map[y][x] = pow(map[y][x], 1.2);
                map[y][x] = (map[y][x] >= 0) ? pow(map[y][x], 1.3) : 0;
                //            map[y][x] = -std::fabsf(map[y][x]);
                
                scale = pow(2, 5);
                map[y][x] += map[y][x]/scale * (0.4+noise.eval(scale * nx*5, scale * ny*5));
                
                // Increase height with distance from center
//                float dist = sqrt(pow(y-size/2, 2) + pow(x-size/2, 2));
//                map[y][x] = (map[y][x]+dist/10.0f);// * (1 + dist/256.0f);
//                map[y][x] += -fabs(50.0f * noise.eval(nx / 2.0f, ny / 2.0f));
            }
        }
        
        return map;
    }

static std::vector<std::vector<float>> GenerateSimplex(int size) {
    OSN::Noise<2> noise(time(0));
    std::vector<std::vector<float>> map;
    
    // Initialize map to all 0s
    for (int row = 0; row < size; row++) {
        std::vector<float> rowVector;
        for (int col = 0; col < size; col++) {
            rowVector.push_back(0);
        }
        map.push_back(rowVector);
    }
    
    float range = 256.0f;
    float halfScale = (size / range / 2);
    
    // Set map values
    float frequency = 5.0f;
    float initialHeight = 40.0f / frequency;
    float octaves = 5;
    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            float nx = ((float)x/range - halfScale) * frequency;
            float ny = ((float)y/range - halfScale) * frequency;
            
            float scale = pow(2, 1);
            map[y][x] += initialHeight/scale * (0.4+noise.eval(scale * nx, scale * ny));
            
            // Increase octaves for more detailed terrain
            for (int oct = 1; oct < octaves; oct++) {
                float scale = pow(2, oct);
                map[y][x] += map[y][x]/scale * (0.4+noise.eval(scale * nx, scale * ny));
            }
//            map[y][x] = pow(map[y][x], 1.2);
            map[y][x] = (map[y][x] >= 0) ? pow(map[y][x], 1.3) : 0;
//            map[y][x] = -std::fabsf(map[y][x]);
            
            scale = pow(2, 5);
            map[y][x] += map[y][x]/scale * (0.4+noise.eval(scale * nx*5, scale * ny*5));
        }
    }
    
    return map;
}
    
    static std::vector<std::vector<float>> SmoothTerrain(std::vector<std::vector<float>> map, int size, int iterations, int kernelSize) {
        OSN::Noise<2> noise(time(0));
        std::vector<std::vector<float>> newMap;
        
        // Initialize map to all 0s
        for (int row = 0; row < size; row++) {
            std::vector<float> rowVector;
            for (int col = 0; col < size; col++) {
                rowVector.push_back(0);
            }
            newMap.push_back(rowVector);
        }
        
        // Set map values
        for (int steps = 0; steps < iterations; steps++) {
            for (int y = 0; y < size; y++) {
                for (int x = 0; x < size; x++) {
                    newMap[y][x] = 0;
                    int count = 0;
                    for (int ky = -kernelSize/2; ky < kernelSize/2; ky++) {
                        for (int kx = -kernelSize/2; kx < kernelSize/2; kx++) {
                            if (y-ky >= 0 && y-ky < size && x-kx >= 0 && x-kx < size) {
                                count++;
                                newMap[y][x] += map[y-ky][x-kx];
                            }
                        }
                    }
                    newMap[y][x] /= count;
                }
            }
            map = newMap;
        }
        
        return newMap;
    }
};

#endif /* Noise_h */
