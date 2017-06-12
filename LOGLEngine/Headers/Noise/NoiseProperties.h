//
//  NoiseProperties.h
//  Shepherd
//
//  Created by Tyler Dahl on 3/4/17.
//
//

#ifndef NoiseProperties_h
#define NoiseProperties_h

#include <time.h>

struct NoiseProperties {
    time_t seed = 0;
    float frequency = 5.0f;
    int octaves = 5;
    float octaveHeight = 40.0f;
};

#endif /* NoiseProperties_h */
