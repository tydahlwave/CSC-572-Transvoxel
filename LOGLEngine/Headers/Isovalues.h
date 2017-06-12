//
//  Isovalues.h
//  LOGLEngine
//
//  Created by Tyler Dahl on 6/12/17.
//
//

#ifndef Isovalues_h
#define Isovalues_h

float getIsovalueFor3DSimplexNoise(OSN::Noise<3> noise, int x, int y, int z) {
    float range = 256.0f;
    float halfScale = (VOLUME_SIZE / range / 2);
    
    // Set map values
    float frequency = 5.0f;
    float initialHeight = 40.0f / frequency;
    float octaves = 5;
    //    for (int z = 0; z < VOLUME_SIZE; z++) {
    //        for (int y = 0; y < VOLUME_SIZE; y++) {
    //            for (int x = 0; x < VOLUME_SIZE; x++) {
    float nx = ((float)x/range - halfScale) * frequency;
    float ny = ((float)y/range - halfScale) * frequency;
    float nz = ((float)z/range - halfScale) * frequency;
    
    float scale = pow(2, 1);
    float isovalue = initialHeight/scale * (0.4+noise.eval(scale * nx, scale * ny, scale * nz));
    
    // Increase octaves for more detailed terrain
    for (int oct = 1; oct < octaves; oct++) {
        float scale = pow(2, oct);
        isovalue += isovalue/scale * (0.4+noise.eval(scale * nx, scale * ny, scale * nz));
    }
    //            map[y][x] = pow(map[y][x], 1.2);
    isovalue = (isovalue >= 0) ? pow(isovalue, 1.3) : 0;
    //            map[y][x] = -std::fabsf(map[y][x]);
    
    scale = pow(2, 5);
    isovalue += isovalue/scale * (0.4+noise.eval(scale * nx*5, scale * ny*5, scale * nz*5));
    //            }
    //        }
    //    }
    return isovalue;
}

float getIsovalueForQ1(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    return pow(dx,4)+pow(dy,4)+pow(dz,4)-4*dx*dx-4*dy*dy*dz*dz-4*dy*dy-4*dz*dz*dx*dx-4*dz*dz-4*dx*dx*dy*dy+20.7846*dx*dy*dz+1;
}

float getIsovalueForTetrahedral(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 8.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 8.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 8.0f;
    return pow(dx,4)+2*dx*dx*dy*dy+2*dx*dx*dz*dz+pow(dy,4)+2*dy*dy*dz*dz+pow(dz,4)+8*dx*dy*dz-10*dx*dx-10*dy*dy-10*dz*dz+25;
}

float getIsovalueForPillowTooth(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    return pow(dx,4)+pow(dy,4)+pow(dz,4)-dx*dx-dy*dy-dz*dz;
}

float getIsovalueForKampyleOfEudoxus(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    return dy*dy+dz*dz-pow(dx,4)+0.04*dx*dx;
}

float getIsovalueForBarthSextic(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    return 67.77708776*dx*dx*dy*dy*dz*dz-27.41640789*pow(dx,4)*dy*dy-27.41640789*dx*dx*pow(dz,4)+10.47213596*pow(dx,4)*dz*dz-27.41640789*pow(dy,4)*dz*dz+10.47213596*pow(dy,4)*dx*dx+10.47213596*dy*dy*pow(dz,4)-4.236067978*pow(dx,4)-8.472135956*dx*dx*dy*dy-8.472135956*dx*dx*dz*dz+8.472135956*dx*dx-4.236067978*pow(dy,4)-8.472135956*dy*dy*dz*dz+8.472135956*dy*dy-4.236067978*pow(dz,4)+8.472135956*dz*dz-4.236067978;
}

float getIsovalueForA1Singluarity(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    return dx*dx-dy*dy-dz*dz;
}

float getIsovalueForMonoid(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 5.0f;
    return -pow(dx,3)+dy*dy*dz-dy*pow(dz,3);
}

float getIsovalueForStarWarsShip(int x, int y, int z) {
    float dx = (x-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 2.1f;
    float dy = (y-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 2.1f;
    float dz = (z-VOLUME_SIZE/2.0f) / (float)VOLUME_SIZE * 2.1f;
    float dx2 = dx*dx;
    float dy2 = dy*dy;
    float dz2 = dz*dz;
    return (dx2+dy2+dz2<0.2)+((dy2+dz2<0.08)*(dx<0.4)*(dx>0))+(dx2+4*dy2<(1-abs(dz))*0.12)+((abs(dz)<0.95)*(abs(dz)>0.9)*(abs(dx)+abs(dy)*0.3<1))+((abs(dz)<1)*(abs(dz)>0.89))*((abs(dx)<0.7)*(abs(dy)>0.9)+(abs(dy)<0.035)+(dx>dy*0.7-0.05)*(dx<dy*0.7+0.05)+(-dx>dy*0.7-0.05)*(-dx<dy*0.7+0.05)+((abs(dx)+abs(dy)*0.3<1.05)*(abs(dx)+abs(dy)*0.3>0.95)));
}

float getIsovalueForSphere(int x, int y, int z) {
    return pow(x-VOLUME_SIZE/2, 2) + pow(y-VOLUME_SIZE/2, 2) + pow(z-VOLUME_SIZE/2, 2) - pow(SPHERE_RADIUS, 2);
}

typedef float (*IsovalueFunction)(int, int, int);
float getIsovalue(OSN::Noise<3> noise, int x, int y, int z, int shapeIndex) {
    std::vector<IsovalueFunction> isovalueFunctions = {
        getIsovalueForSphere,
        getIsovalueForMonoid,
        getIsovalueForA1Singluarity,
        getIsovalueForQ1,
        getIsovalueForTetrahedral,
        getIsovalueForPillowTooth,
        getIsovalueForKampyleOfEudoxus,
        getIsovalueForBarthSextic,
        getIsovalueForStarWarsShip
    };
    if (shapeIndex % (isovalueFunctions.size()+1) == 0) {
        return getIsovalueFor3DSimplexNoise(noise, x, y, z);
    } else {
        return isovalueFunctions[(shapeIndex % (isovalueFunctions.size()+1)) - 1](x, y, z);
    }
}

#endif /* Isovalues_h */
