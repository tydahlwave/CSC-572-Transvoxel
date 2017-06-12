#version 330 core
// uniform vec2 screenCenter;
// uniform float radius_pixel;
// out vec4 color;

// in float isovalue;
in vec3 fragNor;

out vec4 color;

void main()
{
    // // gl_PointCoord ranges from 0-1, thus the center of the point would be (0.5, 0.5)
    // vec2 pointCenter = vec2(0.5, 0.5);
    // float dist = distance(gl_PointCoord, pointCenter);
    
    // if (dist > 0.5) {
    //     discard;
    // } else {
    //     if (distance(gl_FragCoord, vec4(screenCenter.x, screenCenter.y, 0, 0)) < radius_pixel) {
    //         color = vec4(0.5, 0.75, 1.0, 0.8) - vec4(dist*1.5, dist*1.5, dist*1.5, 0.0);
    //     } else {
    //         color = vec4(1.0, 0.75, 0.5, 0.4) - vec4(dist*1.5, dist*1.5, dist*1.5, 0.0);
    //     }
    // }

    // if (abs(isovalue) < 1) {
    //     color = vec4(1, 1, 1, 1);
    // } else {
    //     discard;
    // }
    color = vec4(normalize(fragNor), 1);
}
