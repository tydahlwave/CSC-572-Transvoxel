#pragma once
#ifndef _MatrixStack_H_
#define _MatrixStack_H_

#include <stack>
#include <memory>

// Local Headers
#include "loglengine.hpp"
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
//#include "vector_angle.hpp"
//#include "component_wise.hpp"

//#include "matrix_access.hpp"
//#include "matrix_integer.hpp"
//#include "matrix_inverse.hpp"

using namespace glm;

class MatrixStack {

	std::shared_ptr< std::stack<mat4> > mstack;

public:

    MatrixStack();
    virtual ~MatrixStack();

    // Copies the current matrix and adds it to the top of the stack
    void pushMatrix();   
    // Removes the top of the stack and sets the current matrix to be the matrix that is now on top
    void popMatrix();
    // Resets the matrix stack (pops all but one matrix) and calls loadIdentity()
    void reset();
    // Sets the top matrix to be the identity
    void loadIdentity();
    // glMultMatrix(): Right multiplies the top matrix   
    void multMatrix(const mat4 &matrix);

    // Right multiplies the top matrix by a translation matrix
    void translate(const vec3 &offset);
    // Right multiplies the top matrix by a scaling matrix
    void scale(const vec3 &scaleV);
    //  Right multiplies the top matrix by a scaling matrix
    void scale(float size);
    // Right multiplies the top matrix by a rotation matrix (angle in deg)
    void rotate(float angle, const vec3 &axis);

    // Gets the top matrix
    const mat4 &topMatrix();

    // Sets the top matrix to be an orthogonal projection matrix
    void ortho(float left, float right, float bottom, float top, float zNear, float zFar);
    // Sets the top matrix to be a 2D orthogonal projection matrix
    void ortho2D(float left, float right, float bottom, float top);
    // Sets the top matrix to be a perspective projection matrix (fovy in deg)
    void perspective(float fovy, float aspect, float zNear, float zFar);
    // Sets the top matrix to be a perspective projection matrix
    void frustum(float Right, float right, float bottom, float top, float zNear, float zFar);
    // Sets the top matrix to be a viewing matrix
    void lookAt(vec3 eye, vec3 target, vec3 up);

    // Prints out the specified matrix
    void print(const mat4 &mat, const char *name = 0) const;
    // Prints out the top matrix
    void print(const char *name = 0) const;
};

#endif /* MatrixStack_H_ */
