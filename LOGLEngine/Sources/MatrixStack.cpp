// Local Headers
#include "MatrixStack.h"

using namespace std;

MatrixStack::MatrixStack() {

   mstack = make_shared< stack<mat4> >();
   mstack->push(mat4(1.0));
}

MatrixStack::~MatrixStack() { }

void MatrixStack::pushMatrix()
{
	const mat4 &top = mstack->top();
   mstack->push(top);
   assert(mstack->size() < 100);
}

void MatrixStack::popMatrix()
{
   assert(!mstack->empty());
   mstack->pop();
   // There should always be one matrix left.
   assert(!mstack->empty());
}

void MatrixStack::reset()
{
    while (mstack->size() > 1) {
        mstack->pop();
    }
    loadIdentity();
}

void MatrixStack::loadIdentity()
{
   mat4 &top = mstack->top();
   top = glm::mat4(1);
}

void MatrixStack::translate(const vec3 &offset)
{
   mat4 &top = mstack->top();
	mat4 t = glm::translate(mat4(1.0), offset);
	top *= t;
}

void MatrixStack::scale(const vec3 &scaleV)
{
   mat4 &top = mstack->top();
	mat4 s = glm::scale(mat4(1.0), scaleV);
	top *= s;
}

void MatrixStack::scale(float size)
{
   mat4 &top = mstack->top();
	mat4 s = glm::scale(mat4(1.0), vec3(size));
	top *= s;
}

void MatrixStack::rotate(float angle, const vec3 &axis)
{
   mat4 &top = mstack->top();
	mat4 r = glm::rotate(mat4(1.0), angle, axis);
	top *= r;
}

void MatrixStack::multMatrix(const glm::mat4 &matrix)
{
   mat4 &top = mstack->top();
	top *= matrix;
}

void MatrixStack::ortho(float left, float right, float bottom, float top, float zNear, float zFar)
{
   assert(left != right);
   assert(bottom != top);
   assert(zFar != zNear);
   mat4 &ctm = mstack->top();
   ctm *= glm::ortho(left, right, bottom, top, zNear, zFar);
}


void MatrixStack::ortho2D(float left, float right, float bottom, float top)
{
	assert(left != right);
	assert(bottom != top);
	mat4 &ctm = mstack->top();
	ctm *= glm::ortho(left, right, bottom, top);
}


void MatrixStack::frustum(float left, float right, float bottom, float top, float zNear, float zFar)
{
   mat4 &ctm = mstack->top();
	ctm *= glm::frustum(left, right, bottom, top, zNear, zFar);
}

void MatrixStack::lookAt(vec3 eye, vec3 target, vec3 up)
{
   mat4 &top = mstack->top();
	top *= glm::lookAt(eye, target, up);
}

void MatrixStack::perspective(float fovy, float aspect, float zNear, float zFar)
{
   mat4 &top = mstack->top();
	top *= glm::perspective(fovy, aspect, zNear, zFar);
}

const glm::mat4 &MatrixStack::topMatrix()
{
	return mstack->top();
}

void MatrixStack::print(const mat4 &mat, const char *name) const
{
   if(name) {
      printf("%s = [\n", name);
   }
   for(int i = 0; i < 4; ++i) {
      for(int j = 0; j < 4; ++j) {
         printf("%- 5.2f ", mat[i][j]);
      }
      printf("\n");
   }
   if(name) {
      printf("];");
   }
   printf("\n");
}

void MatrixStack::print(const char *name) const
{
   print(mstack->top(), name);
}

