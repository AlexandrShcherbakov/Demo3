#ifndef GLUTILITY_H_INCLUDED
#define GLUTILITY_H_INCLUDED

#include "VectorMath.h"
#include "Geometry.h"
#include "SceneGeometry.h"

#include <map>
#include <string>
#include <cstdio>
#include <vector>

#include "GL\glew.h"
#include "GL\freeglut.h"
#include <GL\gl.h>

///Shader programs
//static std::map<std::string, GLuint> progGL;
class ShaderProgram {
private:
    GLuint programId;
    GLuint VAO;
    GLuint pointsVBO;
    GLuint indicesVBO;
    uint elementsCount;
public:
    ShaderProgram(const std::string &shaderName);

    void Draw();
    void setGeometry(Geometry &geom);
    void setGeometry(SceneGeometry &geom);
    void setPoints(const std::vector<vec4> &points);
    void setPoints(const uint size, const float* points);
    void setIndices(const std::vector<uint>& indices);
    void setIndices(const uint size, const uint * indices);
    void setUniform(const std::string& name, mat4 value);
    void setUniform(const std::string& name, vec3 value);
};

mat4 Matrix4Rotate(vec3 v);

static GLuint greenSkeleton;

static std::map<std::string, GLuint> VAO;

static std::map<std::string, GLuint> VBO;


GLuint LoadShader(const std::string& name);

template <class T>
void SetUniform(const std::string& shader, const std::string& uniformName, const T& value);

GLuint CreateGLBuffer(const std::string& vbo, const std::vector<vec4>& points);
GLuint CreateGLBuffer(const std::string& vbo, const int type, const uint size, const void *data);

void CreateBuffersFromGeometry(const std::string& progName, const Geometry& geom);

void DrawGeometry(const std::string& prog, const Geometry& geom);

void Matrix4Perspective(mat4 &m, float fovy, float aspect, float znear, float zfar);

mat4 Matrix4Shift(vec3 shift);

void ThrowExceptionOnGLError(int line, const char *file);

#define CHECK_GL_ERRORS ThrowExceptionOnGLError(__LINE__,__FILE__);

static GLenum g_OpenGLError;
#define OPENGL_CHECK_FOR_ERRORS() \
        if ((g_OpenGLError = glGetError()) != GL_NO_ERROR) { \
			fprintf(stderr, "OpenGL error 0x%X\n", (unsigned)g_OpenGLError); \
                exit(1); \
        }



#endif // GLUTILITY_H_INCLUDED
