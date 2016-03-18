#include "GLUtility.h"

void ThrowExceptionOnGLError(int line, const char *file) {
  static char errMsg[512];

  GLenum gl_error = glGetError();

  if(gl_error == GL_NO_ERROR)
    return;

  switch(gl_error)
  {
  case GL_INVALID_ENUM:
    sprintf(errMsg, "GL_INVALID_ENUM file %s line %d\n", file, line);
    break;

  case GL_INVALID_VALUE:
    sprintf(errMsg, "GL_INVALID_VALUE file %s line %d\n",  file, line);
    break;

  case GL_INVALID_OPERATION:
    sprintf(errMsg, "GL_INVALID_OPERATION file %s line %d\n",  file, line);
    break;

  case GL_STACK_OVERFLOW:
    sprintf(errMsg, "GL_STACK_OVERFLOW file %s line %d\n",  file, line);
    break;

  case GL_STACK_UNDERFLOW:
    sprintf(errMsg, "GL_STACK_UNDERFLOW file %s line %d\n",  file, line);
    break;

  case GL_OUT_OF_MEMORY:
    sprintf(errMsg, "GL_OUT_OF_MEMORY file %s line %d\n",  file, line);
    break;

  /*case GL_TABLE_TOO_LARGE:
    sprintf(errMsg, "GL_TABLE_TOO_LARGE file %s line %d\n",  file, line);
    break;*/

  case GL_NO_ERROR:
    break;

  default:
    sprintf(errMsg, "Unknown error @ file %s line %d\n",  file, line);
    break;
  }

  if(gl_error != GL_NO_ERROR)
    fprintf(stderr, "!!!ERROR BUGURT\n%s", errMsg);
}


GLint ShaderStatus(GLuint shader, GLenum param) {
	GLint status, length;
	GLchar buffer[1024];

	glGetShaderiv(shader, param, &status);

	if (status != GL_TRUE) {
		glGetShaderInfoLog(shader, 1024, &length, buffer);
		fprintf(stderr, "Shader: %s\n", (const char*)buffer);
	}
	OPENGL_CHECK_FOR_ERRORS();
	return status;
}

int LoadSource(const char *shaderName, char **textOut, int *textLen) {
	FILE *input;

	input = fopen(shaderName, "r");

	fseek(input, 0, SEEK_END);
	*textLen = ftell(input);
	rewind(input);

	*textOut = (char*)calloc(sizeof(*textOut), *textLen + 1);

	*textLen = fread(*textOut, sizeof(**textOut), *textLen, input);
	fclose(input);
	return 1;
}

void CompileShader(const std::string &name, GLuint &shader, GLenum shaderType) {
 	//Initialization of shader
	shader = glCreateShader(shaderType);                                         CHECK_GL_ERRORS

	//Variables for shader texts
    char *shaderSource;
    int sourceLength;

    //Read vertex shader
    LoadSource(name.c_str(), &shaderSource, &sourceLength);

    //Compile vertex shader
    glShaderSource(shader, 1, (const GLchar**)&shaderSource,
                (const GLint*)&sourceLength);                                    CHECK_GL_ERRORS
    glCompileShader(shader);                                                     CHECK_GL_ERRORS

    //Release source memory
    free(shaderSource);

    //Check whether the compilation is successful
    if (ShaderStatus(shader, GL_COMPILE_STATUS) != GL_TRUE) {
        fprintf(stderr, "BUGURT!!! Line: %d\n", __LINE__);
        exit(1);
    }
}

GLint ShaderProgramStatus(GLuint program, GLenum param)
{
	GLint status, length;
	GLchar buffer[1024];
	glGetProgramiv(program, param, &status);

	if (status != GL_TRUE)
	{
		glGetProgramInfoLog(program, 1024, &length, buffer);
		fprintf(stderr, "Shader program: %s\n", (const char*)buffer);
	}

	OPENGL_CHECK_FOR_ERRORS();
	return status;
}

ShaderProgram::ShaderProgram(const std::string &shaderName) {
    GLuint vertSh, fragSh;
    CompileShader("shaders/" + shaderName + ".vert", vertSh, GL_VERTEX_SHADER);
    CompileShader("shaders/" + shaderName + ".frag", fragSh, GL_FRAGMENT_SHADER);

    programId = glCreateProgram();                                               CHECK_GL_ERRORS
    glAttachShader(programId, vertSh);                                           CHECK_GL_ERRORS
    glAttachShader(programId, fragSh);                                           CHECK_GL_ERRORS

    glLinkProgram(programId);                                                    CHECK_GL_ERRORS

    if (ShaderProgramStatus(programId, GL_LINK_STATUS) != GL_TRUE) {
        fprintf(stderr, "BUGURT!!! Line: %d\n", __LINE__);
        exit(1);
    }
    glGenVertexArrays(1, &VAO);                                                  CHECK_GL_ERRORS
    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    glUseProgram(0);                                                             CHECK_GL_ERRORS
}

void ShaderProgram::setPoints(const uint size, const float *points) {
    glUseProgram(programId);                                                     CHECK_GL_ERRORS
    glBindVertexArray(VAO);                                                      CHECK_GL_ERRORS

    glDeleteBuffers(1, &pointsVBO);                                              CHECK_GL_ERRORS
    glGenBuffers(1, &pointsVBO);                                                 CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);                                    CHECK_GL_ERRORS
    glBufferData(GL_ARRAY_BUFFER, size * sizeof(*points), points, GL_STATIC_DRAW); CHECK_GL_ERRORS

    auto pointsLocation = glGetAttribLocation(programId, "points");              CHECK_GL_ERRORS
    glEnableVertexAttribArray(pointsLocation);                                   CHECK_GL_ERRORS
    glVertexAttribPointer(pointsLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);          CHECK_GL_ERRORS

    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, 0);                                            CHECK_GL_ERRORS
}

void ShaderProgram::setPoints(const std::vector<vec4>& points) {
    std::vector<float> data;
    for (auto i: points)
        for (uint j = 0; j < 4; ++j)
			data.push_back(i[j]);
    setPoints(data.size(), data.data());
}

void ShaderProgram::setIndices(const uint size, const uint * indices) {
    glUseProgram(programId);                                                     CHECK_GL_ERRORS
    glBindVertexArray(VAO);                                                      CHECK_GL_ERRORS

    glDeleteBuffers(1, &indicesVBO);                                             CHECK_GL_ERRORS
    glGenBuffers(1, &indicesVBO);                                                CHECK_GL_ERRORS
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indicesVBO);                           CHECK_GL_ERRORS
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(*indices), indices, GL_STATIC_DRAW); CHECK_GL_ERRORS

    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);                                    CHECK_GL_ERRORS
    elementsCount = size;
}

void ShaderProgram::setIndices(const std::vector<uint>& indices) {
    setIndices(indices.size(), indices.data());
}

void ShaderProgram::setGeometry(Geometry &geom) {
    setPoints(geom.getPoints());
    setIndices(geom.getEdgeListIndices());
}

void ShaderProgram::setGeometry(SceneGeometry &geom) {
    setPoints(geom.getPoints());
    setIndices(geom.getEdges());
}

void ShaderProgram::Draw() {
    glUseProgram(programId);                                                     CHECK_GL_ERRORS
    glBindVertexArray(VAO);                                                      CHECK_GL_ERRORS
    glPointSize(3);                                                              CHECK_GL_ERRORS
    glDrawElements(GL_POINTS, elementsCount, GL_UNSIGNED_INT, 0);                CHECK_GL_ERRORS
    glDrawElements(GL_LINES, elementsCount, GL_UNSIGNED_INT, 0);                CHECK_GL_ERRORS
    glUseProgram(0);                                                             CHECK_GL_ERRORS
    glBindVertexArray(0);                                                        CHECK_GL_ERRORS
}

void ShaderProgram::setUniform(const std::string& name, mat4 value) {
    glUseProgram(programId);                                                     CHECK_GL_ERRORS
	GLint location = glGetUniformLocation(programId, name.c_str());              CHECK_GL_ERRORS
    if (location != -1) {
        float M[16];
		for (uint i = 0; i < 4; ++i)
			for (uint j = 0; j < 4; ++j)
			M[4 * i + j] = value[i][j];
		glUniformMatrix4fv(location, 1, GL_TRUE, M);                             CHECK_GL_ERRORS
    } else {
        fprintf(stderr, "Location: %s not found!\n", name.c_str());
    }
    glUseProgram(0);                                                             CHECK_GL_ERRORS
}

void ShaderProgram::setUniform(const std::string& name, vec3 value) {
    glUseProgram(programId);                                                     CHECK_GL_ERRORS
    GLint location = glGetUniformLocation(programId, name.c_str());              CHECK_GL_ERRORS
    if (location != -1) {
		glUniform3f(location, value.x, value.y, value.z);                        CHECK_GL_ERRORS
    } else {
        fprintf(stderr, "Location: %s not found!\n", name.c_str());
    }
    glUseProgram(0);                                                             CHECK_GL_ERRORS
}

GLuint LoadShader(const std::string& name) {
	GLuint vertSh, fragSh;
    CompileShader("shaders/" + name + ".vert", vertSh, GL_VERTEX_SHADER);
    CompileShader("shaders/" + name + ".frag", fragSh, GL_FRAGMENT_SHADER);
    GLuint greenSkeleton;
    greenSkeleton = glCreateProgram();                                            CHECK_GL_ERRORS
    glAttachShader(greenSkeleton, vertSh);                                        CHECK_GL_ERRORS
    glAttachShader(greenSkeleton, fragSh);                                        CHECK_GL_ERRORS

    glLinkProgram(greenSkeleton);                                                 CHECK_GL_ERRORS

    if (ShaderProgramStatus(greenSkeleton, GL_LINK_STATUS) != GL_TRUE) {
        fprintf(stderr, "BUGURT!!! Line: %d\n", __LINE__);
        exit(1);
    }
    return greenSkeleton;
}

void SetUniformToLocation(const GLint& location, const mat4& value) {
}


void SetUniformToLocation(const GLint& location, const int& value) {
    glUniform1i(location, value);                                                CHECK_GL_ERRORS
}

template <>
void SetUniform<mat4>(const std::string& shader, const std::string& uniformName, const mat4& value) {
	glUseProgram(greenSkeleton);
	GLint location = glGetUniformLocation(greenSkeleton, uniformName.c_str());  CHECK_GL_ERRORS
    if (location != -1) {
        float M[16];
		for (uint i = 0; i < 4; ++i)
			for (uint j = 0; j < 4; ++j)
			M[4 * i + j] = value[i][j];
		glUniformMatrix4fv(location, 1, GL_TRUE, M);                                 CHECK_GL_ERRORS
    } else {
        fprintf(stderr, "Location: %s not found!\n", uniformName.c_str());
    }
    glUseProgram(0);
}


template <class T>
void SetUniform(const std::string& shader, const std::string& uniformName, const T& value) {
	glUseProgram(greenSkeleton);
	GLint location = glGetUniformLocation(greenSkeleton, uniformName.c_str());  CHECK_GL_ERRORS
    if (location != -1) {
        SetUniformToLocation(location, value);
    } else {
        fprintf(stderr, "Location: %s not found!\n", uniformName.c_str());
    }
    glUseProgram(0);
}

GLuint CreateGLBuffer(const std::string& vbo, const int type, const uint size, const void *data) {
    GLuint vb;
	glDeleteBuffers(1, &vb);                                               CHECK_GL_ERRORS
    glGenBuffers(1, &vb);                                                  CHECK_GL_ERRORS
    glBindBuffer(type, vb);                                                CHECK_GL_ERRORS
	glBufferData(type, size, data, GL_STATIC_DRAW);                              CHECK_GL_ERRORS
    glBindBuffer(type, 0);                                                       CHECK_GL_ERRORS
    return vb;
}

GLuint CreateGLBuffer(const std::string& vbo, const std::vector<vec4>& points) {
    GLuint vb;
    glGenBuffers(1, &vb);                                                  CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, vb);                                     CHECK_GL_ERRORS
    std::vector<float> data;
    for (auto i: points)
        for (uint j = 0; j < 4; ++j)
			data.push_back(i[j]);
	glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(data[0]), data.data(), GL_STATIC_DRAW); CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, 0);                                            CHECK_GL_ERRORS
    return vb;
}


void CreateBuffersFromGeometry(const std::string& progName, const Geometry& geom) {
    CreateGLBuffer("points", geom.getPoints());
    CreateGLBuffer("indices", GL_ELEMENT_ARRAY_BUFFER, geom.getTriangles().size() * sizeof(geom.getTriangles()[0]), geom.getTriangles().data());

    glGenVertexArrays(1, &VAO[progName]);                                        CHECK_GL_ERRORS
    glBindVertexArray(VAO[progName]);                                            CHECK_GL_ERRORS
    auto pointsLocation = glGetAttribLocation(greenSkeleton, "points");       CHECK_GL_ERRORS
    glBindBuffer(GL_ARRAY_BUFFER, VBO["points"]);                                CHECK_GL_ERRORS
    glEnableVertexAttribArray(pointsLocation);                                   CHECK_GL_ERRORS
    glVertexAttribPointer(pointsLocation, 4, GL_FLOAT, GL_FALSE, 0, 0);          CHECK_GL_ERRORS

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBO["indices"]);                       CHECK_GL_ERRORS
    glBindVertexArray(VAO[progName]);                                            CHECK_GL_ERRORS
}


void DrawGeometry(const std::string& prog, const Geometry& geom) {
    glBindVertexArray(VAO[prog]);                                                CHECK_GL_ERRORS
    glUseProgram(greenSkeleton);                                                     CHECK_GL_ERRORS
    glDrawElements(greenSkeleton, geom.getTriangles().size(), GL_UNSIGNED_INT, 0); CHECK_GL_ERRORS
}

void BindBufferToVAO(const std::string& vao, const std::string& vbo, const int bufType, const uint index, const uint size, const int type) {
	glBindVertexArray(VAO[vao]);                                                 CHECK_GL_ERRORS
	glBindBuffer(bufType, VBO[vbo]);                                             CHECK_GL_ERRORS
	glVertexAttribPointer(index, size, type, GL_FALSE, 0, 0);                    CHECK_GL_ERRORS
	glBindBuffer(bufType, 0);                                                    CHECK_GL_ERRORS
	glBindVertexArray(0);                                                        CHECK_GL_ERRORS
}


void Matrix4Perspective(mat4 &M, float fovy, float aspect, float znear, float zfar) {
	//Convert fovy from graduses to radians
	float f = 1 / tanf(fovy * 3.1415926535f / 360.0f),
			A = (zfar + znear) / (znear - zfar),
			B = (2 * zfar * znear) / (znear - zfar);

	M[0][0] = f / aspect; M[0][1] =  0; M[0][2] =  0; M[0][3] =  0;
	M[1][0] = 0;          M[1][1] =  f; M[1][2] =  0; M[1][3] =  0;
	M[2][0] = 0;          M[2][1] =  0; M[2][2] =  A; M[2][3] =  B;
	M[3][0] = 0;          M[3][1] =  0; M[3][2] = -1; M[3][3] =  0;
}

mat4 Matrix4Shift(vec3 shift) {
	mat4 M(1.0f);
    M[0][3] = shift.x;
    M[1][3] = shift.y;
    M[2][3] = shift.z;
    return M;
}

mat4 GenRotateMatrix(vec4 v, float a) {
    mat4 M(1.0f);
    float ca = cos(a), sa = sin(a);
    float nca = 1 - ca;
    M[0][0] =       ca + nca * v.x * v.x; M[0][1] = nca * v.x * v.y - sa * v.z; M[0][2] = nca * v.x * v.y + sa * v.y;
    M[1][0] = nca * v.y * v.x + sa * v.z; M[1][1] =       ca + nca * v.y * v.y; M[1][2] = nca * v.y * v.z - sa * v.x;
    M[2][0] = nca * v.z * v.x - sa * v.y; M[2][1] = nca * v.z * v.y + sa * v.x; M[2][2] = ca + nca * v.z * v.z;
    return M;
}

mat4 Matrix4Rotate(vec3 v) {
    vec4 y(0, 1, 0, 1);
    vec4 z(0, 0, 1, 1);
    mat4 Mx = GenRotateMatrix(vec4(1, 0, 0, 1), v.x);
    y = Mx * y;
    mat4 My = GenRotateMatrix(y, v.y);
    mat4 Mxy = My * Mx;
    z = Mxy * z;
    mat4 Mz = GenRotateMatrix(z, v.z);
    return Mz * Mxy;
}
