#ifndef DEMO3_VECTORMATH
#define DEMO3_VECTORMATH

#include <cmath>
#include <array>
#include <iostream>

const float VEC_EPS = 1e-5;

///Convenient typedefs for long typenames
typedef unsigned uint;

///Convenient small functions
template <typename T>
T sqr(T t);

///*****************************************************************************
///Classes for vectors and matrices
///*****************************************************************************

///Class for 3D-vectors

class vec3 {
public:
	///Coordinates
    float x, y, z;

    ///Constructors
    vec3();
    vec3(float x, float y, float z);
    vec3(float* coords);
    vec3(const vec3 &v);

    ///Operators
    vec3 operator+(const vec3 &v) const;
    vec3 operator-(const vec3 &v) const;
    vec3 operator*(const vec3 &v) const;
    vec3 operator*(const float k) const;
    vec3 operator/(const float k) const;
    vec3 operator*=(const float k);
    float& operator[](const uint index);
    const float& operator[](const uint index) const;
    bool operator==(const vec3 &v) const;

    friend std::ostream& operator<<(std::ostream& os, const vec3& v);
};

///Class for 4D-vectors
class vec4 {
public:
	///Coordinates
    float x, y, z, w;

    ///Constructors
    vec4();
    vec4(const float x, const float y, const float z, const float w);
    vec4(const float* coords);
    vec4(const vec4 &v);
    vec4(const vec3 &v, const float w);

	///Operators
    vec4 operator+(const vec4 &v) const;
    vec4 operator-(const vec4 &v) const;
    vec4 operator*(const vec4 &v) const;
    vec4 operator/(const vec4 &v) const;
    vec4 operator+=(const vec4 &v);
    vec4 operator-=(const vec4 &v);
    vec4 operator*=(const vec4 &v);
    vec4 operator/=(const vec4 &v);
    vec4 operator*(const float k) const;
    vec4 operator/(const float k) const;
    vec4 operator/=(const float k);
    float& operator[](const uint index);
    const float& operator[](const uint index) const;
    bool operator==(const vec4 &v) const;
    bool operator!=(const vec4 &v) const;

    friend std::ostream& operator<<(std::ostream& os, const vec4& v);

    vec3 xyz() const;
};

class mat4 {
protected:
    ///Rows
    vec4 rows[4];
public:

    ///Constructors
    mat4();
    mat4(float diagNum);
    mat4(float* elements);
    mat4(float** elements);
    mat4(const mat4 &m);

    ///Operators
    vec4 operator*(const vec4 &v) const;
    mat4 operator*(const mat4 &m) const;
    mat4 operator+(const mat4 &m) const;
    mat4 operator+=(const mat4 &m);
    mat4 operator*=(const mat4 &m);
    vec4& operator[](const uint index);
    const vec4& operator[](const uint index) const;

    friend std::ostream& operator<<(std::ostream& os, const mat4& v);

    ///Other methods
    vec4 col(const uint index) const;
    mat4 unmatrixN3() const;
};


///*****************************************************************************
///External functions for types
///*****************************************************************************

///Class vec<X>

///class vec3

///Length of vector for similarity to openGL syntax
float length(const vec3 &v);

vec3 normalize(const vec3 &v);

///Dot production
float dot(const vec3 &v, const vec3 &w);

///Cross production
vec3 cross(const vec3 &v, const vec3 &w);

///Cosine between two vectors
float cos(const vec3 &v, const vec3 &w);

///class vec4

///Length of vector for similarity to openGL syntax
float length(const vec4 &v);

vec4 normalize(const vec4& v);

///Dot production
float dot(const vec4 &v, const vec4 &w);


///Class mat4

///Transpose matrix
mat4 transpose(const mat4& m);

#endif
