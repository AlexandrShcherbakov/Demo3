#include "VectorMath.h"
#include <cmath>
#include <cstdlib>

///Convenient small functions
template <typename T>
inline T sqr(T t) {
    return t * t;
}

///*****************************************************************************
///Classes for vectors and matrices
///*****************************************************************************

///Class for 3D-vectors

///Constructors
vec3::vec3() {}

vec3::vec3(float x, float y, float z): x(x), y(y), z(z) {}

vec3::vec3(float* coords): x(coords[0]), y(coords[1]), z(coords[2]) {}

vec3::vec3(const vec3 &v): x(v.x), y(v.y), z(v.z) {}

///Operators
vec3 vec3::operator+(const vec3 &v) const {
    return vec3(x + v.x, y + v.y, z + v.z);
}

vec3 vec3::operator-(const vec3 &v) const {
    return vec3(x - v.x, y - v.y, z - v.z);
}

vec3 vec3::operator*(const vec3 &v) const {
    return vec3(x * v.x, y * v.y, z * v.z);
}

vec3 vec3::operator*(const float k) const {
    return vec3(x * k, y * k, z * k);
}

vec3 vec3::operator*=(const float k) {
    return (*this) = (*this) * k;
}

vec3 vec3::operator/(const float k) const {
    return vec3(x / k, y / k, z / k);
}

float& vec3::operator[](const uint index) {
    if (index >= 3) throw "Too big index for vec3";
    if (index == 0) return x;
    if (index == 1) return y;
    return z;
}

const float& vec3::operator[](const uint index) const {
    if (index >= 3) throw "Too big index for vec3";
    if (index == 0) return x;
    if (index == 1) return y;
    return z;
}

bool vec3::operator==(const vec3 &v) const {
    return std::abs(x - v.x) < VEC_EPS && std::abs(y - v.y) < VEC_EPS && std::abs(z - v.z) < VEC_EPS;
}

/*bool vec3::operator<(const vec3 &v) const {
    if (x != v.x) return x < v.x;
    if (y != v.y) return y < v.y;
    return z < v.z;
}*/

std::ostream& operator<<(std::ostream& os, const vec3& v)
{
    os << "(" << v.x << "; " << v.y << "; " << v.z << ")";
    return os;
}


///Class for 4D-vectors
///Constructors
vec4::vec4() {}

vec4::vec4(const float x, const float y, const float z, const float w): x(x), y(y), z(z), w(w) {}

vec4::vec4(const float* coords): x(coords[0]), y(coords[1]), z(coords[2]), w(coords[3]) {}

vec4::vec4(const vec4 &v): x(v.x), y(v.y), z(v.z), w(v.w) {}

vec4::vec4(const vec3 &v, const float w): x(v.x), y(v.y), z(v.z), w(w) {}

///Operators
vec4 vec4::operator+(const vec4 &v) const {
    return vec4(x + v.x, y + v.y, z + v.z, w + v.w);
}

vec4 vec4::operator-(const vec4 &v) const {
    return vec4(x - v.x, y - v.y, z - v.z, w - v.w);
}

vec4 vec4::operator*(const vec4 &v) const {
    return vec4(x * v.x, y * v.y, z * v.z, w * v.w);
}

vec4 vec4::operator/(const vec4 &v) const {
    return vec4(x / v.x, y / v.y, z / v.z, w / v.w);
}

vec4 vec4::operator+=(const vec4 &v) {
    return (*this) = (*this) + v;
}

vec4 vec4::operator-=(const vec4 &v) {
    return (*this) = (*this) - v;
}

vec4 vec4::operator*=(const vec4 &v) {
    return (*this) = (*this) * v;
}

vec4 vec4::operator/=(const vec4 &v) {
    return (*this) = (*this) / v;
}

vec4 vec4::operator*(const float k) const {
    return vec4(x * k, y * k, z * k, w * k);
}

vec4 vec4::operator/(const float k) const {
    return vec4(x / k, y / k, z / k, w / k);
}

vec4 vec4::operator/=(const float k) {
    return (*this) = (*this) / k;
}

float& vec4::operator[](const uint index) {
    if (index >= 4) throw "Too big index for vec4";
    if (index == 0) return x;
    if (index == 1) return y;
    if (index == 2) return z;
    return w;
}

const float& vec4::operator[](const uint index) const {
    if (index >= 4) throw "Too big index for vec4";
    if (index == 0) return x;
    if (index == 1) return y;
    if (index == 2) return z;
    return w;
}

bool vec4::operator==(const vec4 &v) const {
    return std::abs(x - v.x) < VEC_EPS && std::abs(y - v.y) < VEC_EPS &&
           std::abs(z - v.z) < VEC_EPS && std::abs(w - v.w) < VEC_EPS;
}

bool vec4::operator!=(const vec4 &v) const {
    return ! ((*this) == v);
}

/*bool vec4::operator<(const vec4 &v) const {
    if (x != v.x) return x < v.x;
    if (y != v.y) return y < v.y;
    if (z != v.z) return z < v.z;
    return w < v.w;
}*/

std::ostream& operator<<(std::ostream& os, const vec4& v) {
    os << "(" << v.x << "; " << v.y << "; " << v.z << "; " << v.w << ")";
    return os;
}

vec3 vec4::xyz() const {
    return vec3(x, y, z);
}


///Class for 4D-matrices

///Constructors
mat4::mat4() {};

mat4::mat4(float diagNum) {
	for (uint i = 0; i < 4; ++i)
		for (uint j = 0; j < 4; ++j)
			rows[i][j] = 0.0f;
    for (uint i = 0; i < 4; ++i)
        rows[i][i] = diagNum;
}

mat4::mat4(float* elements) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            rows[i][j] = elements[i * 4 + j];
}

mat4::mat4(float** elements) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            rows[i][j] = elements[i][j];
}

mat4::mat4(const mat4 &m) {
    for (int i = 0; i < 4; ++i)
        rows[i] = m.rows[i];
}

///Operators
vec4 mat4::operator*(const vec4 &v) const {
    vec4 res;
    for (uint i = 0; i < 4; ++i)
        res[i] = dot(rows[i], v);
    return res;
}

mat4 mat4::operator*(const mat4 &m) const {
    mat4 res;
    for (uint i = 0; i < 4; ++i)
        for (uint j = 0; j < 4; ++j)
            res[i][j] = dot(rows[i], m.col(j));
    return res;
}

mat4 mat4::operator+(const mat4 &m) const{
    mat4 res;
    for (uint i = 0; i < 4; ++i)
        for (uint j = 0; j < 4; ++j)
            res[i][j] = rows[i][j] + m[i][j];
    return res;
}

mat4 mat4::operator+=(const mat4 &m) {
    return (*this) = (*this) + m;
}

mat4 mat4::operator*=(const mat4 &m) {
	return (*this) = (*this) * m;
}

vec4& mat4::operator[](const uint index) {
    if (index >= 4) throw "Too big row index for matrix";
    return rows[index];
}

const vec4& mat4::operator[](const uint index) const{
    if (index >= 4) throw "Too big row index for matrix";
    return rows[index];
}

std::ostream& operator<<(std::ostream& os, const mat4& v) {
    os.precision(8);
    os << std::fixed;
    for (uint i = 0; i < 4; ++i) {
        for (uint j = 0; j < 4; ++j)
            os << "\t" << v[i][j];
        os << "\n";
    }
    return os;
}


///Other methods

///Get column
vec4 mat4::col(const uint index) const {
    vec4 res;
    for (uint i = 0; i < 4; ++i) {
        res[i] = rows[i][index];
        //std::cout << rows[i][index] << ' ' << i << ' ' << index << std::endl;
    }
    return res;
}

mat4 mat4::unmatrixN3() const {
    mat4 m = *this;
    mat4 um(1.0f);
    for (uint i = 0; i < 4; ++i) {
        ///Choose row with non-null i-coord
        uint j;
        for (j = i; j < 4 && m[j][i] == 0; ++j);
        std::swap(m[i], m[j]);
        std::swap(um[i], um[j]);

        ///Substract row from last rows
        for (uint j = i + 1; j < 4; ++j) {
            um[j] -= um[i] / m[i][i] * m[j][i];
            m[j] -= m[i] / m[i][i] * m[j][i];

        }
    }
    float determinant = 1.0;
    for (uint i = 0; i < 4; ++i)
        determinant *= m[i][i];
    //std::cout << determinant << std::endl;
    for (int i = 3; i >= 0; --i) {
        um[i] /= m[i][i];
        m[i] /= m[i][i];
        for (int j = i - 1; j >= 0; --j) {
            um[j] -= um[i] * m[j][i];
            m[j] -= m[i] * m[j][i];
        }
    }

    return um;
}


///*****************************************************************************
///External functions for types
///*****************************************************************************

///Class vec<X>

///Length of vector for similarity to openGL syntax
template <uint X>
float length(const vec<X> &v) {
    return v.length();
}

///Dot production
template <uint X>
float dot(const vec<X> &v, const vec<X> &w) {
    float res = 0.0f;
    for (uint i = 0; i < X; ++i)
        res += v[i] * w[i];
    return res;
}

template <uint X>
vec<X> normalize(const vec<X>& v) {
    float len = length(v);
    if (len == 0.0f) throw "Too big index for vector";
    return v / len;
}

///class vec3
vec3 normalize(const vec3& v) {
    float len = length(v);
    if (len == 0.0f) throw "Null-length vector normalization.";
    return v / len;
}

float length(const vec3 &v) {
    return sqrt(dot(v, v));
}

///Dot production
float dot(const vec3 &v, const vec3 &w) {
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

///Cross production
vec3 cross(const vec3 &v, const vec3 &w) {
    return vec3(v[1] * w[2] - v[2] * w[1],
                v[2] * w[0] - v[0] * w[2],
                v[0] * w[1] - v[1] * w[0]);
}

///Cosine between two vectors
float cos(const vec3 &v, const vec3 &w) {
    return dot(v, w) / length(v) / length(w);
}

///class vec4
vec4 normalize(const vec4& v) {
    float len = length(v);
    if (len < VEC_EPS) throw "Null-length vector normalization.";
    return v / len;
}

float length(const vec4 &v) {
    return sqrt(dot(v, v));
}

///Dot production
float dot(const vec4 &v, const vec4 &w) {
    return v.x * w.x + v.y * w.y + v.z * w.z + v.w * w.w;
}


///Class mat4

///Transpose matrix
mat4 transpose(const mat4& m) {
    mat4 res;
    for (uint i = 0; i < 4; ++i)
        for (uint j = 0; j < 4; ++j)
            res[i][j] = m[j][i];
    return res;
}
