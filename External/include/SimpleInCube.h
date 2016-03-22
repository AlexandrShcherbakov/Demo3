#ifndef SIMPLEINCUBE_H
#define SIMPLEINCUBE_H

#include "VectorMath.h"
#include "Polygon.h"

#include <vector>
#include <set>
#include <utility>
#include <queue>
#include <tuple>

std::vector<vec4> SimpleInCube(std::vector<vec4>& triangles, const uint maxDeep);

#endif // SIMPLEINCUBE_H
